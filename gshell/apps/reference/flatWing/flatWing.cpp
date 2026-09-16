//
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//
//                 Computational Structural Mechanics Lab
//                         University of Cambridge
//
// This software is copyrighted and all rights are retained by the CSMLab. All
// use, disclosure, and/or reproduction of any part not expressly authorized by
// F Cirak is prohibited. (C) 2010.
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//

// includes from stl
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <functional>
#include <ctime>

// includes from boost
#include <boost/timer/timer.hpp>

// includes from corlib
#include <corlib/PropertiesParser.hpp>
#include <corlib/verify.hpp>
#include <corlib/Accessor.hpp>
#include <corlib/VTUwriter.hpp>
#include <corlib/fuzzyEqual.hpp>
#include <corlib/Quadrature.hpp>
#include <corlib/Integrator.hpp>
#include <corlib/DofNumberer.hpp>
#include <corlib/Distributor.hpp>
#include <corlib/MatrixAssembler.hpp>
#include <corlib/VectorAssembler.hpp>
#include <corlib/BodyForce.hpp>
#include <corlib/SetData.hpp>
#include <corlib/Constraints.hpp>
#include <corlib/AccumulateQuantity.hpp>
#include <corlib/MeshSizeInfo.hpp>

#if   defined(DYNAMICS_BACKWARDEULER)
#include <corlib/BackwardEulerAssembler.hpp>
#elif defined(DYNAMICS_NEWMARK)
#include <corlib/NewmarkAssembler.hpp>
#endif

#if   defined SUPERLU
#include <corlib/SystemSolveLU.hpp>
#elif defined PARDISO_MKL
#include <corlib/SystemSolvePardiso.hpp>
#else
#include <corlib/SystemSolveEigenSparse.hpp>
#endif

// includes from subdiv::surf
#include <subdiv/surf/Vertex.hpp>
#include <subdiv/surf/Edge.hpp>
#include <subdiv/surf/Facet.hpp>
#include <subdiv/surf/FacetTree.hpp>
#include <subdiv/surf/Subdivision.hpp>
#include <subdiv/surf/Mesh.hpp>
#include <subdiv/surf/ShapeFunSubdivision.hpp>

// includes from solid::fem
#include <solid/fem/NodalForces.hpp>

// includes from solid::material
#include <solid/material/MaterialContainer.hpp>
#include <solid/material/svenant/SVenant.hpp>

// includes from gshell::fem
#include <gshell/fem/NodeShell.hpp>
#include <gshell/fem/NodeDynamic.hpp>
#include <gshell/fem/ShapeFunCache.hpp>
#include <gshell/fem/ElementShell.hpp>
#define GSHELL_KL_SYMMETRICSTIFFNESS
#include <gshell/fem/ElementKLDynamic.hpp>
#ifdef PARDISO_MKL
#include <gshell/fem/miscThreaded.hpp>
#endif
#include <gshell/fem/MeshShell.hpp>
#include <gshell/fem/Monitor.hpp>

// local includes
#include "Parameters.hpp"
#include "DirichletFun.hpp"
#include "RotationFlapping.hpp"
#include "Load.hpp"
#include "Write.hpp"

//==============================================================================
// everything starts here
int main( int argc, char * args[] ){

    // Eigen
    namespace eigenX = corlib::eigenX;
    // timer
    boost::timer::cpu_timer timer;
;

    //--------------------------------------------------------------------------
#if LOOP_
    const subdiv::surf::method    subdivMethod   = subdiv::surf::LOOP;
#elif CATMULL_CLARK_
    const subdiv::surf::method    subdivMethod   = subdiv::surf::CATMULL_CLARK;
#endif
    // simplex shape
    const corlib::shape myShape = subdiv::surf::SubdivisionShape< subdivMethod >::myShape;

    typedef subdiv::surf::Vertex                                       Vertex;
    typedef subdiv::surf::Facet< myShape, Vertex >                     Facet;
    typedef subdiv::surf::Edge< Vertex, Facet >                        Edge;
    typedef subdiv::surf::FacetTree< myShape, Vertex >                 FacetTree;
    typedef subdiv::surf::Mesh< Vertex, Edge, Facet, FacetTree >       Mesh;

    typedef subdiv::surf::Subdivision< subdivMethod, FacetTree >       Subdivision;

    typedef subdiv::surf::ShapeFunSubdivision< Facet, Subdivision >    ShapeFun;

    // read input data using Properties Parser
    const std::string inputData = ( argc > 1 ) ? args[ 1 ] : "input.dat";
    std::cout << "Using input file" << std::endl
              << "    \"" << inputData << "\"" << std::endl;
    std::ifstream inp( inputData.c_str() );
    app::Parameters parameters( inp );
    inp.close();

    //read input mesh and tags
    std::ifstream inputStream( parameters.inputMeshFileName.c_str() );
    FTL_VERIFY_DESCRIPTIVE( inputStream.is_open(),
                            "Could not open file %s\n",
                            parameters.inputMeshFileName.c_str() );
    std::ifstream secondaryInputStream( parameters.secondaryInputFileName.c_str() );
    FTL_VERIFY_DESCRIPTIVE( secondaryInputStream.is_open(),
                            "Could not open file %s\n",
                            parameters.secondaryInputFileName.c_str() );

    // Surface
    Mesh sMesh( inputStream );
    sMesh.setVerbose( );
    sMesh.readTags( secondaryInputStream ); // prior to buildFacetTopo
    sMesh.buildFacetTopology( );  // OK

    // Subdivision scheme
    Subdivision * subdiv = new Subdivision();

    //--------------------------------------------------------------------------
    // create g-shell mesh
    //const enum corlib::shape gShape = myShape
    const unsigned dof = 3;  // Kirchhoff-Love kinematics
    const corlib::shape gShape = myShape;

    typedef gshell::fem::NodeShell                                     GBasisNode;
    typedef gshell::fem::NodeDynamic< GBasisNode, dof >                GNode;
    typedef gshell::fem::ShapeFunCache< gShape >                       GShapeFunCache;
    typedef gshell::fem::ElementShell< GNode, GShapeFunCache >         GBasisElement;
    typedef solid::material::SVenant                                   GMaterial;
    typedef gshell::fem::ElementKLDynamic< GBasisElement, GMaterial >  GElement;
    typedef gshell::fem::MeshShell< GElement, ShapeFun >               GMesh;

#ifdef LOOP_
    typedef corlib::Quadrature< myShape, 3 >                           GQuad;
#elif CATMULL_CLARK_
    typedef corlib::Quadrature< myShape, 4 >                           GQuad;
#endif
    typedef GNode::VecDim                                              VecDim;
    typedef GNode::VecDof                                              VecDof;
    typedef GElement::VecLDim                                          VecLDim;

    typedef boost::function< void( GElement*, const VecLDim &, 
                                   const double & ) >                  GIntegrand;
    typedef boost::function< void( const GElement *, const VecLDim &, 
                                   const double &,
                                   Eigen::MatrixXd & ) >               GIntegrand2;
    typedef corlib::ComputeElementMatrix< GQuad, GIntegrand2 >         GMatrixComputer;

    typedef corlib::AccumulateQuantity< GElement, double >             GESumDouble;
    typedef corlib::AccumulateQuantity< GElement, unsigned,
                                        corlib::Max >                  GEMaxUnsigned;

#if   defined SUPERLU
    typedef corlib::SystemSolveLU                                      GSysSolve;
#elif defined PARDISO_MKL
    typedef corlib::SystemSolvePardiso                                 GSysSolve;
#else
    typedef corlib::SystemSolveEigenSparse                             GSysSolve;
#endif

    // the quadrature
    GQuad quadrature;

    // the g-mesh
    GMesh gMesh( sMesh, quadrature );

    // create material container
    solid::material::MaterialContainer * mc = solid::material::MaterialContainer::instance();
    // read material
    std::ifstream mf( parameters.materialFileName.c_str() );
    FTL_VERIFY_DESCRIPTIVE( mf.is_open(),
                            "Cannot open file %s\n", parameters.materialFileName.c_str() );
    mc->readMaterialStream( mf );
    mf.close();

    // set materials
    GMaterial * material = dynamic_cast< GMaterial * >( mc->getMaterial( 0 ) );

    // cache material to elements
    gMesh.iterateOverElements( std::bind( &GElement::cacheMaterial, 
                                          std::placeholders::_1, material ) );

    // set thickness
    const double thickness = parameters.thickness;
    gMesh.iterateOverElements( std::bind( &GElement::setThickness, 
                                          std::placeholders::_1, thickness ) );

    // set unique tangents
    gMesh.iterateOverElements( std::bind( &GElement::setLimitCoefficientsAtNodes, 
                                          std::placeholders::_1 ) );
    gMesh.iterateOverElements( std::bind( &GElement::setUniqueTangentCoefficientsAtNodes, 
                                          std::placeholders::_1 ) );

    // a word to the user
#if   defined(DYNAMICS_BACKWARDEULER)
    std::cout << "Running backward Euler" << std::endl;
#elif defined(DYNAMICS_NEWMARK)
    std::cout << "Running Newmark" << std::endl;
    const double newmarkBeta = parameters.newmarkBeta;
    const double newmarkGamma = parameters.newmarkGamma;
#else
    FTL_VERIFY_DESCRIPTIVE( false, "Don't know what to do\n" );
#endif

    // concentrated loads
    solid::fem::NodalForces< GNode > gNodalForce;  // empty
    // add forces from file
    if ( parameters.nodalForcesFileName != "" ) {
        std::ifstream nodalForcesStream( parameters.nodalForcesFileName.c_str() );
        FTL_VERIFY_DESCRIPTIVE( nodalForcesStream.is_open(),
                                "Could not open file %s\n",
                                parameters.nodalForcesFileName.c_str() );
        gNodalForce.add( nodalForcesStream );
    }

    // number dofs
    const unsigned numdofs = gMesh.iterateOverNodes( corlib::dofNumberFun( &GNode::numberDOFs, 0 ) );
    std::cout << numdofs << " degrees of freedom" << std::endl;

    // set up an animation writer
    std::cout << "Writing VTU file" << std::endl;
    const std::string makeAnimationFolder = "mkdir -p " + parameters.animationFolder;
    FTL_VERIFY( !system( makeAnimationFolder.c_str() ) );
    corlib::VTUanim animation( parameters.animationBaseFileName, "vtu" );

    // set-up history file
    const std::string histFileName = parameters.animationBaseFileName + ".history";
    std::ofstream histFile( histFileName.c_str() );
    FTL_VERIFY_DESCRIPTIVE( histFile.is_open(),
                            "Failed to open file %s\n", histFileName.c_str() );

    // allocate linear system solver
    GSysSolve * sysmat = new GSysSolve( numdofs );


    // body force
    typedef boost::function< VecDim( VecDim ) >                 FunType2;
    typedef corlib::BodyForce< GElement, FunType2 >             GBodyForce;
    GBodyForce bodyForce( std::bind( &app::pressure, std::placeholders::_1, thickness ),
                          &GElement::bodyForce< FunType2 > );

    // mesh size info
    corlib::MeshSizeInfo< GElement > gMeshSizeInfo;
    gMeshSizeInfo = gMesh.iterateOverElements( gMeshSizeInfo );
    gMeshSizeInfo.writeStats( std::cout );

    //==========================================================================
    // load control parameters (load or time loop control)
    const unsigned timeStepMax = parameters.timeStepMax;
    const double period = 1./150.;  // stroke period
    const double timeMax = period;  //0.006666667
    const double timeStepSize = parameters.timeStepSize; //period / static_cast<double>( timeStepMax );
    const unsigned writeEveryStep = parameters.writeEveryStep;
    double time = 0.;
    unsigned timeStep = 0;

    // Initialise mesh for flapping motion.
    {
        typedef boost::function< VecDim( const VecDim & ) >         GNGetVecDimFun;
        typedef boost::function< void( GNode *, VecDim & ) >        GNSetVecDimFun;
        typedef corlib::SetData< GNGetVecDimFun, GNSetVecDimFun >   GNSetData;
        GNSetData initialiseDisplacements( std::bind( &app::flappingDisplacements, 
                                                      std::placeholders::_1, period, time ),
                                           std::bind( &GNode::setDisplacements, 
                                                      std::placeholders::_1,
                                                      std::placeholders::_2 ) );
        gMesh.iterateOverNodes( initialiseDisplacements );
        GNSetData initialiseVelocities( std::bind( &app::flappingVelocities, 
                                                   std::placeholders::_1, period, time ),
                                        std::bind( &GNode::setVelocities, 
                                                   std::placeholders::_1,
                                                   std::placeholders::_2 ) );
        gMesh.iterateOverNodes( initialiseVelocities );
    }

    // Dirichlet boundary conditions
    typedef boost::function< corlib::NodalConstraint( VecDim ) >  GFunType;
    GFunType gSuppFun = std::bind( &app::supportFlapping, std::placeholders::_1, period, time, 0. );
    
    typedef corlib::ConstraintsFromFunction< GNode, GFunType >             GDirichlet;
    GDirichlet dirichlet( gSuppFun );
    gMesh.iterateOverNodes( dirichlet );

    // write initial state to VTU file
    const std::string vtuFile = animation.snapshotName( 0., 0 );
    app::writeMeshData( vtuFile, &gMesh );

    // write initial history
    app::writeHistory( histFile, 0., 0, gMesh );

    // create monitor
    std::ifstream monitorIn( parameters.monitorInputFileName.c_str() );
    FTL_VERIFY_DESCRIPTIVE( monitorIn.is_open(),
                            "Could not open file %s\n", parameters.monitorInputFileName.c_str() );
    std::ofstream monitorOut( parameters.monitorOutputFileName.c_str() );
    FTL_VERIFY_DESCRIPTIVE( monitorOut.is_open(),
                            "Could not open file %s\n", parameters.monitorOutputFileName.c_str() );
    typedef gshell::fem::Monitor< GNode, GElement >            GMonitor;
    GMonitor monitor( gMesh, &GElement::giveDisplacementAtNode,
                      monitorIn, monitorOut, 1.e-6 );
    gMesh.iterateOverElements( std::bind( monitor, std::placeholders::_1, 0, 0. ) );

    //==========================================================================
    // the time loop
    while ( time < (timeMax + 0.5*timeStepSize) and timeStep < timeStepMax ) {

        // a word to the user
        std::cout << std::endl << "time=" << time << ", timeStep=" << timeStep
                  << ", timeStepSize=" << timeStepSize << std::endl;
        const clock_t timerBegin = clock();

        // update Dirichlet boundary conditions --- non-constant
        gMesh.iterateOverNodes( std::bind( &GNode::clearConstraints, 
                                           std::placeholders::_1 ) );
        GFunType gSuppIncrFun = std::bind( &app::supportFlapping, 
                                           std::placeholders::_1,
                                           period, time, timeStepSize );
        GDirichlet incrementDirichlet( gSuppIncrFun );
        gMesh.iterateOverNodes( incrementDirichlet );

        //----------------------------------------------------------------------
        // Newton-Raphson iteration
        const double tolerance = 1.e-3;
        const unsigned maxiter = 10;
        bool converged = false;

        for ( unsigned i=0; i<maxiter and not converged; ++i ) {

            // a word to the user
            std::cout << "Iteration number " << i << std::endl;

            // Clear the element matrices and nodal forces
            std::cout << "    Clear nodal vectors ... " << std::flush;  timer.start();
#ifdef PARDISO_MKL
            gshell::fem::for_each_threaded( gMesh.nodesBegin(), gMesh.nodesEnd(),
                                            std::bind( &GNode::clearForce, std::placeholders::_1 ) );
#else
            gMesh.iterateOverNodes( std::bind( &GNode::clearForce, std::placeholders::_1 ) );
#endif
            std::cout << app::secondsToPosixTime( timer.elapsed().wall * 1e-9 ) << std::endl;

            // compute residual forces
            std::cout << "    Compute internal forces ... " << std::flush;  timer.start();
            GIntegrand internalForce( &GElement::internalForceIntegrand );
            corlib::Integrator< GQuad, GIntegrand > integrator2( quadrature, internalForce );
#ifdef PARDISO_MKL
            gshell::fem::for_each_threaded( gMesh.elementsBegin(), gMesh.elementsEnd(),
                                            integrator2 );
#else
            gMesh.iterateOverElements( integrator2 );
#endif
            std::cout << app::secondsToPosixTime( timer.elapsed().wall * 1e-9 ) << std::endl;

            // compute body forces
            {
                std::cout << "    Compute body forces ... " << std::flush;  timer.start();
                bodyForce.setFactor( 1. );
                corlib::Integrator< GQuad, GBodyForce > bodyForceIntegrator( quadrature, bodyForce );
                gMesh.iterateOverElements( bodyForceIntegrator );
                std::cout << app::secondsToPosixTime( timer.elapsed().wall * 1e-9 ) << std::endl;
            }

            // apply concentrated loads (and their tangent)
            //concentratedLoad.setFactor( time );
            //gMesh.iterateOverElements( concentratedLoad );

            // clear the stiffness matrix and RHS
            sysmat->clearMatrix();
            sysmat->clearRhs();

            // assemble residual forces
            std::cout << "    Assemble residual ... " << std::flush;  timer.start();

            // assemble static residual forces to the system matrices rhs
            gMesh.iterateOverNodes( corlib::vectorAssemblerFun( &GNode::giveForce,
                                                                &GNode::copyDofArray,
                                                                sysmat ) );

            // add inertia forces
            GIntegrand2 massIntegrand( &GElement::massIntegrand );
            GMatrixComputer mass( quadrature, massIntegrand );
#if   defined(DYNAMICS_BACKWARDEULER)
            corlib::BackwardEulerInertiaForceAssembler< GSysSolve, GElement, GMatrixComputer >
                nfa( timeStepSize, mass );
#elif defined(DYNAMICS_NEWMARK)
            corlib::NewmarkInertiaForceAssembler< GSysSolve, GElement, GMatrixComputer >
                nfa( newmarkBeta, newmarkGamma, timeStepSize, mass );
#endif
#ifdef PARDISO_MKL
            gshell::fem::for_each_threaded( gMesh.elementsBegin(), gMesh.elementsEnd(),
                                            std::bind( nfa, std::placeholders::_1, sysmat ) );
#else
            gMesh.iterateOverElements( std::bind( nfa, std::placeholders::_1, sysmat ) );
#endif

            // finish assemble residual
            std::cout << app::secondsToPosixTime( timer.elapsed().wall * 1e-9 ) << std::endl;

            // assemble the element stiffness matrices to the system matrix
            std::cout << "    Compute and assemble matrix ... " << std::flush;  timer.start();
            GIntegrand2 stiffIntegrand( &GElement::stiffnessIntegrand );
            GMatrixComputer stiff( quadrature, stiffIntegrand );
#if   defined(DYNAMICS_BACKWARDEULER)
            corlib::BackwardEulerMatrixAssembler< GSysSolve, GElement, GMatrixComputer, GMatrixComputer >
                nma( timeStepSize, mass, stiff );
#elif defined(DYNAMICS_NEWMARK)
            corlib::NewmarkMatrixAssembler< GSysSolve, GElement, GMatrixComputer, GMatrixComputer >
                nma( newmarkBeta, newmarkGamma, timeStepSize, mass, stiff );
#endif
#ifdef PARDISO_MKL
            gshell::fem::for_each_threaded( gMesh.elementsBegin(), gMesh.elementsEnd(),
                                            std::bind( nma, std::placeholders::_1, sysmat ) );
#else
            gMesh.iterateOverElements( std::bind( nma, std::placeholders::_1, sysmat ) );
#endif

            // finish assemble LHS
            std::cout << app::secondsToPosixTime( timer.elapsed().wall * 1e-9 ) << std::endl;

            // process the matrix and right hand side
            sysmat->finishAssembly();

            // collect and apply constraints
            std::cout << "    Applying constraints ... " << std::flush;  timer.start();
            double factor = ( i==0 ? 1. : 0. );
            corlib::ConstraintFunctor< GNode, GSysSolve > constraint( sysmat, factor );
            constraint = gMesh.iterateOverNodes( constraint );
            constraint.applyConstraints();
            std::cout << app::secondsToPosixTime( timer.elapsed().wall * 1e-9 ) << std::endl;

            // norm of the residual
            const double normResidual = sysmat->normRhs();
            std::cout << "    Norm of residual  " << std::scientific << normResidual << std::endl;
            if ( normResidual < tolerance )  converged = true;

            // solve the linear system
            std::cout << "    Solve linear system ... " << std::flush;  timer.start();
            sysmat->solveSystem();
            std::cout << app::secondsToPosixTime( timer.elapsed().wall * 1e-9 ) << std::endl;

            // add nodal solutions to increment
            std::cout << "    Distribute solution ... " << std::flush;  timer.start();
            gMesh.iterateOverNodes( corlib::distributorFun( sysmat,
                                                            &GNode::copyDofArray,
                                                            &GNode::addToIncrement ) );
            std::cout << app::secondsToPosixTime( timer.elapsed().wall * 1e-9 ) << std::endl;

            const double normDeltaU = sysmat->normSol();
            std::cout << "    Norm of delta U  " << std::scientific << normDeltaU << std::endl;
            if ( normDeltaU < tolerance )  converged = true;

        } // Newton

        // inform user
        if ( not converged ) {
            std::cout << "WARNING: Not converged in " << maxiter << " steps. Continuing." << std::endl;
        }

        // update velocities and accelerations
#if   defined(DYNAMICS_BACKWARDEULER)
        corlib::BackwardEulerSolutionDistributor< GNode > nmc( timeStepSize );
#elif defined(DYNAMICS_NEWMARK)
        corlib::NewmarkSolutionDistributor< GNode > nmc( newmarkBeta, newmarkGamma,
                                                         timeStepSize );
#endif
        gMesh.iterateOverNodes( nmc );

        // update displacements
        std::cout << "Update nodal displacements" << std::endl;
        gMesh.iterateOverNodes( std::bind( &GNode::updateDisplacements, std::placeholders::_1 ) );

        // monitor
        gMesh.iterateOverElements( std::bind( monitor, std::placeholders::_1, timeStep+1, time+timeStepSize ) );

        // write data to VTU file
        if ( (timeStep+1) % writeEveryStep == 0 ) {
            std::cout << "Writing VTU file" << std::endl;
            const std::string vtuFile = animation.snapshotName( time+timeStepSize, timeStep+1 );
            app::writeMeshData( vtuFile, &gMesh );
        }

        // write history
        app::writeHistory( histFile, time+timeStepSize, timeStep+1, gMesh );

        // increment
        time += timeStepSize;
        timeStep += 1;

        // time step timing
        const clock_t timerEnd = clock();
        const double timerDiff = timerEnd - timerBegin;
        const double timerDiffMs = (timerDiff*1000)/CLOCKS_PER_SEC;
        std::cout << "Time elapsed " << timerDiffMs << "ms" << std::endl;
    }

    //------------------------------------------------------------
    // write the animation file
    std::string animationFile( "animation.pvd" );
    std::ofstream anim( animationFile.c_str() );
    animation.writeAnimationFile( anim );
    anim.close();

    histFile.close();

    //--------------------------------------------------------------------------
    // clean-up
    delete sysmat;
    mc->destroy();
    delete subdiv;

    return 0;
}


