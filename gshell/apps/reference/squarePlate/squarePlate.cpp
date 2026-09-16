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

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <functional>

// includes from Eigen
#include <Eigen/Core>

// includes from corlib
#include <corlib/eigenX.hpp>
#include <corlib/PropertiesParser.hpp>
#include <corlib/verify.hpp>
#include <corlib/Accessor.hpp>
#include <corlib/VTUwriter.hpp>
#include <corlib/fuzzyEqual.hpp>
#include <corlib/Quadrature.hpp>
#include <corlib/Integrator.hpp>
#include <corlib/SystemSolve.hpp>
#include <corlib/SystemSolveEigenSparse.hpp>
#include <corlib/DofNumberer.hpp>
#include <corlib/Distributor.hpp>
#include <corlib/MatrixAssembler.hpp>
#include <corlib/VectorAssembler.hpp>
#include <corlib/BodyForce.hpp>
#include <corlib/AccumulateQuantity.hpp>

// includes from subdiv::surf
#include <subdiv/surf/Vertex.hpp>
#include <subdiv/surf/Edge.hpp>
#include <subdiv/surf/Facet.hpp>
#include <subdiv/surf/FacetTree.hpp>
#include <subdiv/surf/Subdivision.hpp>
#include <subdiv/surf/Mesh.hpp>
#include <subdiv/surf/ShapeFunSubdivision.hpp>

// includes from solid::material
#include <solid/material/MaterialContainer.hpp>
#include <solid/material/svenant/SVenant.hpp>

// includes from gshell::fem
#include <gshell/fem/NodeShell.hpp>
#include <gshell/fem/NodeStatic.hpp>
#include <gshell/fem/ShapeFunCache.hpp>
#include <gshell/fem/ElementShell.hpp>
#ifdef KINEMATICS_KIRCHHOFFLOVE
#include <gshell/fem/ElementKLStatic.hpp>
#else
#include <gshell/fem/ElementStatic.hpp>
#endif
#include <gshell/fem/MeshShell.hpp>
//#include <gshell/fem/DebugPlot.hpp>
#include <gshell/fem/Monitor.hpp>

#ifdef KINEMATICS_KIRCHHOFFLOVE
const unsigned DOF = 3;
#else
const unsigned DOF = 5;
#endif

// local includes
#include "Parameters.hpp"
#include "DirichletFun.hpp"
#include "write.hpp"

// helper functions
corlib::eigenX::VectorSd<3> 
pressure( const corlib::eigenX::VectorSd<3> x );


//==============================================================================
// everything starts here
int main( int argc, char * args[] ) {

    // uBLAS
    namespace eigenX = corlib::eigenX;

    //-----------------------------------------------<---------------------------
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
#if   defined(LOOP_)          && defined(MESH_A)
    const std::string inputData = "inputLoopA.dat";
#elif defined(LOOP_)          && defined(MESH_B)
    const std::string inputData = "inputLoopB.dat";
#elif defined(LOOP_)          && defined(MESH_D)
    const std::string inputData = "inputLoopD.dat";
#elif defined(CATMULL_CLARK_) && defined(MESH_CC)
    const std::string inputData = "inputCC.dat";
#else
    const std::string inputData = "input.dat";
#endif
    std::ifstream inp( inputData.c_str() );
    FTL_VERIFY_DESCRIPTIVE( inp.is_open(), 
                            "Could not open file %s\n",
                            inputData.c_str() );
    app::Parameters parameters( inp );
    inp.close();

    //read input mesh and tags
    std::ifstream inputStream( parameters.inputMeshFileName.c_str() );
    std::ifstream secondaryInputStream( parameters.secondaryInputFileName.c_str() );
    FTL_VERIFY_DESCRIPTIVE( inputStream.is_open(), 
                            "Could not open file %s\n",
                            parameters.inputMeshFileName.c_str() );
    FTL_VERIFY_DESCRIPTIVE( secondaryInputStream.is_open(),
                            "Could not open file %s\n",
                            parameters.secondaryInputFileName.c_str() );

    // Surface
    Mesh sMesh( inputStream );
    sMesh.readTags( secondaryInputStream );  // before buildFacetTopo
    sMesh.buildFacetTopology();  // OK

    // Subdivision scheme
    Subdivision * subdiv = new Subdivision();

//#define PLAIN_SUBDIVISION
#ifdef PLAIN_SUBDIVISION
    // Make output folder
    const std::string makeOutputFolder = "mkdir -p " + parameters.outputFolder;
    FTL_VERIFY( !system( makeOutputFolder.c_str() ) );
    // create facet trees
    sMesh.implantFacetTrees();
    // Make complete subdivision
    const unsigned numSubdivisions = 4;
    for ( unsigned s = 0; s <= numSubdivisions; ++s ) {
        // subdivide
        std::cout << s << ". subdivision" << std::endl;
        if ( s > 0 ) {
            sMesh.subdivide< Subdivision >( subdiv );
        }
        // write resulting file
        const std::string outputSmfFileName = parameters.outputBaseFileName +
            "." + std::to_string( s ) + ".smf";
        const std::string outputTgFileName = parameters.outputBaseFileName +
            "." + std::to_string( s ) + ".tg";
        std::ofstream outputSmfStream( outputSmfFileName.c_str() );
        std::ofstream outputTgStream( outputTgFileName.c_str() );
        sMesh.writeMeshWithIndex( outputSmfStream, &outputTgStream );
        outputSmfStream.close();
        outputTgStream.close();
    }
#endif

#ifndef PLAIN_SUBDIVISION
    //--------------------------------------------------------------------------
    // create g-shell mesh
    const enum corlib::shape gShape = myShape;
    typedef gshell::fem::NodeShell                                     GBasisNode;
    typedef gshell::fem::NodeStatic< GBasisNode, DOF >                 GNode;
    typedef gshell::fem::ShapeFunCache< gShape >                       GShapeFunCache;
    typedef gshell::fem::ElementShell< GNode, GShapeFunCache >         GBasisElement;
    typedef solid::material::SVenant                                   GMaterial;
#ifdef KINEMATICS_KIRCHHOFFLOVE
    typedef gshell::fem::ElementKLStatic< GBasisElement, GMaterial >   GElement;
#else
    typedef gshell::fem::ElementStatic< GBasisElement, GMaterial >     GElement;
#endif
    typedef gshell::fem::MeshShell< GElement, ShapeFun >               GMesh;
#ifdef LOOP_
    typedef corlib::Quadrature< myShape, 3 >                           GQuad;
#elif CATMULL_CLARK_
    typedef corlib::Quadrature< myShape, 4 >                           GQuad;
#endif

    typedef GNode::VecDim                                              VecDim;
    typedef GNode::VecDof                                              VecDof;
    typedef GElement::VecLDim                                          VecLDim;

    //! keep integrand with boost - c.f. corlib/integrator
    typedef boost::function< void( GElement *, const VecLDim &,
                                   const double & ) >                  GIntegrand;
    
    typedef corlib::AccumulateQuantity< GElement, double >             GESumDouble;

    typedef corlib::SystemSolveEigenSparse                             GSysSolve;

    // the quadrature points and weights
    GQuad quadrature;

    // the g-mesh
    GMesh gMesh( sMesh, quadrature );

    // read material
    solid::material::MaterialContainer * mc = solid::material::MaterialContainer::instance();
    std::ifstream mf( parameters.materialFileName.c_str() );
    FTL_VERIFY_DESCRIPTIVE( mf.is_open(),
                            "Cannot open file %s\n", parameters.materialFileName.c_str() );
    mc->readMaterialStream( mf );
    mf.close();
    // set material
    GMaterial * material = dynamic_cast< GMaterial * >( mc->getMaterial( 0 ) );
    // cache material to elements
    gMesh.iterateOverElements( std::bind( &GElement::cacheMaterial, std::placeholders::_1, material ) );

    // set thickness
#ifdef KINEMATICS_KIRCHHOFFLOVE
    const double thickness = parameters.thicknessKL;
#else
    const double thickness = parameters.thicknessSF;
#endif
    gMesh.iterateOverElements( std::bind( &GElement::setThickness, std::placeholders::_1, thickness ) );

    // set unique tangents
    gMesh.iterateOverElements( std::bind( &GElement::setLimitCoefficientsAtNodes, std::placeholders::_1 ) );
    gMesh.iterateOverElements( std::bind( &GElement::setUniqueTangentCoefficientsAtNodes, std::placeholders::_1 ) );

    // number dofs
    const unsigned numdofs = gMesh.iterateOverNodes( corlib::dofNumberFun( &GNode::numberDOFs, 0 ) );
    std::cout << numdofs << " degrees of freedom" << std::endl;

    // set up an animation writer
    const std::string makeAnimationFolder = "mkdir -p " + parameters.animationFolder;
    FTL_VERIFY( !system( makeAnimationFolder.c_str() ) );
    corlib::VTUanim animation( parameters.animationBaseFileName, "vtu" );

    // allocate linear system solver
    GSysSolve * sysmat = new GSysSolve( numdofs );

    // body force
    typedef std::function< VecDim( VecDim ) >                   FunType2;
    typedef corlib::BodyForce< GElement, FunType2 >             GBodyForce;
    GBodyForce bodyForce( std::bind( &pressure, std::placeholders::_1 ) );

    // Dirichlet boundary conditions
    typedef std::function< corlib::NodalConstraint( VecDim ) >           GFunType;
    typedef corlib::ConstraintsFromFunction< GNode, GFunType >           GDirichlet;
    GDirichlet dirichlet( std::bind( &app::dirichletSimply, std::placeholders::_1 ) );
    gMesh.iterateOverNodes( dirichlet );
    //GDirichlet dirichlet2( boost::bind( &app::enforceZeroShears, _1 ) );
    //gMesh.iterateOverNodes( dirichlet2 );

    // integrate area
    {
        GIntegrand areaIntegrand( &GElement::areaIntegrand );
        corlib::Integrator<GQuad, GIntegrand, GElement> integratorArea( quadrature, areaIntegrand );
        gMesh.iterateOverElements( integratorArea );
        GESumDouble sumArea( &GElement::giveArea, 0.0 );
        double area = gMesh.iterateOverElements( sumArea );
        std::cout << "area=" << area << std::endl;
    }

    // write initial state to VTU file
    const std::string vtuFile = animation.snapshotName( 0.0, 0 );
    app::writeMeshData( vtuFile, &gMesh );

    // create monitor
    std::ifstream monitorIn( parameters.monitorInputFileName.c_str() );
    FTL_VERIFY_DESCRIPTIVE( monitorIn.is_open(),
                            "Could not open file %s\n", parameters.monitorInputFileName.c_str() );
    std::ofstream monitorOut( parameters.monitorOutputFileName.c_str() );
    FTL_VERIFY_DESCRIPTIVE( monitorOut.is_open(),
                            "Could not open file %s\n", parameters.monitorOutputFileName.c_str() );
    typedef gshell::fem::Monitor< GNode, GElement >   GMonitor;
    GMonitor monitor( gMesh, &GElement::giveDisplacementAtNode,
                      monitorIn, monitorOut, 1.e-6 );

    //==========================================================================
    // load control parameters
    const double loadMax = 1.0;
    const double loadStepSize = 1.0;
    double loadFactor = 1.0;
    unsigned loadStep = 1;

    // load control loop
    while ( loadFactor < (loadMax+0.5*loadStepSize) ) {

        // a word to the user
        std::cout << std::endl << "loadFactor=" << loadFactor << std::endl;
        
        //----------------------------------------------------------------------
        // Newton-Raphson iteration
        const double tolerance = 1.e-10;
        const unsigned maxiter = 1;
        bool converged = false;

        for ( unsigned i=0; i<maxiter and not converged; ++i ) {
            // a word to the user
            std::cout << "Iteration number " << i << std::endl;

            // Clear the element matrices and nodal forces
            gMesh.iterateOverNodes( std::bind( &GNode::clearForce, std::placeholders::_1 ) );

            // compute residual forces
            GIntegrand internalForce( &GElement::internalForceIntegrand );
            corlib::Integrator< GQuad, GIntegrand > integrator2( quadrature, internalForce );
            gMesh.iterateOverElements( integrator2 );

            // apply body forces
            bodyForce.setFactor( loadFactor );
            corlib::Integrator< GQuad, GBodyForce > bodyForceIntegrator( quadrature, bodyForce );
            gMesh.iterateOverElements( bodyForceIntegrator );

            // apply concentrated loads (and their tangent)
            //concentratedLoad.setFactor( loadFactor );
            //gMesh.iterateOverElements( concentratedLoad );

            // clear the stiffness matrix and RHS
            sysmat->clearMatrix();
            sysmat->clearRhs();

            // assemble the residual forces to the system matrices rhs
            gMesh.iterateOverNodes( corlib::vectorAssemblerFun( &GNode::giveForce, 
                                                                &GNode::copyDofArray, 
                                                                sysmat ) );

            // compute and assemble the element stiffness matrices to the system matrix
            corlib::MatrixComputeAndAssembleFun< const GElement, GQuad, GSysSolve >
                mafStiffness( &GElement::stiffnessIntegrand, 
                              &GElement::getDofIndices,
                              &GElement::getDofIndices,
                              quadrature, sysmat );
            gMesh.iterateOverElements( mafStiffness );

            // process the matrix and right hand side
            sysmat->finishAssembly();

            // collect and apply constraints
            double factor = ( i==0 ? 1. : 0. );
            corlib::ConstraintFunctor< GNode, GSysSolve > constraint( sysmat, factor );
            constraint = gMesh.iterateOverNodes( constraint );
            constraint.applyConstraints();

            // norm of the residual
            const double normResidual = sysmat->normRhs();
            std::cout << "Norm of residual  " << std::scientific << normResidual << std::endl;
            //if ( normResidual < tolerance ) converged = true;

            // solve the linear system
            sysmat->solveSystem();

            // add nodal solutions to increment
            gMesh.iterateOverNodes( corlib::distributorFun( sysmat, 
                                                            &GNode::copyDofArray, 
                                                            &GNode::addToIncrement ) );

            const double normDeltaU = sysmat->normSol();
            std::cout << "Norm of delta U  " << std::scientific << normDeltaU << std::endl;
            if ( normDeltaU < tolerance ) {
                converged = true;
            }

        } // Newton

        // inform user
        if ( not converged ) {
            std::cout << "WARNING: Not converged in " << maxiter << " steps. Continuing." << std::endl;
        }

        // update displacements
        gMesh.iterateOverNodes( std::bind( &GNode::updateDisplacements, std::placeholders::_1 ) );

        // monitor
        gMesh.iterateOverElements( std::bind( monitor, std::placeholders::_1, loadStep, loadFactor ) );

        // write data to VTU file
        const std::string vtuFile = animation.snapshotName( loadFactor, loadStep );
        app::writeMeshData( vtuFile, &gMesh );

        // increment
        loadFactor += loadStepSize;
        loadStep   += 1;
    }

    //------------------------------------------------------------
    // write the animation file
    std::string animationFile( "animation.pvd" );
    std::ofstream anim( animationFile.c_str() );
    animation.writeAnimationFile( anim );
    anim.close();

    //--------------------------------------------------------------------------
    // clean-up
    delete sysmat;
    mc->destroy();
#endif  // PLAIN_SUBDIVISION
    delete subdiv;

    return 0;
}

//==============================================================================
corlib::eigenX::VectorSd<3>
pressure( const corlib::eigenX::VectorSd<3> x )
{
    //const double E = 1.092e6;  // N/m^2
    //const double l = 10.;  // m
    //const double nu = 0.3;  // 1
    const double t = 1.e-5;   // m
    //const double pRes = E*t*t*t/(12.*(1.-nu*nu)*l*l*l*l); // N/m^3
    const double pRes = 1.e-15;  // N/m^3

    corlib::eigenX::VectorSd<3> p;
    p.setZero();
    p[ 2 ] = pRes / t;  // N/m^4
    return p;
}

