// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//
//                 Computational Structural Mechanics Lab
//                         University of Cambridge
//
// This software is copyrighted and all rights are retained by the University of
// Cambridge. All use, disclosure, and/or reproduction of any part not expressly
// authorised by the University of Cambridge is prohibited. (C) 2010
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#include <Eigen/Core>

// corlib includes
#include <corlib/eigenX.hpp>
#include <corlib/Shapefun.hpp>
#include <corlib/Quadrature.hpp>
#include <corlib/SystemSolve.hpp>
#include <corlib/SystemSolveEigenSparse.hpp>
#include <corlib/Integrator.hpp>
#include <corlib/Accessor.hpp>
#include <corlib/MatrixAssembler.hpp>
#include <corlib/VectorAssembler.hpp>
#include <corlib/Distributor.hpp>
#include <corlib/DofNumberer.hpp>
#include <corlib/TensorSpline.hpp>

//#define BACKWARDEULER
#ifdef BACKWARDEULER
#include <corlib/BackwardEulerAssembler.hpp>
#else
#include <corlib/NewmarkAssembler.hpp>
#endif
#include <corlib/Constraints.hpp>
#include <corlib/BodyForce.hpp>
#include <corlib/PropertiesParser.hpp>
#include <corlib/VTUwriter.hpp>
#include <corlib/Sensor.hpp>

// beam::fem includes
#include <beam/fem/NodeDynamic.hpp>
#include <beam/fem/ElementPlanarDynamic.hpp>
#include <beam/fem/Mesh.hpp>
#include <beam/fem/Supports.hpp>
#include <beam/fem/ConcentratedLoad.hpp>
#include <beam/fem/HighResVTU.hpp>

// solid::material includes 
#include <solid/material/MaterialContainer.hpp>
#include <solid/material/MaterialBase.hpp>


// local includes
#include "Loads.hpp"

// system includes
#include <iostream>
#include <fstream>
#include <string>
#include <boost/function.hpp>

namespace eigenX = corlib::eigenX;

//------------------------------------------------------------------------------
template< typename MESH>
void writeMeshData( const std::string & vtuFile, MESH * mesh );

template< typename MESH>
void writeMeshDataHiRes( const std::string & vtuFile, MESH * mesh );

//==============================================================================
int main( int argc, const char* argv[] )
{ 
    std::string meshFile, materialFile, constraintsFile, supportsFile, loadsFile;
    double beamArea, beamSecMomArea;
    double newmarkBeta = 0.25, newmarkGamma = 0.5;
    double timeStepSize;
    unsigned numTimeSteps;

    // feed properties parser with the variables to read
    corlib::PropertiesParser *prop = new corlib::PropertiesParser;
    prop->registerPropertiesVar( "meshFile",          meshFile   );
    prop->registerPropertiesVar( "materialFile",      materialFile );
    prop->registerPropertiesVar( "supportsFile",      supportsFile );
    prop->registerPropertiesVar( "loadsFile",         loadsFile );
    prop->registerPropertiesVar( "beamArea",          beamArea );
    prop->registerPropertiesVar( "beamSecMomArea",    beamSecMomArea );
    prop->registerPropertiesVar( "newmarkBeta",       newmarkBeta );
    prop->registerPropertiesVar( "newmarkGamma",      newmarkGamma );
    prop->registerPropertiesVar( "timeStepSize",      timeStepSize );
    prop->registerPropertiesVar( "numTimeSteps",      numTimeSteps );

    // read variables from the input.dat file
    std::ifstream inputFile( "./input.dat" );
    FTL_VERIFY_DESCRIPTIVE( inputFile.is_open(), "Cannot open input file" );
    prop->readValues( inputFile );    
    delete prop;
    inputFile.close( );

    // read material
    solid::material::MaterialContainer * mc = solid::material::MaterialContainer::instance();
    std::ifstream mf( materialFile.c_str( ) );
    FTL_VERIFY_DESCRIPTIVE( mf.is_open( ), "Cannot open material file" );    
    mc->readMaterialStream( mf );
    mf.close( );

    //! set some typedefs
    const unsigned numGaussPoints = 2;
    typedef corlib::Quadrature< corlib::LINE, numGaussPoints >       QuadType;
    const unsigned splineDegree = 3;
    typedef corlib::TensorSpline<splineDegree,1>                     SFunType;
    const unsigned dimEmbSpace = 2;
    const unsigned numDof = 2;
    typedef beam::fem::NodeDynamic< dimEmbSpace, numDof >            Node;
    typedef solid::material::MaterialBase                            MaterialType;
    typedef beam::fem::ElementPlanarDynamic< Node,
                                             SFunType,
                                             MaterialType >          Element;
    typedef beam::fem::LinkBasic< Node >                             LinkType;
    typedef beam::fem::Mesh< Element >                               MeshType;
    typedef corlib::SystemSolveEigenSparse                           SysSolve;
    typedef boost::function< void( Element*,
                                   const Element::VecLDim &,
                                   const double & ) >                Integrand;

    typedef eigenX::VectorSd< numDof >                               Vec2;
    typedef std::function<Vec2( Vec2 )>                              FunType;
    typedef corlib::BodyForce< Element, FunType >                    BodyForce;
    typedef beam::fem::ConcentratedLoad< Element >                   ConcentratedLoad;

    std::ifstream smf( meshFile.c_str( ) );
    FTL_VERIFY( smf.is_open() );
    MeshType mesh( smf );
    smf.close( );

    // obtain the 1st material
    solid::material::MaterialBase* material = mc->getMaterial( 0 );
    //! cache material to elements
    mesh.iterateOverElements( std::bind( &Element::cacheMaterial, std::placeholders::_1, material ) );

    // distribute cross section geometry
    mesh.iterateOverElements( std::bind( &Element::setArea,       std::placeholders::_1, beamArea ) );
    mesh.iterateOverElements( std::bind( &Element::setSecMomArea, std::placeholders::_1, beamSecMomArea ) );

    //! read links
    std::ifstream suppFile( supportsFile.c_str( ) );
    FTL_VERIFY( suppFile.is_open( ) );
    beam::fem::SupportsFromFile suppReader( suppFile );
    suppFile.close();
    mesh.setLinksOfSupports( suppReader );

    //! body force
    BodyForce bodyForce( &gravity );

    //! concentrated loads
    std::ifstream ldsFile( loadsFile.c_str( ) );
    FTL_VERIFY( ldsFile.is_open() );
    ConcentratedLoad concentratedLoad( ldsFile, mesh );
    ldsFile.close();

    //! Set up a sensor listing on displacements
    typedef std::function< Node::VecDof( const Node* )> SensorFun;
    SensorFun sensorFun = std::bind( &Node::giveDisplacements, std::placeholders::_1 );
    corlib::Accessor< SensorFun > accessor( sensorFun );
    corlib::Sensor< corlib::Accessor<SensorFun> > nodeSensor( mesh.getNodePointer( 8 ), accessor );
    
    // quadrature object
    QuadType quadrature;

    // volume of ring
    Integrand volumeIntegrand( &Element::volumeIntegrand );
    corlib::Integrator< QuadType, Integrand, Element > integratorVolume( quadrature, volumeIntegrand );
    mesh.iterateOverElements( integratorVolume );
    double volume = 0.0;
    mesh.iterateOverElements( std::bind2nd( std::mem_fun( &Element::addVolume ), &volume ) );
    std::cout << "volume=" << volume << std::endl;
    //std::cout << "volume=" << beamArea*6.0*M_PI << std::endl;

    //! number dofs
    unsigned numdofs = mesh.iterateOverNodes( corlib::dofNumberFun( &Node::numberDOFs, 0 ) );
    std::cout << numdofs << " degrees of freedom" << std::endl;
    // number Lagrange multipliers
    numdofs = mesh.iterateOverLinks( corlib::dofNumberFun( &LinkType::numberDOFs, numdofs ) );

    std::cout << numdofs << " degrees of freedom" << std::endl;

    // set up an animation writer
    int check = system( "mkdir -p animation" ); FTL_VERIFY( !check );
    std::string baseName = meshFile.substr( 0, meshFile.rfind( "." ) );
    corlib::VTUanim animation( "./animation/" + baseName, "vtu" );

#ifdef D_HIRES
    int checkHiRes = system( "mkdir -p animation-hr" ); FTL_VERIFY( !checkHiRes );
    corlib::VTUanim animationHiRes( "./animation-hr/" + baseName + "-hr", "vtu" );
#endif
    
    // write data to VTU file
    std::string vtuFile = animation.snapshotName( 0., 0 );
    writeMeshData( vtuFile, &mesh );

#ifdef D_HIRES
    std::string vtuFileHiRes = animationHiRes.snapshotName( 0., 0 );
    writeMeshDataHiRes( vtuFileHiRes, &mesh );
#endif

    // Compute element mass matrices. The matrices are stored and are constant
    // during the computation
    {
        // Clear the element matrices and nodal forces
        mesh.iterateOverElements( std::bind( &Element::clearMass, std::placeholders::_1 ) );
        
        //! compute element mass matrices
        Integrand massIntegrand( &Element::massIntegrand );
        corlib::Integrator< QuadType, Integrand, Element >
            massIntegrator( quadrature, massIntegrand );
        mesh.iterateOverElements( massIntegrator );
    }

    // system matrix object
    SysSolve* sysmat = new SysSolve( numdofs );
    const double tolerance = 1e-8;
    const unsigned maxiter = 30;

    // time control parameters
    const double loadFactor = 1.0;
    const double timeMax = numTimeSteps * timeStepSize;
    double time = timeStepSize;

    // a word to the user
#ifdef BACKWARDEULER
    std::cout << "Runnning Backward-Euler" << std::endl;
#else
    std::cout << "Runnning Newmark (beta=" << newmarkBeta
              << ", gamma=" << newmarkGamma << ")" << std::endl;
#endif

    // time loop
    unsigned stepNum = 0;
    while ( time < (timeMax+0.5*timeStepSize) ) {
        // a word to the user
        std::cout << std::endl << "time=" << time << std::endl;

        // Newton-Raphson iteration
        bool converged = false;
        for ( unsigned i=0; i<maxiter and not converged; ++i ) {
            // a word to the user
            std::cout << "Iteration number " << i << std::endl;
        
            // Clear the element matrices and nodal forces
            mesh.iterateOverElements( std::bind( &Element::clearStiffness, std::placeholders::_1 ) );

            //! compute element stiffness matrices
            Integrand stiffness( &Element::stiffnessIntegrand );
            corlib::Integrator< QuadType, Integrand, Element > integrator( quadrature, stiffness );
            mesh.iterateOverElements( integrator );

            // Blank system matrix
            sysmat->clearMatrix();

            //! assemble the element mass and stiffness matrices to the system matrix
            typedef std::function< void( const Element*,
                                         Eigen::MatrixXd & )> MatrixDonator;
            MatrixDonator massMat( &Element::addMassMatrix );
            MatrixDonator stiffnessMat( &Element::addStiffnessMatrix );

#ifdef BACKWARDEULER
            corlib::BackwardEulerMatrixAssembler<SysSolve,Element,MatrixDonator,MatrixDonator> 
                nma( timeStepSize, massMat, stiffnessMat );
#else
            corlib::NewmarkMatrixAssembler<SysSolve,Element,MatrixDonator,MatrixDonator> 
                nma( newmarkBeta, newmarkGamma, timeStepSize, massMat, stiffnessMat );
#endif
            mesh.iterateOverElements( std::bind( nma, std::placeholders::_1, sysmat )  );

            // Clear element nodal forces
            mesh.iterateOverNodes( std::bind( &Node::clearForce, std::placeholders::_1 ) );

            // compute internal forces
            Integrand internalForce( &Element::internalForceIntegrand );
            corlib::Integrator<QuadType,Integrand,Element> integrator2( quadrature, internalForce );
            mesh.iterateOverElements( integrator2 );

            // apply body forces (external forces)
            bodyForce.setFactor( loadFactor );
            corlib::Integrator<QuadType,BodyForce> bodyForceIntegrator( quadrature, bodyForce );
            mesh.iterateOverElements( bodyForceIntegrator );

            // apply concentrated loads (and their tangent)
            //concentratedLoad.setFactor( loadFactor );
            //mesh.iterateOverElements( concentratedLoad );

            //! clear the global RHS
            sysmat->clearRhs( );

            //! assemble the static residual forces to the system matrices rhs
            mesh.iterateOverNodes( corlib::vectorAssemblerFun( &Node::giveForce,
                                                               &Node::copyDofArray,
                                                               sysmat ) );

            // add negative inertia forces to RHS
#ifdef BACKWARDEULER
            corlib::BackwardEulerInertiaForceAssembler<SysSolve,Element,MatrixDonator> 
                nfa( timeStepSize, massMat );
#else
            corlib::NewmarkInertiaForceAssembler<SysSolve,Element,MatrixDonator> 
                nfa( newmarkBeta, newmarkGamma, timeStepSize, massMat );
#endif
            mesh.iterateOverElements( std::bind( nfa, std::placeholders::_1, sysmat ) );

            //! assemble link matrices to the system matrix
            corlib::MatrixGiveAndAssembleFun<const LinkType,SysSolve>
                couplingMat( &LinkType::giveCouplingMatrix, 
                             &LinkType::getDofIndices, 
                             &LinkType::getDofIndicesP,
                             sysmat );
            mesh.iterateOverLinks( couplingMat );

            corlib::MatrixGiveAndAssembleFun<const LinkType,SysSolve>
                couplingMatT( &LinkType::giveCouplingMatrixT, 
                              &LinkType::getDofIndicesP, 
                              &LinkType::getDofIndices,
                              sysmat );
            mesh.iterateOverLinks( couplingMatT );

            mesh.iterateOverLinks( corlib::vectorAssemblerFun( &LinkType::giveMultiplierRhs,
                                                               &LinkType::getDofIndices,
                                                               sysmat ) );
            mesh.iterateOverLinks( corlib::vectorAssemblerFun( &LinkType::giveLinkResiduum,
                                                               &LinkType::getDofIndicesP,
                                                               sysmat ) );

            //! collect and apply constraints
            sysmat->finishAssembly();

            //! norm of the residual
            const double normResidual = sysmat->normRhs( );
            std::cout << "Norm of residual  " << std::scientific << normResidual << std::endl;
            //if ( normResidual < tolerance ) converged = true;

            //! solve the linear system
            sysmat->solveSystem();

            //! add nodal solutions to increment
            mesh.iterateOverNodes( corlib::distributorFun( sysmat, 
                                                           &Node::copyDofArray, 
                                                           &Node::addToIncrement ) );
            mesh.iterateOverLinks( corlib::distributorFun( sysmat, 
                                                           &LinkType::getDofIndicesP, 
                                                           &LinkType::addToIncrement ) );
            // check convergence
            const double normDeltaU = sysmat->normSol( );
            std::cout << "Norm of delta U  " << std::scientific << normDeltaU << std::endl;
            if ( normDeltaU < tolerance ) {
                converged = true;
            }

        }

        // collect nodal solutions and update
#ifdef BACKWARDEULER
        corlib::BackwardEulerSolutionDistributor<Node> 
            collector( timeStepSize  );
#else
        corlib::NewmarkSolutionDistributor< Node> 
            collector( newmarkBeta, newmarkGamma, timeStepSize  );
#endif
        mesh.iterateOverNodes( collector );

        //! update displacements
        mesh.iterateOverNodes( std::bind( &Node::updateDisplacements, std::placeholders::_1 ) );
        // update Lagrange multipliers
        mesh.iterateOverLinks( std::bind( &LinkType::updateMultiplier, std::placeholders::_1 ) );

        // inform user
        if ( not converged ) {
            std::cout << "WARNING: Not converged in " << maxiter << " steps. Continuing." << std::endl;
        }

        // let sensor store time and 0-th component
        nodeSensor.listenComponent( time, 0 );

        // write data to VTU file
        std :: string vtuFile = animation.snapshotName( time, stepNum+1 );
        writeMeshData( vtuFile, &mesh );

#ifdef D_HIRES
        std :: string vtuFileHiRes = animationHiRes.snapshotName( time, stepNum+1 );
        writeMeshDataHiRes( vtuFileHiRes, &mesh );
#endif

        // increment time
        time += timeStepSize;
        stepNum ++;

    }

    //------------------------------------------------------------
    //! let the sensor write its data
    const std::string gpFile = baseName + ".gp";
    std::ofstream gp( gpFile.c_str() );
    FTL_VERIFY( gp.is_open( ) );
    nodeSensor.write( gp );
    gp.close( );

    //------------------------------------------------------------
    // write the animation file
    std::string animationFile( "animation.pvd" );
    std::ofstream anim( animationFile.c_str( ) );
    animation.writeAnimationFile( anim );
    anim.close();

#ifdef D_HIRES
    std::string animationFileHiRes( "animation-hr.pvd" );
    std::ofstream animHiRes( animationFileHiRes.c_str( ) );
    animationHiRes.writeAnimationFile( animHiRes );
    animHiRes.close();
#endif

    //------------------------------------------------------------
    //! free the memory
    mc->destroy( );
    delete sysmat;

    return 0;

}

//------------------------------------------------------------------------------
// write solution data to a file
template< typename MESH>
void writeMeshData( const std::string & vtuFile, MESH * mesh )
{
    typedef          MESH       Mesh;
    typedef typename Mesh::Node Node;

    // set up the VTU writer and open a file
    std::ofstream vtu( vtuFile.c_str( ) );
    FTL_VERIFY( vtu.is_open() );         
    corlib::VTUwriter< Mesh > vtuwriter( mesh, vtu );

    //! write mesh data
    vtuwriter.writeMesh( );
    
    //-- write point data
    vtuwriter.openPointData( );
    vtuwriter.writeNodalQuantity( corlib::accessorFun( &Node::giveDisplacements, "Displacements" ) );
    vtuwriter.writeNodalQuantity( corlib::accessorFun( &Node::giveVelocities   , "Velocities"    ) );
    vtuwriter.writeNodalQuantity( corlib::accessorFun( &Node::giveAccelerations, "Accelerations" ) );
    vtuwriter.closePointData( );

    // finish writing
    vtuwriter.finalize( );
    vtu.close( );

    return;

}

//------------------------------------------------------------------------------
// write solution data to a file
template< typename MESH>
void writeMeshDataHiRes( const std::string & vtuFile, MESH * mesh )
{
    const unsigned resolution = 10;

    std::ofstream vtu(  vtuFile.c_str() );
    beam::fem::HighResVTU<MESH> vtuWriter( mesh, vtu, resolution );
    vtuWriter.writeMesh();

    vtuWriter.openPointData();
    vtuWriter.writePointQuantity( corlib::accessorFun( &MESH::Node::giveDisplacements, "Displacements" ) );
    vtuWriter.writePointQuantity( corlib::accessorFun( &MESH::Node::giveVelocities   , "Velocities"    ) );
    vtuWriter.writePointQuantity( corlib::accessorFun( &MESH::Node::giveAccelerations, "Accelerations" ) );

    vtuWriter.closePointData();


    //{
    //    vtuWriter.openCellData();
    //    boost::function<typename MESH::Element::VecDim( const typename MESH::Element*, 
    //                                                    const typename MESH::Element::VecLDim &) > getNormal 
    //        = boost::bind( &MESH::Element::giveUnitNormalCur, _1, _2 );
    //    vtuWriter.writeCellQuantity( getNormal, "Normal" );
    //    vtuWriter.closeCellData();
    //}


    vtuWriter.finalize();
    vtu.close();

    return;

}
