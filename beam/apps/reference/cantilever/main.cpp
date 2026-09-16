// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//
//                 Computational Structural Mechanics Lab
//                         University of Cambridge
//
// This software is copyrighted and all rights are retained by the University of
// Cambridge. All use, disclosure, and/or reproduction of any part not expressly
// authorised by the University of Cambridge is prohibited. (C) 2010
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

//! @file main.cpp

// corlib includes
#include <corlib/eigenX.hpp>
#include <corlib/Shapefun.hpp>
#include <corlib/Quadrature.hpp>
#include <corlib/SystemSolveEigenSparse.hpp>
#include <corlib/Integrator.hpp>
#include <corlib/Accessor.hpp>
#include <corlib/MatrixAssembler.hpp>
#include <corlib/VectorAssembler.hpp>
#include <corlib/Constraints.hpp>
#include <corlib/BodyForce.hpp>
#include <corlib/Distributor.hpp>
#include <corlib/DofNumberer.hpp>
#include <corlib/PropertiesParser.hpp>
#include <corlib/SmfHead.hpp>
#include <corlib/VTUwriter.hpp>
#include <corlib/TensorSpline.hpp>

// solid::material includes 
#include <solid/material/MaterialContainer.hpp>
#include <solid/material/MaterialBase.hpp>

// beam::fem includes
#include <beam/fem/NodeStatic.hpp>
#include <beam/fem/ElementPlanarStatic.hpp>
#include <beam/fem/Mesh.hpp>
#include <beam/fem/Supports.hpp>
#include <beam/fem/ConcentratedLoad.hpp>
#include <beam/fem/HighResVTU.hpp>

// local includes
#include "Loads.hpp"

// system includes
#include <iostream>
#include <fstream>
#include <string>

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

    // feed properties parser with the variables to read
    corlib::PropertiesParser *prop = new corlib::PropertiesParser;
    prop->registerPropertiesVar( "meshFile",        meshFile   );
    prop->registerPropertiesVar( "materialFile",    materialFile );
    prop->registerPropertiesVar( "constraintsFile", constraintsFile );
    prop->registerPropertiesVar( "supportsFile",    supportsFile );
    prop->registerPropertiesVar( "loadsFile",       loadsFile );
    prop->registerPropertiesVar( "beamArea",        beamArea );
    prop->registerPropertiesVar( "beamSecMomArea",  beamSecMomArea );

    // read variables from the input.dat file
    std::ifstream inputFile( "./input.dat" );
    FTL_VERIFY_DESCRIPTIVE( inputFile.is_open(), "Cannot open input file" );
    prop->readValues( inputFile );
    delete prop;
    inputFile.close( );

    // read material
    solid::material::MaterialContainer* mc = solid::material::MaterialContainer::instance();
    std::ifstream mf( materialFile.c_str( ) );
    FTL_VERIFY_DESCRIPTIVE( mf.is_open(), "Cannot open material file" );
    mc->readMaterialStream( mf );
    mf.close( );

    // set some typedefs
    const unsigned numGaussPoints = 2;
    typedef corlib::Quadrature< corlib::LINE, numGaussPoints >       QuadType;
    const unsigned splineDegree = 3;
    typedef corlib::TensorSpline<splineDegree,1>                     SFunType;
    const unsigned dimEmbSpace = 2;
    const unsigned numDof = 2;
    typedef beam::fem::NodeStatic< dimEmbSpace, numDof >             NodeType;
    typedef solid::material::MaterialBase                            MaterialType;
    typedef beam::fem::ElementPlanarStatic< NodeType,
                                            SFunType,
                                            MaterialType >           ElementType;
    typedef beam::fem::LinkBasic< NodeType >                         LinkType;
    typedef beam::fem::Mesh< ElementType >                           MeshType;
    typedef corlib::SystemSolveEigenSparse                           SysSolve;
    typedef boost::function<void( ElementType *,
                                  const ElementType::VecLDim &,
                                  const double & ) >                 Integrand;

    typedef eigenX::VectorSd< numDof >                               Vec2;
    typedef std::function< Vec2 ( Vec2 ) >                           FunType;
    typedef corlib::BodyForce< ElementType, FunType >                BodyForce;
    typedef beam::fem::ConcentratedLoad< ElementType >               ConcentratedLoad;

    std::ifstream smf( meshFile.c_str( ) );
    FTL_VERIFY( smf.is_open() );
    corlib::SmfHead smfHead;
    smfHead.readValidated( smf, corlib::LINE, 2 );
    MeshType mesh( smf );
    smf.close( );

    QuadType quadrature;

    // obtain the 1st material and cache to elements
    solid::material::MaterialBase* material = mc->getMaterial( 0 );
    mesh.iterateOverElements( boost::bind( &ElementType::cacheMaterial, _1, material ) );

    // cache cross section geometry
    mesh.iterateOverElements( boost::bind( &ElementType::setArea,  _1, beamArea ) );
    mesh.iterateOverElements( boost::bind( &ElementType::setSecMomArea, _1, beamSecMomArea ) );

    // read the constraints
    std::ifstream cstrFile( constraintsFile.c_str( ) );
    FTL_VERIFY( cstrFile.is_open( ) );
    corlib::ConstraintsFromFile< NodeType > constraintsHandler( cstrFile );
    cstrFile.close();
    mesh.iterateOverNodes( constraintsHandler );

    // read links or supports
    std::ifstream suppFile( supportsFile.c_str( ) );
    FTL_VERIFY( suppFile.is_open( ) );
    beam::fem::SupportsFromFile suppReader( suppFile );
    suppFile.close();
    mesh.setLinksOfSupports( suppReader );

    // body force
    BodyForce bodyForce( &gravity );

    // volume of ring
    Integrand volumeIntegrand( &ElementType::volumeIntegrand );
    corlib::Integrator< QuadType, Integrand, ElementType > integratorVolume( quadrature, volumeIntegrand );
    mesh.iterateOverElements( integratorVolume );
    double volume = 0.0;
    mesh.iterateOverElements( boost::bind( &ElementType::addVolume, _1, &volume ) );

    std::cout << "volume=" << volume << std::endl;
    std::cout << "volume=" << beamArea*6.0*M_PI << std::endl;

    // number dofs
//    unsigned numdofs  = 0;
//    mesh.iterateOverNodes( &NodeType::numberDOFs,  std::ref( numdofs ) );
    unsigned numdofs = mesh.iterateOverNodes( corlib::dofNumberFun( &NodeType::numberDOFs, 0 ) );
    std::cout << numdofs << " degrees of freedom" << std::endl;
    // number Lagrange multipliers
    numdofs = mesh.iterateOverLinks( corlib::dofNumberFun( &LinkType::numberDOFs, numdofs ) );  
    std::cout << numdofs << " degrees of freedom" << std::endl;

    // set up an animation writer
    int check = system( "mkdir -p animation" ); FTL_VERIFY( !check );
    std::string basename = meshFile.substr( 0, meshFile.rfind( "." ) );
    corlib::VTUanim animation( "./animation/" + basename, "vtu" );
#ifdef D_HIRES
    int checkHiRes = system( "mkdir -p animation-hr" ); FTL_VERIFY( !checkHiRes );
    corlib::VTUanim animationHiRes( "./animation-hr/" + basename + "-hr", "vtu" );
#endif
    
    // write data to VTU file
    std::string vtuFile = animation.snapshotName( 0., 0 );
    writeMeshData( vtuFile, &mesh );
#ifdef D_HIRES
    std::string vtuFileHiRes = animationHiRes.snapshotName( 0., 0 );
    writeMeshDataHiRes( vtuFileHiRes, &mesh );
#endif

    // system matrix object
    SysSolve* sysmat = new SysSolve( numdofs );
    const double tolerance = 1e-12;
    const unsigned maxiter = 30;

    // load control parameters
    const double loadMax = 5;
    const double loadStepSize = 1;
    double loadFactor = loadStepSize;

    // load control loop
    unsigned loadStep = 0;
    while ( loadFactor < (loadMax+0.5*loadStepSize) ) {
        // a word to the user
        std::cout << std::endl << "loadFactor=" << loadFactor << std::endl;

        // Newton-Raphson iteration
        bool converged = false;
        for ( unsigned i=0; i<maxiter; ++i ) {
            // a word to the user
            std::cout << "Iteration number " << i << std::endl;
        
            // Clear the element matrices and nodal forces
            mesh.iterateOverElements( boost::bind( &ElementType::clearStiffness, _1 ) );
            mesh.iterateOverNodes( boost::bind( &NodeType::clearForce, _1 ) );

            // compute element stiffness matrices
            Integrand stiffness( &ElementType::stiffnessIntegrand );
            corlib::Integrator< QuadType, Integrand, ElementType > integrator( quadrature, stiffness );
            mesh.iterateOverElements( integrator );

            // compute residual forces
            Integrand internalForce( &ElementType::internalForceIntegrand );
            corlib::Integrator< QuadType, Integrand > integrator2( quadrature, internalForce );
            mesh.iterateOverElements( integrator2 );

            // apply body forces
            bodyForce.setFactor( loadFactor );
            corlib::Integrator< QuadType, BodyForce > bodyForceIntegrator( quadrature, bodyForce );
            mesh.iterateOverElements( bodyForceIntegrator );

            // assemble the element stiffness matrices to the system matrix
            mesh.iterateOverElements(
                    corlib::MatrixGiveAndAssembleFun<const ElementType, SysSolve>(
                            &ElementType::addStiffnessMatrix, &ElementType::getDofIndices,
                            &ElementType::getDofIndices,
                            sysmat ) );

            // assemble the residual forces to the system matrices rhs
            mesh.iterateOverNodes( corlib::vectorAssemblerFun( &NodeType::giveForce,
                                                               &NodeType::copyDofArray,
                                                               sysmat ) );

            // assemble link (i.e. Lagrange multiplier) matrices to the system matrix
            mesh.iterateOverLinks(
                    corlib::MatrixGiveAndAssembleFun<const LinkType, SysSolve>(
                            &LinkType::giveCouplingMatrix, &LinkType::getDofIndices,
                            &LinkType::getDofIndicesP, sysmat ) );
            mesh.iterateOverLinks(
                    corlib::MatrixGiveAndAssembleFun<const LinkType, SysSolve>(
                            &LinkType::giveCouplingMatrixT, &LinkType::getDofIndicesP,
                            &LinkType::getDofIndices, sysmat ) );
            mesh.iterateOverLinks(
                    corlib::vectorAssemblerFun(
                            &LinkType::giveMultiplierRhs,
                            &LinkType::getDofIndices, sysmat ) );
            mesh.iterateOverLinks(
                    corlib::vectorAssemblerFun(
                            &LinkType::giveLinkResiduum,
                            &LinkType::getDofIndicesP, sysmat ) );

            // collect and apply constraints
            sysmat->finishAssembly();
            double factor = ( i==0 ? 1. : 0. );
            corlib::ConstraintFunctor< NodeType, SysSolve > constraint( sysmat, factor );
            constraint = mesh.iterateOverNodes( constraint );
            constraint.applyConstraints();

            // norm of the residual
            const double normResidual = sysmat->normRhs( );
            std::cout << "Norm of residual  " << std::scientific << normResidual << std::endl;

            // solve the linear system
            sysmat->solveSystem( );

            // add nodal solutions to increment
            mesh.iterateOverNodes( corlib::distributorFun( sysmat, 
                                                           &NodeType::copyDofArray, 
                                                           &NodeType::addToIncrement ) );
            // add Lagrange solution to multiplier
            mesh.iterateOverLinks( corlib::distributorFun( sysmat, 
                                                           &LinkType::getDofIndicesP, 
                                                           &LinkType::addToIncrement ) );

            const double normDeltaU = sysmat->normSol( );
            std::cout << "Norm of delta U  " << std::scientific << normDeltaU << std::endl;
            if ( normDeltaU < tolerance ) {
                converged = true;
                break;
            }

            // clear the stiffness matrix and RHS
            sysmat->clearMatrix( );
            sysmat->clearRhs( );

        }

        // inform user
        if ( not converged ) {
            std::cerr << "WARNING: Not converged in " << maxiter << " steps. Continuing." << std::endl;
        }

        // update displacements
        mesh.iterateOverNodes( boost::bind( &NodeType::updateDisplacements, _1 ) );
        // update Lagrange multipliers
        mesh.iterateOverLinks( boost::bind( &LinkType::updateMultiplier, _1 ) );

        // write data to VTU file
        std::string vtuFile = animation.snapshotName( loadFactor, loadStep );
        writeMeshData( vtuFile, &mesh );
#ifdef D_HIRES
        std::string vtuFileHiRes = animationHiRes.snapshotName( loadFactor, loadStep );
        writeMeshDataHiRes( vtuFileHiRes, &mesh );
#endif

        // increment load factor
        loadFactor += loadStepSize;
        loadStep ++;

    }

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
    // free the memory
    mc->destroy( );
    delete sysmat;

    return 0;

}

//------------------------------------------------------------------------------
// write solution data to a file
template< typename MESH>
void writeMeshData( const std::string & vtuFile, MESH * mesh )
{
    typedef MESH MeshType;
    typedef typename MESH::Node NodeType;

    // set up the VTU writer and open a file
    std::ofstream vtu( vtuFile.c_str( ) );
    FTL_VERIFY( vtu.is_open( ) );

    corlib::VTUwriter< MeshType > vtuwriter( mesh, vtu );

    // write mesh data
    vtuwriter.writeMesh( );
    
    //-- write point data
    vtuwriter.openPointData( );
    vtuwriter.writeNodalQuantity( corlib::accessorFun( &NodeType::giveDisplacements, "Displacements" ) );
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
    // types
    typedef MESH MeshType;

    // set up the VTU writer and open a file
    std::ofstream vtu( vtuFile.c_str( ) );
    FTL_VERIFY( vtu.is_open( ) );

    const unsigned highRes = 10;

    beam::fem::HighResVTU< MeshType > vtuwriter( mesh, vtu, highRes );

    // write mesh data
    vtuwriter.writeMesh( );

    vtuwriter.openPointData( );
    vtuwriter.writePointQuantity( corlib::accessorFun( &MESH::Node::giveDisplacements, "Displacements" ) );
    vtuwriter.writePointQuantity( corlib::accessorFun( &MESH::Node::giveVelocities   , "Velocities"    ) );
    vtuwriter.writePointQuantity( corlib::accessorFun( &MESH::Node::giveAccelerations, "Accelerations" ) );

    vtuwriter.closePointData();

    // finish writing
    vtuwriter.finalize( );
    vtu.close( );

    return;

}
