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

//------------------------------------------------------------------------------
// system includes
#include <iostream>
#include <fstream>
#include <string>
#include <cstdlib>

// boost includes
#include <boost/function.hpp>

// corlib includes
#include <corlib/Mesh.hpp>
#include <corlib/Shapefun.hpp>
#include <corlib/Quadrature.hpp>
#include <corlib/SystemSolve.hpp>
#include <corlib/Constraints.hpp>
#include <corlib/Accessor.hpp>
#include <corlib/Distributor.hpp>
#include <corlib/DofNumberer.hpp>
#include <corlib/MeshSizeInfo.hpp>
#include <corlib/Integrator.hpp>
#include <corlib/BodyForce.hpp>
#include <corlib/VectorAssembler.hpp>
#include <corlib/NewmarkAssembler.hpp>
#include <corlib/BackwardEulerAssembler.hpp>
#include <corlib/ComputeElementMatrix.hpp> 
#include <corlib/PropertiesParser.hpp>
#include <corlib/Sensor.hpp>
#include <corlib/VTUwriter.hpp>
#include <corlib/NodeBasic.hpp>
#include <corlib/ElementBasic.hpp>

// solid::material
#include <solid/material/MaterialContainer.hpp>
#include <solid/material/MaterialBase.hpp>

// solid::fem includes
#include <solid/fem/NodeDynamic.hpp>
#include <solid/fem/ElementDynamic.hpp>

// local includes
#include "Loads.hpp"

//------------------------------------------------------------------------------
//!helper functions
template< typename MESH >
void writeMeshData( const std::string & vtuFile, MESH * mesh );

//------------------------------------------------------------------------------
int main( int argc, const char* argv[] )
{ 
    // Variable attributes
    const corlib::shape elemShape  = corlib::TRIANGLE;
    const unsigned      numNodesPE = 3;
    const unsigned      numGaussP  = 3;

    // Fixed attributes (some derived from the above)
    const unsigned dim = corlib::ShapeTraits<elemShape>::dim;
    typedef corlib::Shapefun<elemShape,numNodesPE>               Sfun;
    typedef corlib::Quadrature<elemShape,numGaussP>              Quad;
    typedef corlib::NodeBasic<dim>                               BasisNode;
    typedef solid::fem::NodeDynamic<BasisNode>                   Node;
    typedef solid::material::MaterialBase                        Material;
    typedef corlib::ElementBasic<Node,Sfun>                      BasisElement;
    typedef solid::fem::ElementDynamic<BasisElement,Material>    Element;
    typedef corlib::Mesh<Element>                                Mesh;
    typedef boost::function<void( Element *, 
                                  const Node::VecDim &, 
                                  const double & ) >             Integrand;
    typedef corlib::SystemSolve                                  SysSolve;


    // Variables to be read from input file
    std::string meshFile, materialFile, constraintsFile;
    double newmarkBeta, newmarkGamma, rayleigh1, rayleigh2;
    unsigned numberOfSteps;
    bool doEulerBackward;

    // Feed properties parser with the variables to read
    corlib::PropertiesParser *prop = new corlib::PropertiesParser;
    prop -> registerPropertiesVar( "meshFile",        meshFile     );
    prop -> registerPropertiesVar( "materialFile",    materialFile );
    prop -> registerPropertiesVar( "constraintsFile", constraintsFile );
    prop -> registerPropertiesVar( "newmarkBeta",     newmarkBeta   );
    prop -> registerPropertiesVar( "newmarkGamma",    newmarkGamma  );
    prop -> registerPropertiesVar( "numberOfSteps",   numberOfSteps );
    prop -> registerPropertiesVar( "rayleigh1",       rayleigh1 );
    prop -> registerPropertiesVar( "rayleigh2",       rayleigh2 );
    prop -> registerPropertiesVar( "doEulerBackward", doEulerBackward );

    // Read variables from the input.dat file
    std::ifstream inputFile( "./input.dat" );
    FTL_VERIFY( inputFile.is_open() );
    prop -> readValues( inputFile );
    delete prop;
    inputFile.close( );

    // Read material
    solid::material::MaterialContainer * mc = 
        solid::material::MaterialContainer::instance( );
    std::ifstream mf( materialFile.c_str( ) );
    FTL_VERIFY( mf.is_open( ) );
    mc -> readMaterialStream( mf );
    mf.close( );

    // Read and construct mesh
    std::ifstream smf( meshFile.c_str( ) );
    FTL_VERIFY( smf.is_open() );
    Mesh mesh( smf );
    smf.close( );

    // Obtain the 1st material
    Material * material = mc -> getMaterial( 0 );

    // cache material to elements
    mesh.iterateOverElements( boost::bind( &Element::cacheMaterial, _1, material ) );

    // mesh size info
    corlib::MeshSizeInfo<Element> meshSizeInfo;
    meshSizeInfo = mesh.iterateOverElements( meshSizeInfo );

    // compute time step according to CFL
    const double dX  = meshSizeInfo.giveAverage();
    const double c1  = material->estimatedWaveVelocity();
    const double CFL = 0.2;
    const double dt  = dX / c1 * CFL;
    std::cout << "Computing " << numberOfSteps << "  steps of size " << dt << std::endl;

    // read the constraints
    std::ifstream cstrFile( constraintsFile.c_str( ) );
    FTL_VERIFY( cstrFile.is_open( ) );
    corlib::ConstraintsFromFile<Node>  constraintsHandler( cstrFile );
    cstrFile.close();
    mesh.iterateOverNodes( constraintsHandler );
  
    // number dofs
    const unsigned numdofs = mesh.iterateOverNodes( corlib::dofNumberFun( &Node::numberDOFs, 0 ) );
    std :: cout << numdofs << " degrees of freedom" << std :: endl;

    // system matrix object
    SysSolve * sysmat  = new SysSolve( numdofs );

    // the quadrature
    Quad quadrature;

    // mass and stiffness matrix computation
    typedef boost::function<void( const Element *, 
                                  const Element::VecLDim &, 
                                  const double &, 
                                  boost::numeric::ublas::matrix<double> & ) > Integrand2;
    typedef corlib::ComputeElementMatrix<Quad,Integrand2>         MatrixComputer;
    
    Integrand2 massIntegrand(      &Element::massIntegrand );
    Integrand2 stiffnessIntegrand( &Element::stiffnessIntegrand );
    
    MatrixComputer mass(  quadrature, massIntegrand );
    MatrixComputer stiff( quadrature, stiffnessIntegrand );

    //--------------------------------------------------------------------------
    //! @name Possible time integration routines
    //@{

    //! Newmark
    corlib::NewmarkMatrixAssembler<SysSolve, Element, MatrixComputer, MatrixComputer>
        newmarkMatrix( newmarkBeta, newmarkGamma, dt, mass, stiff, rayleigh1, rayleigh2 );

    corlib::NewmarkInertiaForceAssembler<SysSolve,Element,MatrixComputer,MatrixComputer>
        newmarkForce(  newmarkBeta, newmarkGamma, dt, mass, stiff, rayleigh1, rayleigh2 );

    corlib::NewmarkSolutionDistributor<Node> 
        newmarkUpdater( newmarkBeta, newmarkGamma, dt  );

    //! Backward Euler
    corlib::BackwardEulerMatrixAssembler<SysSolve, Element, MatrixComputer, MatrixComputer>
        ebMatrix( dt, mass, stiff, rayleigh1, rayleigh2 );

    corlib::BackwardEulerInertiaForceAssembler<SysSolve,Element,MatrixComputer,MatrixComputer>
        ebForce( dt, mass, stiff, rayleigh1, rayleigh2 );

    corlib::BackwardEulerSolutionDistributor<Node>
        ebUpdater( dt );
    //@}


    // body force
    typedef boost::function< Node::VecDof( const Node::VecDim ) > VecDim2VecDof;
    typedef corlib::BodyForce<Element,VecDim2VecDof>                  BodyForce;
    BodyForce bodyForce( boost::bind( gravity, _1 ) );
    corlib::Integrator<Quad,BodyForce>  bodyForceIntegrator( quadrature, bodyForce );

    // Set up a sensor listing on displacements
    typedef boost::function<Node::VecDof( const Node* )> Node2VecDof;
    Node2VecDof sensorFun = boost::bind( &Node::giveDisplacements, _1 );
    corlib::Accessor<Node2VecDof> accessor( sensorFun );
    corlib::Sensor<corlib::Accessor<Node2VecDof> > nodeSensor( mesh.getNodePointer(1), accessor );

    // set up an animation writer
    FTL_VERIFY( !system( "mkdir -p animation" ) );
    std::string basename = "./animation/" + meshFile;
    corlib::VTUanim animation( basename, "vtu" );

    // write data to VTU file
    std::string vtuFile = animation.snapshotName( 0., 0 );
    writeMeshData( vtuFile, &mesh );

    //--------------------------------------------------------------------------
    // do time marching
    for ( unsigned n = 0; n < numberOfSteps ; n ++ ) {
        // a word to the user
        std::cout << std::endl << "time step=" << n << " at time=" << dt*n << std::endl;

        // Newton-Raphson iteration
        const double tolerance = 1e-12;
        const unsigned maxiter = 20;
        bool converged = false;
        for ( unsigned i = 0; i < maxiter; ++i ) {

            // assembly of dynamic stiffness matrix
            sysmat -> clearMatrix( );

            if ( doEulerBackward )
                mesh.iterateOverElements( boost::bind( ebMatrix, _1, sysmat ) );
            else 
                mesh.iterateOverElements( boost::bind( newmarkMatrix, _1, sysmat ) );

            // Clear element nodal forces
            mesh.iterateOverNodes( boost::bind( &Node::clearForce, _1 ) );

            // compute internal forces
            Integrand internalForce( &Element::internalForceIntegrand );
            corlib::Integrator<Quad,Integrand> integrator( quadrature, internalForce );
            mesh.iterateOverElements( integrator );

            // apply body forces (external forces)
            bodyForce.setFactor( 1.0 );
            mesh.iterateOverElements( bodyForceIntegrator );

            // assemble the internal (strain) forces to the system matrices rhs
            sysmat -> clearRhs( );
            mesh.iterateOverNodes( corlib::vectorAssemblerFun( &Node::giveForce, 
                                                               &Node::copyDofArray, 
                                                               sysmat ) );

            // add inertia forces
            if ( doEulerBackward )
                mesh.iterateOverElements( boost::bind( ebForce, _1, sysmat ) );
            else 
                mesh.iterateOverElements( boost::bind( newmarkForce, _1, sysmat ) );
            
            // collect and apply constraints
            sysmat -> finishAssembly();

            // apply constraints
            const double constrFactor = ( n==0 ) ? 1. : 0.;      // user's choice
            const double factor = ( i==0 ) ? constrFactor : 0.;  // for Newton-Raphson iteration
            corlib::ConstraintFunctor< Node, SysSolve > constraint( sysmat, factor );
            constraint = mesh.iterateOverNodes( constraint );
            constraint.applyConstraints();

            // norm of the residual
            const double normResidual = sysmat -> normRhs( );
            std::cout << "Norm of residual  " << std::scientific << normResidual << std::endl;
            //if ( normResidual < tolerance ) {
            //    converged = true;
            //    break;
            //}

            // solve the system
            sysmat -> solveSystem( );

            // add nodal solutions to increment
            mesh.iterateOverNodes( corlib::distributorFun( sysmat, 
                                                           &Node::copyDofArray,
                                                           &Node::addToIncrement ) );

            // check norm of residual displacement to determine convergence
            const double normDeltaU = sysmat->normSol( );
            std::cout << "Norm of delta U  " << std::scientific << normDeltaU << std::endl;
            if ( normDeltaU < tolerance ) {
                converged = true;
                break;
            }

        }

        // inform user
        if ( not converged ) {
            std::cout << "WARNING: Not converged in " << maxiter << " steps. Continuing." << std::endl;
        }

        // collect nodal solutions and update
        if ( doEulerBackward ) 
            mesh.iterateOverNodes( ebUpdater );
        else
            mesh.iterateOverNodes( newmarkUpdater );

        // add increments onto nodal displacements
        mesh.iterateOverNodes( boost::bind( &Node::updateDisplacements, _1 ) );

        // let sensor store time and 0-th component
        nodeSensor.listenComponent( n * dt, 0);
     
        // write data to VTU file
        std::string vtuFile = animation.snapshotName( (n+1) * dt, n+1 );
        writeMeshData( vtuFile, &mesh );

    }

    //------------------------------------------------------------
    // write animation file
    std::string animationFile( "animation.pvd" );
    std::ofstream anim( animationFile.c_str( ) );
    animation.writeAnimationFile( anim );


    //------------------------------------------------------------
    // let the sensor write its data
    const std::string gpFile = meshFile.substr( 0, meshFile.rfind( "." ) ) + ".dat";
    std::ofstream gp( gpFile.c_str() );
    FTL_VERIFY( gp.is_open( ) );
    nodeSensor.write( gp );
    gp.close( );
    
    // free the memory
    mc -> destroy( );
    delete sysmat;

    return 0;
}

//------------------------------------------------------------------------------
template< typename MESH >
void writeMeshData( const std::string & vtuFile, MESH * mesh )
{
    // set up the VTU writer and open a file
    std::ofstream vtu( vtuFile.c_str( ) );
    FTL_VERIFY( vtu.is_open( ) );
    corlib::VTUwriter<MESH> vtuwriter( mesh, vtu );


    // write mesh data
    vtuwriter.writeMesh();
  
    // write point data
    vtuwriter.openPointData();
    vtuwriter.writeNodalQuantity( corlib::accessorFun( &MESH::Node::giveDisplacements, "Displacements" ) );
    vtuwriter.writeNodalQuantity( corlib::accessorFun( &MESH::Node::giveVelocities   , "Velocities" ) );
    vtuwriter.closePointData( );

    // write element data
    vtuwriter.openCellData( );
    typedef boost::function< typename MESH::Element::Mat3x3( const typename MESH::Element * ) > Mat3x3Fun;
    Mat3x3Fun giveStress = boost::bind( &MESH::Element::firstPiolaKirchhoff, _1,
                                        corlib::ShapeTraits< MESH::Element::myShape >::centroid( ) );
    corlib::Accessor< Mat3x3Fun > accessStress( giveStress, "Stress" );
    vtuwriter.writeElementQuantity( accessStress );
    vtuwriter.closeCellData( );

    // finish writing
    vtuwriter.finalize();
    vtu.close( );
}
