// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//
//                 Computational Structural Mechanics Lab
//                         University of Cambridge
//
// This software is copyrighted and all rights are retained by the University of
// Cambridge. All use, disclosure, and/or reproduction of any part not expressly
// authorised by the University of Cambridge is prohibited. (C) 2010
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

//! @file NewmarkDriven.cpp

// system includes
#include <iostream>
#include <fstream>
#include <string>
#include <cstdlib>
#include <boost/function.hpp>

// corlib includes
#include <corlib/Mesh.hpp>
#include <corlib/Shapefun.hpp>
#include <corlib/PropertiesParser.hpp>
#include <corlib/Accessor.hpp>
#include <corlib/Sensor.hpp>
#include <corlib/NodeBasic.hpp>
#include <corlib/ElementBasic.hpp>
#include <corlib/MeshSizeInfo.hpp>

#ifdef SUPERLU
    #include <corlib/SystemSolveLU.hpp>
#else
    #include <corlib/SystemSolve.hpp>
#endif

// solid::material
#include <solid/material/MaterialContainer.hpp>
#include <solid/material/MaterialBase.hpp>

// solid::fem includes
#include <solid/fem/NodeDynamic.hpp>
#include <solid/fem/ElementDynamic.hpp>
#include <solid/driver/Dynamic.hpp>

// local includes
#include "helpers.hpp"

//==============================================================================
// everyth
int main( int argc, const char* argv[] )
{ 
    // variable attributes
    const corlib::shape elemShape  = corlib::TRIANGLE;
    const unsigned      numNodesPE = 3;
    const unsigned      numGaussP  = 3;

    // fixed attributes (some derived from the above)
    const unsigned dim = corlib::ShapeTraits<elemShape>::dim;
    typedef corlib::Shapefun<elemShape,numNodesPE>               Sfun;
    typedef corlib::NodeBasic<dim>                               BasisNode;
    typedef solid::fem::NodeDynamic<BasisNode>                   Node;
    typedef corlib::ElementBasic<Node,Sfun>                      BasisElement;
    typedef solid::material::MaterialBase                        Material;
    typedef solid::fem::ElementDynamic<BasisElement,Material>    Element;
    typedef corlib::Mesh<Element>                                Mesh;

    // solver type for linear(ised) system
#ifdef SUPERLU
    typedef corlib::SystemSolveLU                                SysSolve;
#else
    typedef corlib::SystemSolve                                  SysSolve;
#endif

    // the driver --- always use your safety belt
    typedef solid::driver::Newmark<Element,SysSolve,
                                   numGaussP,numGaussP>          Driver;


    // Variables to be read from input file
    std::string meshFile, materialFile, constraintsFile;
    double newmarkBeta, newmarkGamma;
    unsigned numTimeSteps;

    // Feed properties parser with the variables to read
    corlib::PropertiesParser * prop = new corlib::PropertiesParser;
    prop -> registerPropertiesVar( "meshFile",        meshFile        );
    prop -> registerPropertiesVar( "materialFile",    materialFile    );
    prop -> registerPropertiesVar( "constraintsFile", constraintsFile );
    prop -> registerPropertiesVar( "newmarkBeta",     newmarkBeta     );
    prop -> registerPropertiesVar( "newmarkGamma",    newmarkGamma    );
    prop -> registerPropertiesVar( "numberOfSteps",   numTimeSteps    );

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
    Driver driver( smf, newmarkBeta, newmarkGamma );
    smf.close( );

    // Obtain the 1st material
    Material * material = mc -> getMaterial( 0 );
    // cache material to elements
    driver.setMaterial( material );

    // mesh size info
    corlib::MeshSizeInfo<Element> meshSizeInfo;
    meshSizeInfo = driver.accessMesh().iterateOverElements( meshSizeInfo );

    // compute time step according to CFL
    const double dX  = meshSizeInfo.giveAverage();
    const double c1  = material->estimatedWaveVelocity();
    const double CFL = 0.2;
    const double dt  = dX / c1 * CFL;
    std::cout << "Computing " << numTimeSteps << "  steps of size " << dt << std::endl;
    // set step size in driver accordingly
    driver.setTimeStepSize( dt );

    // read the constraints
    std::ifstream cstrFile( constraintsFile.c_str( ) );
    FTL_VERIFY( cstrFile.is_open( ) );
    driver.setNodeConstraints( cstrFile );
    cstrFile.close();
    // Set nodal constraints' factor
    Driver::DoubleToDouble constrFactorFun = 
        boost::bind( &app::scalePresribedDisplacement, _1, dt );
    driver.setNodeConstraintsFactor( constrFactorFun );

    // body force
    driver.setBodyLoad( boost::bind( &app::gravity, _1 ) );  // default factor is 1

    // Set up a sensor listing on displacements
    typedef boost::function<Node::VecDof( const Node* )> NodeToVecDof;
    NodeToVecDof sensorFun = boost::bind( &Node::giveDisplacements, _1 );
    corlib::Accessor<NodeToVecDof> accessor( sensorFun );
    corlib::Sensor<corlib::Accessor<NodeToVecDof> > 
        nodeSensor( driver.accessMesh().getNodePointer(1), accessor );

    // set up an animation writer
    FTL_VERIFY( !system( "mkdir -p animation" ) );
    const std::string basename = "./animation/" + 
        meshFile.substr( meshFile.rfind( "/" ) + 1, meshFile.size() );
    driver.allocateAnimation( basename );

    // write initial data to VTU file
    driver.writeAnimationStep( driver.getTime(), 0 );

    //--------------------------------------------------------------------------
    // do time marching
    for ( unsigned n = 0; n < numTimeSteps; n ++ ) {
        // a word to the user
        std::cout << std::endl << "time step=" << n << " at time=" << driver.getTime() << std::endl;

        // Newton-Raphson iteration
        const unsigned maxiter = 20;
        const double tolerance = 1e-12;
        driver.iterate( maxiter,
                        boost::bind( std::less<double>(), _1, tolerance ),
                        boost::bind( std::less<double>(), _1, tolerance ) );

        // add increments onto nodal displacements, velocities, etc
        driver.updateSolution();

        // let sensor store time and 0-th component
        nodeSensor.listenComponent( driver.getNextTime(), 0 );
     
        // write data to VTU file
        driver.writeAnimationStep( driver.getNextTime(), n+1 );

        // update time
        driver.updateTime();

    }

    //--------------------------------------------------------------------------
    // write animation file
    driver.writeAnimationFile( "animation.pvd" );

    // let the sensor write its data
    const std::string gpFile = meshFile.substr( 0, meshFile.rfind( "." ) ) + ".dat";
    std::ofstream gp( gpFile.c_str() );
    FTL_VERIFY( gp.is_open( ) );
    nodeSensor.write( gp );
    gp.close( );
    
    // free the memory
    mc -> destroy( );

    return 0;
}

