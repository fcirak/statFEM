// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//
//                 Computational Structural Mechanics Lab
//                         University of Cambridge
//
// This software is copyrighted and all rights are retained by the University of
// Cambridge. All use, disclosure, and/or reproduction of any part not expressly
// authorised by the University of Cambridge is prohibited. (C) 2010
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

//! @file   singleSpanDriven.cpp

//! system includes
#include <iostream>
#include <fstream>
#include <string>

//! corlib includes
#include <corlib/eigenX.hpp>
#include <corlib/SystemSolveEigenSparse.hpp>
//! beam includes
#include <beam/driver/Beam.hpp>

namespace eigenX = corlib::eigenX;

//------------------------------------------------------------------------------
//! Applied body load
eigenX::VectorSd< 2 >
forceFun( const eigenX::VectorSd< 2 > coorRef )
{
    eigenX::VectorSd< 2 > weigDens;
    weigDens.setZero();
    weigDens( 1 ) = - 100.;  // Kg/m^3 * N/Kg
    return weigDens;
}

//==============================================================================
int main( int argc, const char* argv[] )
{
    const unsigned degree = 3;
    const unsigned embDim = 2;
    const bool     dynamic = true;
    
    std::string meshFile, materialFile, supportsFile, loadsFile;
    double beamThickness, stepSize; 
    unsigned numSteps, every;
    bool isDynamic, useEuler;
    double newmarkBeta=0.25, newmarkGamma=0.5;
    double tolerance;
    unsigned maxiter;
    double rayleighM, rayleighK;

    // feed properties parser with the variables to read
    corlib::PropertiesParser *prop = new corlib::PropertiesParser;
    prop -> registerPropertiesVar( "meshFile",          meshFile   );
    prop -> registerPropertiesVar( "materialFile",      materialFile );
    prop -> registerPropertiesVar( "supportsFile",      supportsFile );
    prop -> registerPropertiesVar( "beamThickness",     beamThickness );
    prop -> registerPropertiesVar( "isDynamic",         isDynamic );
    prop -> registerPropertiesVar( "useEuler",          useEuler );
    prop -> registerPropertiesVar( "stepSize",          stepSize );
    prop -> registerPropertiesVar( "numSteps",          numSteps );
    prop -> registerPropertiesVar( "every",             every );    
    prop -> registerPropertiesVar( "newmarkBeta",       newmarkBeta );
    prop -> registerPropertiesVar( "newmarkGamma",      newmarkGamma );
    prop -> registerPropertiesVar( "tolerance",         tolerance );
    prop -> registerPropertiesVar( "maxiter",           maxiter );
    prop -> registerPropertiesVar( "rayleighM",         rayleighM );
    prop -> registerPropertiesVar( "rayleighK",         rayleighK );

    // read variables from the input.dat file
    std::ifstream inputFile( "./input.dat" );
    FTL_VERIFY( inputFile.is_open() );
    prop -> readValues( inputFile );
    delete prop;
    inputFile.close();

    // read material
    solid::material::MaterialContainer * mc = solid::material::MaterialContainer::instance();
    std::ifstream mf( materialFile.c_str() );
    FTL_VERIFY( mf.is_open() );
    mc -> readMaterialStream( mf );
    mf.close();

    // define type of driver for the beam
    typedef corlib::SystemSolveEigenSparse Solver;
    typedef beam::driver::Beam<beam::driver::Traits<embDim,degree,dynamic>::Element,Solver> Driver;
    
    // instatiate a driver
    std::ifstream smf( meshFile.c_str() );
    FTL_VERIFY( smf.is_open() );
    Driver driver( smf );
    smf.close();

    // obtain the 1st material
    solid::material::MaterialBase* material = mc -> getMaterial( 0 );

    //! cache material to elements
    driver.setMaterial( material );

    driver.setRayleighDamping( rayleighM, rayleighK );

    // set the thickness of the beam
    driver.setThickness( beamThickness );

    //! set time integrator
    if ( dynamic ) {
        if ( useEuler ) driver.setTimeInterator( Driver::TI_BACKWARDEULER );
        else            driver.setTimeInterator( Driver::TI_NEWMARK, newmarkBeta, newmarkGamma );
    }
    else driver.setTimeInterator( Driver::TI_STATIC );

    //! read links
    std::ifstream suppFile( supportsFile.c_str() );
    FTL_VERIFY( suppFile.is_open() );
    driver.prepareLinks( suppFile );
    suppFile.close();

    //! body force
    driver.setBodyForceFun( &forceFun );

    //! number dofs
    std::pair<unsigned,unsigned> aux = driver.numberDofs();
    std::cout << aux.first << " degrees of freedom and " << aux.second << " Lagrange multipliers "
              << std::endl;

    // output folder
    const std::string animationDir     = "animation/";
    const std::string makeAniamtionDir = "mkdir -p " + animationDir;
    FTL_VERIFY( !system( makeAniamtionDir.c_str() ) );

    // write data to VTU file
    const std::string outFile0 = animationDir + corlib::UniqueFilename()( "output", 0, "vtu" );
    driver.writeMeshData( outFile0 );

    //! system matrix object
    driver.allocateSolver( aux.first + aux.second );

    // time control parameters
    double time = 0;
    if ( not isDynamic ) numSteps = 1;

    const std::string monitorFile = animationDir + "nodeMonitor.dat";
    std::ofstream nm( monitorFile.c_str() );

    // time loop
    for ( unsigned step = 1; step <= numSteps; step ++ ) {

        time += stepSize;

        // a word to the user
        std::cout << std::endl << "time=" << time << std::endl;

        // Newton-Raphson iteration
        const bool converged =
            driver.iterate( stepSize, maxiter, 
                            std::bind( std::less<double>(), std::placeholders::_1, tolerance ),
                            std::bind( std::less<double>(), std::placeholders::_1, tolerance ) );


        // compute latest velocities and acclerations; clear increments
        driver.updateSolution( stepSize );

        // inform user
        if ( not converged ) {
            std::cout << "WARNING: Not converged in " << maxiter << " steps. Continuing." << std::endl;
        }

        if ( step % every == 0 ) {
            // write data to VTU file
            const std::string outFile = animationDir + corlib::UniqueFilename()( "output", step, "vtu" );
            driver.writeMeshData( outFile );
        }
        
        driver.monitorNode( 50, time, nm );
        std::cout << time << "  " << driver.linkValue( 1 ) << std::endl;
    }

    nm.close();

    return 0;
}
