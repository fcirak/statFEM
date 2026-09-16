// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//
//                 Computational Structural Mechanics Lab
//                         University of Cambridge
//
// This software is copyrighted and all rights are retained by the University of
// Cambridge. All use, disclosure, and/or reproduction of any part not expressly
// authorised by the University of Cambridge is prohibited. (C) 2010 --
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

//! @file brasilRoots.cpp

// system headers
#include <iostream>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <iomanip>
#include <limits>

// corlib headers
#include <corlib/verify.hpp>

//------------------------------------------------------------------------------
void readResult ( const std::string & fileName, std::vector<double> & out )
{
    std::ifstream fileStream( fileName );
    FTL_VERIFY( fileStream.is_open( ) );
    for( unsigned i = 0; i < out.size( ); ++ i )
        fileStream >> out[ i ];
    fileStream.close( );

    return;
}

//------------------------------------------------------------------------------
int main( )
{
    // check if baryrat is installed in user's python environment
    int status = std::system( "python3 -c 'import baryrat'" );
    FTL_VERIFY_DESCRIPTIVE( status == 0,
                            "Cannot call python3 -c 'import baryrat'." );

    // input to python script
    const unsigned m = 4;   // polynomial degree
    double exponent = 0.5;  // exponent of power function

    // command line
    const std::string tool  = "python3";
    const std::string file  = "$STATFEMROOT/tools/baryrat/baryrat-1.4.0/roots.py";
    const std::string argv1 = std::to_string( m );

    // test multiple runs
    for ( unsigned k = 0; k < 3; ++ k )
    {
        // convert double to std::string
        std::stringstream stream;
        stream << std::fixed
               << std::setprecision( std::numeric_limits<double>::digits10 + 1 )
               << exponent;

        // call python script from command line
        const std::string argv2 = stream.str( );
        const std::string line = tool + " " + file + " " + argv1 + " " + argv2;
        std::cout << line << "\n";
        std::system( line.c_str( ) );

        // output from python script
        std::vector<double> topRoots( m + 1,  0. );
        std::vector<double> botRoots( m + 1,  0. );

        // read approximation result
        readResult( "top.out", topRoots );
        readResult( "bot.out", botRoots );

        // display approximation result
        std::cout << "exponent = " << exponent << "\n";
        std::cout << "numerator:" << "\t";
        for( unsigned j = 0; j < topRoots.size( ); ++j )
            std::cout << topRoots[ j ] << " ";
        std::cout << "\n";
        std::cout << "denominator:" << "\t";
        for( unsigned j = 0; j < botRoots.size( ); ++j )
            std::cout << botRoots[ j ] << " ";
        std::cout << "\n";

        // different exponent values
        exponent += 0.1;
    }
    
    return 0;
}

