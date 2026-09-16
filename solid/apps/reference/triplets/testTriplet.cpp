// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//
//                 Computational Structural Mechanics Lab
//                         University of Cambridge
//
// This software is copyrighted and all rights are retained by the University of
// Cambridge. All use, disclosure, and/or reproduction of any part not expressly
// authorised by the University of Cambridge is prohibited. (C) 2010
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

//! @file testTriplet.cpp

#include <cstdlib>
#include <iostream>
#include <fstream>
#include <functional>
#include <stdlib.h>

#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_int.hpp>
#include <boost/random/variate_generator.hpp>

// object tested
#include <corlib/Triplet.hpp>

# define NROW 5000000
# define ALLOCATE 10*NROW

//------------------------------------------------------------------------------
// Simple tester for Triplet class
boost::mt19937 gen;

unsigned gen_indices() {
    boost::uniform_int<> dist(0, NROW - 1);
    boost::variate_generator<boost::mt19937&, boost::uniform_int<> > indices(gen, dist);
    return static_cast<unsigned>( indices() );
}
double gen_values() {
    boost::uniform_int<> dist( - NROW, NROW );
    boost::variate_generator<boost::mt19937&, boost::uniform_int<> > indices(gen, dist);
    return static_cast<double> ( indices() );
}

int main( int argc, char *argv[] )
{
    typedef corlib::Triplet<unsigned,double> Triplet;
    Triplet triplet;

    // prepare triplet
    triplet.prepareAssembly( ALLOCATE );


    //triplet.insert( 0, 2, 13.0);
    // insert some values 
    //triplet.insert( 0, 2, 13.0);
    
    // insert values in a loop
    for ( unsigned i = 0; i < ALLOCATE; i++ ) {
        unsigned iRow = gen_indices( );
        unsigned iCol = gen_indices( );
        double   val  = gen_values( );

        triplet.insert( iRow, iCol, val );
     }

    // finalize the assembly
    triplet.finishAssembly( );

    // zero row 1
    triplet.zeroRow( 1, -333.0);

    //std::cout << " *** final state " << std::endl;
    //triplet.write( std::cout );


    std::cout << " inf norm of the matrix " << triplet.infNorm( ) << std::endl;
    std::cout << " Frobenius norm of the matrix " << triplet.frobeniusNorm( ) << std::endl;
    std::cout << " recommended diagonal scalar for fixing BC " << triplet.getDiagScalar( ) << std::endl;
    
}

