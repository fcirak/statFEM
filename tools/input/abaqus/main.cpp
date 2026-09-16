//
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// 
//                         Computational Structural Mechanics Lab
//                             University of Cambridge
//                           (C) 2010 All Rights Reserved
//
// <LicenseText>
//
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//

#include "abaqus.hpp"
#include <iostream>
#include <fstream>
#include <string>
#include <cstdlib>
#include <corlib/SmfHead.hpp>

//------------------------------------------------------------------------------
int main( int argc, const char* argv[] )
{ 
    if ( argc != 2 ) {
        std::cout << "Usage:  abaqus  AbaqusInputFile " << std::endl;
        return 0;
    }

    const std::string abqFile( argv[1] );
    std::ifstream abq( abqFile.c_str( ) );
    assert( abq.is_open( ) );

    std::string tagname;
    // read heading
    if ( abq.peek( ) == '*' ) tagname = readTag( abq );
    std::string HEADING( "Heading" );
    if ( tagname.compare( 0, HEADING.size(), HEADING ) != 0 ) 
        quitTagNotFound( HEADING );
    abq.ignore( std::numeric_limits<std::streamsize>::max() , '\n' );


    // skip the following lines until '*'-char
    while( abq.peek() != '*' ) abq.ignore( );

    //----------------------------------------------------------------------
    // read nodes
    tagname = readTag( abq );
    std::string NODE( "Node" );
    if ( tagname.compare( 0, NODE.size(), NODE ) != 0 ) 
        quitTagNotFound( NODE );
    abq.ignore( std::numeric_limits<std::streamsize>::max() , '\n' );

    typedef std::vector< double >      coordinates;
    typedef std::vector< coordinates >       nodes;
    nodes theNodes;
    readNodes( abq, theNodes );


    //----------------------------------------------------------------------
    // read elements
    tagname = readTag( abq );
    std::string ELEMENT( "Element" );
    if ( tagname.compare( 0, ELEMENT.size(), ELEMENT ) != 0 ) 
        quitTagNotFound( ELEMENT );

    enum abaquselement abqele;
    std::vector< std::vector< unsigned > > connectivity;
    readElements( abq, abqele, connectivity );

    //----------------------------------------------------------------------
    // write smf file
    std::string smfFile(abqFile + ".smf" );
    std::ofstream smf( smfFile.c_str( ) );
    // write header
    corlib::SmfHead smfHead;
    smfHead.setElementShape( abaqusShape( abqele ) );
    smfHead.setElementNumPoints( numAbaqusPoints( abqele ) );
    smfHead.write( smf );
    // write node and element sum
    smf << theNodes.size( ) << "  " << connectivity.size( ) << std::endl;
    // write nodes
    for( unsigned n =0; n < theNodes.size(); n++ ) 
        smf << theNodes.at(n) << std::endl;
    // write elements
    for( unsigned e=0; e < connectivity.size(); e ++ ) {
        std::copy( connectivity[e].begin(),connectivity[e].end(),
                   std::ostream_iterator< unsigned >( smf, "  " ) );
        smf << std::endl;
    }
    smf.close();


    //----------------------------------------------------------------------
    // read nodal constraints (NOTE: ESETS are not yet supported!)
    std::vector< std::string > constraints;
    while( abq.peek() == '*' ) {

        tagname = readTag( abq );
        std::string NSET( "NSET" );
        if ( tagname.compare( 0, NSET.size(), NSET ) == 0 ) {
            readNSET( abq, constraints );
        }
        else {
            // function read ESET yet to be done
            quitTagNotFound( NSET );
        }

    }// finished reading NSETS
    

    std::string constraintFile( abqFile + ".constraints" );
    std::ofstream cf( constraintFile.c_str( ) );

    cf << constraints.size( ) << std::endl;
    std::copy( constraints.begin(), constraints.end(),
               std::ostream_iterator<std::string>( cf ) );
    cf.close();
    
    abq.close( );
    return 0;
}



