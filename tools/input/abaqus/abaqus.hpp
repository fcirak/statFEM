//
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// 
//                            Fehmi Cirak, Thomas Rueberg
//                             University of Cambridge
//                           (C) 2008 All Rights Reserved
//
// <LicenseText>
//
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//
#ifndef tools_input_abaqus_h
#define tools_input_abaqus_h

#include <stdio.h>
#include <iostream>
#include <cstdlib>
#include <vector>
#include <iterator>
#include <limits>
#include <cassert>
#include <sstream>
#include <corlib/Shape.hpp>

//------------------------------------------------------------------------------
enum abaquselement{
    B21,        //!< bar with two nodes
    B22,        //!< bar with three nodes
    S3,         //!< triangle with three nodes
    CPE4R,      //!< quadrilateral with four nodes
    C3D4,       //!< tetrahedron with four nodes
    C3D10,      //!< tetrahedron with ten nodes
    C3D8R,      //!< hexahedron  with eight nodes
    C3D20R      //!< hexahedron  with twenty nodes
};

//------------------------------------------------------------------------------
std::istream & operator>>( std::istream & inp, abaquselement & ae )
{
    std::string buffer;
    getline( inp, buffer, '\n' );
    // stream is expected to begin with "TYPE="
    std::size_t found = buffer.find( "=" );
    found++;

    if(      buffer.compare( found, buffer.size(), "B21"    ) == 0 ) { ae = B21;    return inp; }
    else if( buffer.compare( found, buffer.size(), "B22"    ) == 0 ) { ae = B22;    return inp; }
    else if( buffer.compare( found, buffer.size(), "S3"     ) == 0 ) { ae = S3;     return inp; }
    else if( buffer.compare( found, buffer.size(), "CPE4R"  ) == 0 ) { ae = CPE4R;  return inp; }
    else if( buffer.compare( found, buffer.size(), "C3D4"   ) == 0 ) { ae = C3D4;   return inp; }
    else if( buffer.compare( found, buffer.size(), "C3D10"  ) == 0 ) { ae = C3D10;  return inp; }
    else if( buffer.compare( found, buffer.size(), "C3D8R"  ) == 0 ) { ae = C3D8R;  return inp; }
    else if( buffer.compare( found, buffer.size(), "C3D20R" ) == 0 ) { ae = C3D20R; return inp; }
    else {
        std::cout << "Do not understand" << buffer << std::endl;
        exit( 1 );
    }

    return inp;
}


//------------------------------------------------------------------------------
unsigned numAbaqusPoints( const abaquselement & ae )
{
    switch( ae ) {
    case B21:    return  2; break;
    case B22:    return  3; break;
    case S3:     return  3; break;
    case CPE4R:  return  4; break;
    case C3D4:   return  4; break;
    case C3D10:  return 10; break;
    case C3D8R:  return  8; break;
    case C3D20R: return 20; break;
    default: return 0;
    }
    return 0;
}

//------------------------------------------------------------------------------
enum corlib::shape abaqusShape( const enum abaquselement & ae )
{
    enum corlib::shape es = corlib::UNDEFINED;
    switch( ae ) {
    case B21:    es = corlib::LINE;          break;
    case B22:    es = corlib::LINE;          break;
    case S3:     es = corlib::TRIANGLE;      break;
    case CPE4R:  es = corlib::QUADRILATERAL; break;
    case C3D4:   es = corlib::TETRAHEDRON;   break;
    case C3D10:  es = corlib::TETRAHEDRON;   break;
    case C3D8R:  es = corlib::HEXAHEDRON;    break;
    case C3D20R: es = corlib::HEXAHEDRON;    break;
    default:     es = corlib::UNDEFINED;
    }
    return es;
}



//------------------------------------------------------------------------------
std::string readTag( std::istream & inp )
{
    inp.ignore( ); // skip the '*' character
    std::string tagname;
    inp >> tagname;
    return tagname;
}

//------------------------------------------------------------------------------
void quitTagNotFound( std::string & tagname )
{
    std::cout<< tagname << " not found where expected! " << std::endl;
    exit( 1 );
    return;
}

//------------------------------------------------------------------------------
std::istream & readNodes( std::istream & inp, 
                          std::vector< std::vector< double > > & theNodes )
{
    // coordinate vector
    typedef std::vector< double >  coordinates;

    // read until next appearance of '*'
    unsigned nodeCtr = 0;
    while( inp.peek() != '*' ) {
        // eat number and comma
        unsigned dummy;
        inp >> dummy;
        if( inp.peek() == ',' ) inp.ignore( );
        
        // new coordinates
        coordinates aNode( 3 );
        
        // read all three components
        for ( unsigned d = 0; d < 3; d ++ ) {
            inp >> aNode[ d ];
            if( inp.peek() == ',' ) inp.ignore( );
        }
        theNodes.push_back( aNode );

        // ignore rest of line
        inp.ignore( std::numeric_limits<std::streamsize>::max() , '\n' );
        nodeCtr++;

    }
    std::cout << "Have read " << theNodes.size( ) << " nodes " << std::endl;

    return inp;
}    

//------------------------------------------------------------------------------
std::istream & readElements( std::istream & inp,
                             enum abaquselement & ae,
                             std::vector< std::vector< unsigned > > & connectivity )
{
    // get element type -> number of points
//    abaquselement ae;
    inp >> ae;
    const unsigned       n_ips = numAbaqusPoints( ae );
    //inp.ignore( std::numeric_limits<std::streamsize>::max() , '\n' );

    while( inp.peek() != '*' ) {
        unsigned dummy;
        inp >> dummy;
        if( inp.peek() == ',' ) inp.ignore();
        
        std::vector< unsigned > thisElement;
        for ( unsigned v =0; v < n_ips; v ++ ) {
            unsigned vertex;
            inp >> vertex;
            if( inp.peek() == ',' ) inp.ignore();
            thisElement.push_back( vertex-1 );
        }
        connectivity.push_back( thisElement );

        inp.ignore( std::numeric_limits<std::streamsize>::max() , '\n' );
    }
    std::cout << "Have read " << connectivity.size() << " elements" << std::endl;

    return inp;
}


//------------------------------------------------------------------------------
std::ostream & operator<<( std::ostream & out, std::vector< double > & vec )
{
    std::copy( vec.begin(), vec.end(), 
               std::ostream_iterator<double>( out, "  " ) );
    return out;
}

//------------------------------------------------------------------------------
std::istream & readNSET( std::istream & abq, 
                         std::vector< std::string > & constraints )
{
    // read the name of the NSET
    if( abq.peek() == ',' ) abq.ignore();
    std::string buffer;
    abq >> buffer;

    //------------------------------------------------------------
    std::cout << "Found " << buffer << std::endl;
    std::cout << "Give direction of constraint (x=0|y=1|z=2|all=3)" << std::endl;
    unsigned dir;
    std::cin >> dir; 
    assert( dir < 4 );
    
    char direction[] = {'x', 'y', 'z' };
    std::vector<std::pair<unsigned, double> > constraintValues;
    unsigned firstComponent, lastComponent;
    if ( dir == 0 ) {
        firstComponent = 0;
        lastComponent  = 1;
    }
    if ( dir == 1 ) {
        firstComponent = 1;
        lastComponent  = 2;
    }
    if ( dir == 2 ) {
        firstComponent = 2;
        lastComponent  = 3;
    }
    if ( dir == 3 ) {
        firstComponent = 0;
        lastComponent  = 3;
    }
    for ( unsigned c = firstComponent; c < lastComponent; c ++ ) {
        std::cout << "Give value for " << direction[c] << "-direction: "
                  << std::endl;
        double value;
        std::cin >> value;
        constraintValues.push_back( std::make_pair( c, value ) );
    }
    abq.ignore( std::numeric_limits<std::streamsize>::max() , '\n' );

    //------------------------------------------------------------
    // read constraint nodes
    std::vector< unsigned > constraintNodes;
    while( true ) {
        if( abq.peek() == '*' || abq.peek() == EOF ) break;
        std::string aux;
        getline( abq, aux, '\n' );
        std::stringstream pipe;
        pipe << aux << std::endl;
        while( true ) {
            while( pipe.peek() == ' ' ) pipe.ignore();
            if( pipe.peek() == ',' ) pipe.ignore();
            if( pipe.peek() == EOF || pipe.peek() =='\n' ) break;
            unsigned nodeNo;
            pipe >> nodeNo; 
            if ( !pipe.fail( ) ) constraintNodes.push_back( nodeNo-1 );
        }
        pipe.str( "" );
    }


    // push back the line of constraint (#Node #Direction Value)
    for( unsigned c = 0; c < constraintNodes.size(); c ++ ) {
        for ( unsigned d = 0; d < constraintValues.size(); d ++ ) {
            std::stringstream cstr;
            cstr << constraintNodes[c] << " " 
                 << constraintValues[d].first << " "
                 << constraintValues[d].second << std::endl;
            constraints.push_back( cstr.str() );
        }
    }

    return abq;
}

#endif
