//
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// 
//                         Computational Structural Mechanics Lab
//                             University of Cambridge
//                           (C) 2008 All Rights Reserved
//
// <LicenseText>
//
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//
#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <cstdlib>
#include <string>
#include <cassert>
#include <limits>
#include <map>
#include <vector>
#include <iterator>
#include "SmfHead.hpp"
//------------------------------------------------------------------------------
template< typename T>
void extractValue( const std::string & buffer, T & value )
{
    std::size_t equalSign = buffer.find( "=" );
    std::string rightPart = buffer.substr( equalSign+1, buffer.size() );
    std::stringstream pipe;
    pipe << rightPart;
    pipe >> value;
}

void extractName( const std::string & buffer, std::string & name )
{
    std::size_t pos1 = buffer.find( "\"" );
    std::size_t pos2 = buffer.find( "\"", pos1+1 );
    name = buffer.substr( pos1+1, pos2-3 );
}

// convert element type (string) to utils::shape enum
enum utils::shape elemTypeToShape( const std::string & et )
{
    enum utils::shape es = utils::UNDEFINED;

    if ( (et == "SPHERE") or 
         (et == "SPRING") ) {
        es = utils::POINT;
    }
    else if ( (et == "BAR") or
              (et == "BAR2") or
              (et == "BAR3") or
              (et == "BEAM") or
              (et == "BEAM2") or
              (et == "BEAM3") or
              (et == "TRUSS") or
              (et == "TRUSS2") or
              (et == "TRUSS3") ) {
        es = utils::LINE;
    }
    else if ( (et == "QUAD") or 
              (et == "QUAD4") or 
              (et == "QUAD5") or
              (et == "QUAD8") or
              (et == "QUAD9") ) {
        es = utils::QUADRILATERAL;
    }
    else if ( (et == "SHELL") or
              (et == "SHELL4") or
              (et == "SHELL8") or
              (et == "SHELL9") or 
              (et == "HEXSHELL") ) {
        es = utils::QUADRILATERAL;
    }
    else if ( (et == "TRI") or
              (et == "TRI3") or
              (et == "TRI6") or
              (et == "TRI7") ) {
        es = utils::TRIANGLE;
    }
    else if ( (et == "TRISHELL") or
              (et == "TRISHELL3") or
              (et == "TRISHELL6") or
              (et == "TRISHELL7") ) {
        es = utils::TRIANGLE;
    }
    else if ( (et == "HEX") or 
              (et == "HEX8") or
              (et == "HEX9") or
              (et == "HEX20") or
              (et == "HEX27") ) {
        es = utils::HEXAHEDRON;
    }
    else if ( (et == "TETRA") or
              (et == "TETRA4") or
              (et == "TETRA8") or
              (et == "TETRA10") or
              (et == "TETRA14") ) {
        es = utils::TETRAHEDRON;
    }
    else {
        es = utils::UNDEFINED;
    }

    return es;
}

//------------------------------------------------------------------------------
int main( int argc, const char* argv[] )
{ 
    if ( argc != 2 ) {
        std::cout << "Usage:  genesis2smf  GenesisAsciiFile " << std::endl
                  << "Note: you need to apply ncdump first!" << std::endl;
        return 0;
    }

    const std::string genFile( argv[1] );
    std::ifstream gen( genFile.c_str( ) );
    assert( gen.is_open( ) );
    const std::string baseName = genFile.substr( 0, genFile.rfind( "." ) );

    // ignore two lines
    gen.ignore( std::numeric_limits<std::streamsize>::max() , '\n' );
    gen.ignore( std::numeric_limits<std::streamsize>::max() , '\n' );

    // read the first block "dimensions:"
    std::string buffer;
    unsigned num_dim    = 0;
    unsigned num_elem   = 0;
    unsigned num_nodes  = 0;
    unsigned num_node_sets = 0;
    unsigned num_nod_per_el = 0;
    std::map< unsigned, unsigned > nsets_dim;
    while ( true ) {
        getline( gen, buffer );
        if( buffer.find( "num_dim", 0 ) != std::string::npos ) {
            extractValue( buffer, num_dim );
        }
        if( buffer.find( "num_nodes", 0 ) != std::string::npos ) {
            extractValue( buffer, num_nodes );
        }
        if( buffer.find( "num_elem", 0 ) != std::string::npos ) {
            extractValue( buffer, num_elem );
        }
        if( buffer.find( "num_node_sets", 0 ) != std::string::npos ) {
            extractValue( buffer, num_node_sets );
        }
        if( buffer.find( "num_nod_per_el", 0 ) != std::string::npos ) {
            extractValue( buffer, num_nod_per_el );
        }
        if ( num_node_sets > 0 ) {
            for ( unsigned nset = 1; nset <= num_node_sets; nset++ ) {
                std::stringstream converter;
                converter << nset;
                std::string thisNset( "num_nod_ns" );
                thisNset.append( converter.str( ) );
                if ( buffer.find( thisNset, 0 ) != std::string::npos ) {
                    unsigned num_nod_this_nset = 0;
                    extractValue( buffer, num_nod_this_nset );
                    nsets_dim[ nset ] = num_nod_this_nset;
                }
            }
        }

        if ( buffer. find( "variables:", 0 ) != std::string::npos ) break;

        if ( gen.peek() == '}' || gen.peek() == EOF ) break;
    }

    // read second block "variables:"
    std::string elem_type = "undefined";
    while ( true ) {
        getline( gen, buffer );
        if( buffer.find( "elem_type", 0 ) != std::string::npos ) {
            extractValue( buffer, elem_type );
            elem_type.erase( elem_type.begin(), elem_type.begin()+1 );
            elem_type.erase( elem_type.end()-1, elem_type.end() );
        }

        if ( buffer.find( "data:", 0 ) != std::string::npos ) break;

        if ( gen.peek() == '}' || gen.peek() == EOF ) break;
    }

    // write some information
    std::cout << "num_dim       = "  << num_dim       << std::endl
              << "num_nodes     = "  << num_nodes     << std::endl
              << "num_elem      = "  << num_elem      << std::endl
              << "num_node_sets = "  << num_node_sets << std::endl
              << "elem_type     = "  << elem_type     << std::endl;

    // move until the data block
    // ===> already done in limiting search of second block

    // go through data block
    std::map< unsigned, std::string > nset_names;
    typedef std::vector< unsigned >   VecNumbers;
    std::map< unsigned, VecNumbers >  node_set_container;
    std::vector< VecNumbers >         connectivity;
    std::vector< double >             coordinatesInARow;

    while( true ) {
        getline( gen, buffer );

        // retrieve the names of the node sets
        if ( buffer.find( "ns_names", 0 ) != std::string::npos ) {
            for ( unsigned nset = 1; nset <= num_node_sets; nset ++ ) {
                getline( gen, buffer );
                std::string nset_name;
                extractName( buffer, nset_name );
                nset_names[ nset ] = nset_name;
            }
        }
        
        // get the node numbers in each node set
        if ( buffer.find( "node_ns", 0 ) != std::string::npos ) {
            for ( unsigned nset = 1; nset <= num_node_sets; nset++ ) {
                std::stringstream converter;
                converter << nset;
                std::string thisNset( "node_ns" );
                thisNset.append( converter.str( ) );
                // look for specific nodes set
                if ( buffer.find( thisNset, 0 ) != std::string::npos ) {
                    // cut off the first part until the equal sign
                    std::size_t equalSign = buffer.find( "=" );
                    buffer = buffer.substr( equalSign+1, buffer.size() );
                    // vector of numbers
                    VecNumbers nodes_in_this_set;
                    unsigned valueCtr = 0;
                    while( true && valueCtr < nsets_dim[ nset ] ) {
                        if( gen.peek() == ';' || gen.peek() == EOF ) break;
                        std::stringstream pipe;
                        pipe << buffer << std::endl;
                        
                        while( true ) {
                            while( pipe.peek() == ' ' ) pipe.ignore();
                            if( pipe.peek() == ',' ) pipe.ignore();
                            if( pipe.peek() == EOF || 
                                pipe.peek() == '\n' ||
                                pipe.peek() == ';') break;
                            unsigned number;
                            pipe >> number; 
                            if( !pipe.fail() ) {
                                nodes_in_this_set.push_back( number-1 );
                                valueCtr++;
                            }
                        }
                        getline( gen, buffer );
                    }
                    assert( nsets_dim[ nset ] == nodes_in_this_set.size() );
                    node_set_container[ nset ] = nodes_in_this_set;
                }
            }
        }

        // get the element connectivity
        if( buffer.find( "connect", 0 ) != std::string::npos ) {
            for ( unsigned e = 0; e < num_elem; e ++ ) {
                getline( gen, buffer );
                VecNumbers thisElement;
                std::stringstream converter;
                converter << buffer;
                for ( unsigned v = 0; v < num_nod_per_el; v ++ ) {
                    if( converter.peek() == ',' ) converter.ignore( );
                    if( converter.peek() == '\n' ||
                        converter.peek() == ';'  ||
                        converter.peek() == EOF ) break;
                    unsigned numV;
                    converter >> numV;
                    if( !converter.fail() )  thisElement.push_back( numV-1);
                }
                connectivity.push_back( thisElement );
            }
        }
        
        // get the coordinates (note: they'll be order like:
        //                      x1 x2 ... y1 y2 ... z1 z2 .. )
        if ( buffer.find( "coord", 0 ) != std::string::npos ) {
            while( true ) {
                if( gen.peek() == '}'  ||
                    gen.peek() == EOF ) break;
                getline( gen, buffer );
                std::stringstream converter;
                converter << buffer;
                while( true ) {
                    if( converter.peek() == '\n' ||
                        converter.peek() == ';'  ||
                        converter.peek() == EOF ) break;
                    if( converter.peek() == ',' ) converter.ignore( );
                    double value;
                    converter >> value;
                    if( !converter.fail() ) coordinatesInARow.push_back( value );
                }
            }
            assert( coordinatesInARow.size() == num_dim * num_nodes );
        }
        if ( gen.peek() == '}' || gen.peek() == EOF ) break;
    }


     //----------------------------------------------------------------------
     // write smf file
     std::string smfFile( baseName + ".smf" );
     std::ofstream smf( smfFile.c_str( ) );
     // write header
     utils::SmfHead smfHead;
     smfHead.setElementShape( elemTypeToShape( elem_type ) );
     smfHead.setElementNumPoints( num_nod_per_el );
     smfHead.write( smf );
     // write node element sums
     smf << num_nodes << "  " << num_elem << std::endl;
     // write nodes
     for( unsigned n =0; n < num_nodes; n++ )  {
         smf << coordinatesInARow.at( n             ) <<  "  "
             << coordinatesInARow.at( n + num_nodes ) <<  "  ";
         if ( num_dim > 2 ) 
             smf << coordinatesInARow.at( n + 2*num_nodes ) <<  "  ";
         else smf << "0.";
         smf << std::endl;
     }
     // write elements
     for( unsigned e=0; e < connectivity.size(); e ++ ) {
         std::copy( connectivity[e].begin(),connectivity[e].end(),
                    std::ostream_iterator< unsigned >( smf, "  " ) );
         smf << std::endl;
     }
     smf.close();


     if ( num_node_sets == 0 ) return 0;

     //----------------------------------------------------------------------
     // if NSETs are given, ask for prescribed values and write the file
     std::vector< std::string > constraintsBuffer;

     for ( unsigned nset = 1; nset <= num_node_sets; nset ++ ) {
         std::cout << "Found NSET " << nset_names[ nset ] << std::endl
                   << "Give direction of constraint (x=0|y=1|z=2|all=3)" << std::endl;
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

         // write this nodes set to the constraint file
         for ( unsigned n = 0; n < nsets_dim[ nset ]; n ++ ) {
             for ( unsigned c = 0; c < constraintValues.size(); c ++ ) {
                 std::stringstream converter;
                 converter << node_set_container[ nset ][n]    << "  " 
                           << (constraintValues.at( c )).first << "  "
                           << (constraintValues.at( c )).second;
                 constraintsBuffer.push_back( converter.str( ) );
             }
         }
     }

     //----------------------------------------------------------------------
     // write the constraints
     std::string constraintFile( baseName + ".constraints" );
     std::ofstream cf( constraintFile.c_str( ) );
     cf << constraintsBuffer.size() << std:: endl;
     std::copy( constraintsBuffer.begin(), constraintsBuffer.end(),
                std::ostream_iterator< std::string >( cf, "\n" ) );

     cf.close();
    
     gen.close( );
     return 0;
}



