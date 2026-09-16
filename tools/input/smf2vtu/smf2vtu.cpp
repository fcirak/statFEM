// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//
//                 Computational Structural Mechanics Lab
//                         University of Cambridge
//
// This software is copyrighted and all rights are retained by the University of
// Cambridge. All use, disclosure, and/or reproduction of any part not expressly
// authorised by the University of Cambridge is prohibited. (C) 2010
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

//! @file smfvtu.cpp

/** \brief Convert an SMF file into a VTU file.
 * \info Converts a file in Simple Mesh Format (SMF) to a file into VTU format
 * for ParaView.
 */

#include <iostream>
#include <fstream>
#include <string>

#include <corlib/Shape.hpp>
#include <corlib/verify.hpp>

#include "readwrite.hpp"


int main( int argc, char ** argv )
{
    namespace smf2vtu = tools::input::smf2vtu;

    // usage
    if ( not ( ( argc == 2 ) or ( argc == 3) ) ) {
    	std::cerr << " Usage: smf2vtu <mesh.smf> [point-data] " << std::endl
    			  << " 		Output is written to <mesh>.vtu " << std::endl
    	          << " 		point-data  -- Point data input file name" << std::endl;

    	return -1;
    }

    const std::string smfFilename = argv[1];

    std::string pdatFilename = "";
    if (argc == 3) pdatFilename = argv[2];

    // create output file name
    const std::size_t pos = smfFilename.find_last_of( "." ) ;
    std::string vtuFilename = smfFilename.substr( 0, pos ) + ".vtu";

    // define containers for the mesh
    namespace ublas = boost::numeric::ublas;
    std::vector< ublas::bounded_vector<double,3> > coordinates;
    std::vector< std::vector<unsigned> >           connectivity;
    enum corlib::shape                             eleShape;

    // read the mesh from SMF file
    std::ifstream smf( smfFilename.c_str() );    // open smf file
    FTL_VERIFY_DESCRIPTIVE( smf.is_open(),
                            "Could not find %s\n", smfFilename.c_str() );
    smf2vtu::readSmfFile( smf, eleShape, coordinates, connectivity );
    smf.close();
    FTL_VERIFY( not connectivity.empty() );
    FTL_VERIFY( not coordinates.empty() );

    // read point data file
    std::vector< std::pair< std::string, ublas::matrix< double > > > pointData;
    if ( pdatFilename != "" )
        smf2vtu::readPointDataFile( pdatFilename, pointData );

    // write the mesh in the VTU format
    smf2vtu::writeVtuFile( vtuFilename, eleShape, coordinates, connectivity, pointData );

    std::cout << " File " << vtuFilename << " created or replaced." << std::endl;

    return 0;
}
