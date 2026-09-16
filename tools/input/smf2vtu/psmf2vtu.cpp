// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//
//                 Computational Structural Mechanics Lab
//                         University of Cambridge
//
// This software is copyrighted and all rights are retained by the University of
// Cambridge. All use, disclosure, and/or reproduction of any part not expressly
// authorised by the University of Cambridge is prohibited. (C) 2010
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

//! @file psmf2vtu.cpp

/** \brief Convert a set of PSMF files into VTU files.
 *  \info Convert a set of files in Parallel Simple Mesh Format (PSMF) into 
 *  VTU files. In addition a PVD file is generated with a list of all VTU files. 
 */

#include <iostream>
#include <fstream>
#include <string>

#include <boost/lexical_cast.hpp>

#include <corlib/UniqueFilename.hpp>
#include <corlib/Shape.hpp>
#include <corlib/verify.hpp>

#include "readwrite.hpp"


int main( int argc, char **argv )
{
    // a word to the user
    if ( argc != 3 ) {
        std::cerr << " Usage: psmf2vtu P basename " << std::endl
                  << " Arguments: " << std::endl
                  << "       P        -- number of subdomains " << std::endl
                  << "       basename -- will search for files basename.isub.psmf file with mesh in PSMF " << std::endl
                  << std::endl;
        return -1;
    }

    // namespace aliases
    namespace ublas   = boost::numeric::ublas;
    namespace smf2vtu = tools::input::smf2vtu;    

    // read command line arguments
    const unsigned    numSub   = boost::lexical_cast<int>( argv[1] );
    const std::string basename = argv[2];
    std::vector< std::pair< unsigned, std::string > >  vtuCollection;

    // loop over all subdomains
    for ( int iSub = 0; iSub < numSub; iSub++ ) {

        // define containers for the mesh
        std::vector< ublas::bounded_vector<double,3> > coordinates;
        std::vector< std::vector<unsigned> >           elements;
        enum corlib::shape eleShape;

        // read the mesh from PSMF file
        // unique filename generator
        corlib::UniqueFilename filenameGen( 3 );
        const std::string psmfFilename = filenameGen( basename, iSub, "psmf" );
        std::ifstream psmf( psmfFilename.c_str() );
        FTL_VERIFY_DESCRIPTIVE( psmf.is_open(),
                                "Could not open file %s\n", psmfFilename.c_str() );
        smf2vtu::readSmfFile( psmf, eleShape, coordinates, elements );
        psmf.close();
        FTL_VERIFY( not elements.empty( ) );
        FTL_VERIFY( not coordinates.empty( ) );

        // point data --- dummy remains empty
        std::vector< std::pair< std::string, ublas::matrix< double > > > pointData;
        
        // write the mesh in the VTU format
        const std::string vtuFilename = psmfFilename.substr( 0, psmfFilename.find(".psmf") ) + ".vtu";
        smf2vtu::writeVtuFile( vtuFilename, eleShape, coordinates, elements, pointData, iSub );

        // add file to collection
        vtuCollection.push_back( std::make_pair( iSub, vtuFilename ) );
    }

    // write the pvd file with collection of vtu files to handle them comfortably in Paraview
    const std::string pvdFilename = basename + ".pvd";
    smf2vtu::writePvdFile( pvdFilename, vtuCollection );

    std::cout << " Conversion O.K., file " << pvdFilename << " created." << std::endl;

    return 0;
}
