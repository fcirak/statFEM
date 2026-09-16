//
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//
//                 Computational Structural Mechanics Lab
//                         University of Cambridge
//
// This software is copyrighted and all rights are retained by the CSMLab. All
// use, disclosure, and/or reproduction of any part not expressly authorized by
// F Cirak is prohibited. (C) 2010.
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//

//! @author Burkhard Bornemann
//! @date   02/2012

/** \brief Merger of node-list files (TGL*) to subdivision surface tag file (PTG file)
 *         files used in parallel input.
 *
 *  It is the 'inverse' operation to the methods used by tg2tgl.cpp. The
 *  only difference is that the PTG files contain _global_ nodes IDs rather
 *  _local_ nodes IDs in TG. A detailed description of the file formats can be
 *  found at tg2tgl.cpp.
 */
#include <iostream>
#include <fstream>
#include <string>
#include <map>
#include <vector>
#include <algorithm>
#include <iterator>
#include <boost/lexical_cast.hpp>
#include <boost/bind.hpp>

#include <corlib/misc.hpp>
#include <corlib/verify.hpp>

#include <tools/input/partitioners/tagFile/readwrite.hpp>

//------------------------------------------------------------------------------
namespace tool{

    /// Convert partId to string with padding zeros in front
    /// e.g. '3' is converted to '003'
    std::string paddedWithZeros( const unsigned partId,
                                 const unsigned width = 3 )
    {
        std::ostringstream ss;
        ss << std::setw( width ) << std::setfill( '0' ) << partId;
        return ss.str();
    }
}

//==============================================================================
int main( int argc, char ** argv )
{
    // a word to the user
    if ( argc < 4 ) {
        std::cerr << " Usage: tgl2ptg N basename I,J,..." << std::endl
                  << " Arguments: " << std::endl
                  << "       N        -- number of partitions" << std::endl
                  << "       basename -- will search for files basename.tgl* file with tags in them" << std::endl
                  << "       I,J,...  -- name of TGL* files with *=I,J,..." << std::endl
                  << std::endl;
        return -1;
    }

    // read command line arguments
    const unsigned numParts    = boost::lexical_cast<unsigned>( argv[1] );
    const std::string basename = argv[2];
    const std::string tagList  = argv[3];

    // list tags
    std::vector<unsigned> tags;
    std::string tag;
    for ( std::string::const_iterator c = tagList.begin();
          c != tagList.end();
          ++c ) {
        if ( *c != ',' ) {
            tag.push_back( *c );
        }
        else { 
            tags.push_back( boost::lexical_cast<unsigned>( tag ) );
            tag.clear();
        }
    }
    tags.push_back( boost::lexical_cast<unsigned>( tag ) );
    std::cout << "Available tags=";
    std::copy( tags.begin(), tags.end(),
               std::ostream_iterator<unsigned>( std::cout, " " ) );
    std::cout << std::endl;

    // determine TGL file names
    std::vector<std::vector<std::ifstream *> > tagLists( numParts );
    for ( unsigned p = 0; p < numParts; ++p ) {
        for ( unsigned t = 0; t < tags.size(); ++t ) {
            const std::string tglName = basename + 
                "." + tool::paddedWithZeros( p ) +
                ".tgl" + boost::lexical_cast<std::string>( tags[t] );
            std::cout << "Opening TGL file " << tglName << std::endl;
            std::ifstream * tgl = new std::ifstream( tglName.c_str() );
            FTL_VERIFY_DESCRIPTIVE( tgl->is_open(),
                                    "Failed to open file %s\n", tglName.c_str() );
            tagLists[p].push_back( tgl );
        }
    }

    // abbreviate local helper's namespace
    namespace tagfile = tools::input::partitioners::tagfile;

    // write multiple PTG files
    for ( unsigned p = 0; p < numParts; ++p ) {

        // read multiple TGL files
        typedef std::map<unsigned,unsigned> MapIdToTag;
        MapIdToTag idToTag;
        tagfile::readPartitionedTagNodeLists( tags, tagLists[p], idToTag );

        // close input streams
        std::for_each( tagLists[p].begin(), tagLists[p].end(),
                       boost::bind( &std::ifstream::close, _1 ) );
        std::for_each( tagLists[p].begin(), tagLists[p].end(),
                       corlib::deleteFunctor() );

        // create output stream of single TG file
        const std::string ptgName = basename + "." + tool::paddedWithZeros(p) + ".ptg";
        std::cout << "Writing PTG file:" << ptgName << std::endl;
        std::ofstream ptg( ptgName.c_str() );
        FTL_VERIFY_DESCRIPTIVE( ptg.is_open(),
                            "Failed to open file %s\n", ptgName.c_str() );

        // write PTG file
        tagfile::writePartitionedTagFile( idToTag, ptg );
        
        // close output stream
        ptg.close();
    }

    return 0;
}
