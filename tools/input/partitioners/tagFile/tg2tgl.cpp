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

/** \brief Converter of subdivision surface tag file (TG file) to vertex list
 *         files used in parallel input.
 *
 *  The converter splits the vertex tags in the TG file to separate files;
 *  the edge tags cannot be processed in the moment.
 *
 *  The layout of the TG file is described in <subdiv/surf/Mesh.hpp>. For each
 *  vertex tag a separate vertex list file TGL<int> is generated containing
 *  the total number of vertices subject to the vertex tag VTAG and their IDs.
 *
 *  Example:
 *  The TG file (basename.tg) is
 *  \verbatim
 *     TG
 *      4 0
 *     2 0
 *     7 2
 *     8 2
 *     9 2
 *  \endverbatim
 *  ie there in total 4 tagged vertices of which one has tag '0'
 *  and three (vertices 7 to 9) have tag '2'.
 *  Thus, the result of vertex tag list TGL file (basename.tgl0) is
 *  \verbatim
 *     1
 *    2
 *  \endverbatim
 *  and the second TGL file (basename.tgl2) contains
 *  \verbatim
 *     3
 *    7
 *    8
 *    9
 *  \endverbatim
 */
#include <iostream>
#include <fstream>
#include <string>
#include <map>
#include <vector>
#include <algorithm>
#include <boost/lexical_cast.hpp>
#include <boost/bind.hpp>

#include <corlib/misc.hpp>
#include <corlib/verify.hpp>

#include <tools/input/partitioners/tagFile/readwrite.hpp>

//------------------------------------------------------------------------------
int main( int argc, char ** argv )
{
    // a word to the user
    if ( argc < 2 ) {
        std::cerr << " Usage: tg2tgl basename " << std::endl
                  << " Arguments: " << std::endl
                  << "       basename -- will search for files basename.tg file with tags in them " << std::endl
                  << std::endl;
        return -1;
    }

    // read command line arguments
    const std::string basename = argv[1];

    // open tag file
    const std::string tgName = basename + ".tg";
    std::ifstream tg( tgName.c_str() );
    FTL_VERIFY_DESCRIPTIVE( tg.is_open(),
                            "Failed to open file %s\n", tgName.c_str() );

    // abbreviate local helper's namespace
    namespace tagfile = tools::input::partitioners::tagfile;
    
    // read file
    typedef std::map<unsigned,std::vector<unsigned> > MapTagToIds;
    MapTagToIds tagToIds;
    tagfile::readJoinedTagFile( tg, tagToIds );
    tg.close();
    
    // write files names for node list containing vertex tags
    std::vector<std::ofstream *> tagLists;
    MapTagToIds::const_iterator tIt = tagToIds.begin(); 
    for ( ; tIt != tagToIds.end(); ++tIt ) {
        // convenience
        const unsigned vTag = tIt->first;
        
        // make file name
        const std::string tglName = basename + ".tgl" +
            boost::lexical_cast<std::string>( vTag );
        tagLists.push_back( new std::ofstream( tglName.c_str() ) );
    }

    // write into streams
    tagfile::writeJoinedTagNodeLists( tagToIds, tagLists );

    // close streams
    std::for_each( tagLists.begin(), tagLists.end(),
                   boost::bind( &std::ofstream::close, _1 ) );
    std::for_each( tagLists.begin(), tagLists.end(),
                   corlib::deleteFunctor() );

    return 0;
}
