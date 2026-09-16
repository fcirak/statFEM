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

//! @author Jakub Sistek, Burkhard Bornemann
//! @date   01/2012

#ifndef tools_input_partitioners_utils_readoptions_h
#define tools_input_partitioners_utils_readoptions_h

#include <iostream>
#include <vector>
#include <string>
#include <boost/lexical_cast.hpp>

//==============================================================================
namespace tools{
    namespace input{
        namespace partitioners{
            namespace utils{

                void printUsageMeshBasedPartitioner();

                void readOptionsMeshBasedPartitioner( int argc,
                                                      char * argv[],
                                                      int & numSub,
                                                      unsigned & overlapLayers,
                                                      std::string & basename,
                                                      std::vector<std::string> & listFiles );

            }
        }
    }
}

//==============================================================================
/// Print a usage information
void tools::input::partitioners::utils::printUsageMeshBasedPartitioner()
{
    std::cerr << " SMF - Partitioner " << std::endl << std::endl
              << " Usage: " << std::endl
              << "       ./partitioner P OVL basename [nodeListFile1 nodeListFile2 ... ]" << std::endl
              << std::endl
              << " Arguments: " << std::endl
              << "       P        -- target number of subdomains " << std::endl
              << "       OVL      -- layers of overlap" << std::endl
              << "                   0 - for substructuring (no overlap) " << std::endl
              << "                   1 - e.g. for filtering values and Lagrange elements  " << std::endl
              << "                   N - e.g. sufficient width to assembly matrix rows locally " << std::endl
              << "       basename -- will search for file basename.smf with mesh in SMF " << std::endl
              << "  nodeListFileX -- will search for file basename.nodeListFileX with list of " << std::endl
              << "                   global nodes and partition it to new global numbering " << std::endl
              << std::endl;
    exit( -1 );
    return;
}

//==============================================================================
/// Extract arguments from commmand line options
///
/// \param[in]  argc            Number of command-line arguments
/// \param[in]  argv            List of command-line arguments
/// \param[out] numSub          Number of partitions (or sub-domains)
/// \param[out] overlapLayers   Depth of overlapping of partitions ('ghost' elements and nodes)
/// \param[out] basename        Base name of mesh and node list files
/// \param[out] listFiles       Node list files (optional; this may be returned empty)
void tools::input::partitioners::utils::
readOptionsMeshBasedPartitioner( int argc,
                                 char * argv[],
                                 int & numSub,
                                 unsigned & overlapLayers,
                                 std::string & basename,
                                 std::vector<std::string> & listFiles )
{
    // number of arguments
    const int numArgs = argc;

    if ( numArgs >= 4 ) { // standard call
        numSub                = boost::lexical_cast<int>( argv[1] );
        overlapLayers         = boost::lexical_cast<int>( argv[2] );
        basename              = boost::lexical_cast<std::string>( argv[3] );
        // load files with lists of global nodes to partition
        const unsigned numListFiles = numArgs - 4;
        for ( unsigned i = 0; i < numListFiles; i++ ) {
            listFiles.push_back( std::string( argv[ i + 4 ] ) );
        }
    }
    else {
        partitioners::utils::printUsageMeshBasedPartitioner();
    }
}

#endif
