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

//! @author Kosala Bandara
//! @date   12/2012

/** \brief Convertor to generate a tag file for a mesh
 * \info Compares the coordinates of the given mesh region with that of the full 
 * mesh to identify the tagged node indices.
 */

#include <iostream>
#include <fstream>
#include <vector>

#include <corlib/verify.hpp>
#include <corlib/fuzzyEqual.hpp>
#include <corlib/SmfHead.hpp>

#include <boost/numeric/ublas/vector.hpp>

#include "options.hpp"
#include "helpers.hpp"
//==============================================================================

int main( int argc, char **argv )
{
    namespace smf2tg = tools::helpers::smf2tg;
    namespace ublas = boost::numeric::ublas;
    // read options
    std::string meshFilename, regionFilename, outFilename;
    double tol;
    bool inverse = false, writeTag = false, writeETag = false, writeFix = false;
    const int status = smf2tg::readOptions( argc, argv, meshFilename, 
                                            regionFilename, outFilename, tol, 
                                            inverse, writeTag, writeETag, writeFix);

    if ( status != 0 ) return status;
    
    // read mesh files
    std::vector< ublas::bounded_vector<double,3> >                fullCoord;
    std::vector< ublas::bounded_vector<double,3> >                regionCoord;
    std::set<std::pair<unsigned, unsigned>, smf2tg::CompareEdges> fullEdgeSet;
    std::ifstream full( meshFilename.c_str() );
    std::ifstream region( regionFilename.c_str() );
    FTL_VERIFY_DESCRIPTIVE( full.is_open(),
                            "Could not find %s\n", meshFilename.c_str() );
    FTL_VERIFY_DESCRIPTIVE( region.is_open(),
                            "Could not find %s\n", regionFilename.c_str() );

    smf2tg::readSmfFiles( full, region, fullCoord, regionCoord, fullEdgeSet );
    full.close();
    region.close();
    FTL_VERIFY( fullCoord.size() );
    FTL_VERIFY( regionCoord.size() );

    // find tagged vertices
    // (coordiantes that appear in both fullMesh and region)
    std::set<unsigned> taggedNodes;
    smf2tg::CompareNodes compareNodes(regionCoord, tol);
    for (unsigned i = 0; i < fullCoord.size(); i++){
        ublas::bounded_vector<double,3> node = fullCoord[i];
        unsigned match = compareNodes.findNode(node);
        if ( ( (not inverse) and 
               (match < regionCoord.size() ) ) or 
             ( (inverse) and 
               (match == regionCoord.size() ) ) )
            taggedNodes.insert(i);
    }

    // find tagged edges
    // (edges with coordiantes that appear in both fullMesh and region)
    typedef std::set<std::pair<unsigned, unsigned>, 
                     smf2tg::CompareEdges>::iterator EdgeSetIter;
    std::set<std::pair<unsigned, unsigned>, smf2tg::CompareEdges> taggedEdges;
    EdgeSetIter sBegin =fullEdgeSet.begin();
    EdgeSetIter sEnd =fullEdgeSet.end();
    for (; sBegin != sEnd; sBegin++){
        const unsigned v1 = (*sBegin).first;
        const unsigned v2 = (*sBegin).second;

        unsigned match1 = compareNodes.findNode( fullCoord[v1] );
        unsigned match2 = compareNodes.findNode( fullCoord[v2] );

        if ( ( (not inverse) and 
               (match1 < regionCoord.size() ) and
               (match2 < regionCoord.size() ) ) or
             ( (inverse) and 
               (match1 == regionCoord.size() ) and
               (match2 == regionCoord.size() ) ) )
            taggedEdges.insert(std::make_pair(v1, v2));
    }

    // write pdat file
    if (taggedNodes.empty()) 
        std::cout << "No matching nodes found"<<std::endl;
    else {
        std::cout << taggedNodes.size()<<" matching nodes found"<<std::endl;
        if(writeETag) std::cout << taggedEdges.size()<<" matching edges found"<<std::endl;
        smf2tg::writeOutput(outFilename, fullCoord.size(), taggedNodes, taggedEdges, 
                            writeTag, writeETag, writeFix);
    }
    return 0;
}

