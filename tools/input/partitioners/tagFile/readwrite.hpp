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

#ifndef tools_input_partitioners_tagfile_readwrite_h
#define tools_input_partitioners_tagfile_readwrite_h

#include <iostream>
#include <string>
#include <corlib/verify.hpp>


//==============================================================================
namespace tools{
    namespace input{
        namespace partitioners{
            namespace tagfile{

                template< typename ISTREAM, typename MAPTAGTOVERTEXIDS >
                void readJoinedTagFile( ISTREAM & tg,
                                        MAPTAGTOVERTEXIDS & tagToIds );
            
                template< typename MAPTAGTOVERTEXIDS, typename OSTREAM >    
                void writeJoinedTagNodeLists( const MAPTAGTOVERTEXIDS & tagToIds,
                                              std::vector<OSTREAM *> & tagLists );

                template< typename ISTREAM, typename MAPVERTEXIDTOTAG >
                void readPartitionedTagNodeLists( const std::vector<unsigned> & tags,
                                                  const std::vector<ISTREAM *> & tagLists,
                                                  MAPVERTEXIDTOTAG & idToTag );

                template< typename MAPVERTEXIDTOTAG, typename OSTREAM >
                void writePartitionedTagFile( const MAPVERTEXIDTOTAG & idToTag,
                                              OSTREAM & ptg );

            }
        }
    }
}

//==============================================================================

//------------------------------------------------------------------------------
/// Takes a sequential (united) tag file (TG file) and stores
/// its content in map: tag -> vertexIDs
///
/// \param[in]  tg         Tag file stream (TG file)
/// \param[out] tagToIds   Map: tag -> vertexIDs
template< typename ISTREAM, typename MAPTAGTOVERTEXIDS >
void tools::input::partitioners::tagfile::readJoinedTagFile( ISTREAM & tg,
                                                             MAPTAGTOVERTEXIDS & tagToIds )
{
    // read header
    std::string key;
    unsigned numTaggedVertices, numTaggedEdges;
    tg >> key >> numTaggedVertices >> numTaggedEdges;
    FTL_VERIFY( key == "TG" );

    // inform user about ignorance of edge tags
    if ( numTaggedEdges > 0 )
        std::cout << "WARNING: Edge tags are ignored and lost after translation." << std::endl;

    // read vertex tags
    typedef typename MAPTAGTOVERTEXIDS::key_type    Tag;
    typedef typename MAPTAGTOVERTEXIDS::mapped_type IdVec;
    typedef typename IdVec::value_type              Id;
    for ( unsigned v = 0; v < numTaggedVertices; ++v ) {
        Id  vId;   tg >> vId;
        Tag vTag;  tg >> vTag;
        if ( tagToIds.find( vTag ) == tagToIds.end() )
            tagToIds.insert( std::make_pair( vTag, IdVec() ) );
        tagToIds.find( vTag )->second.push_back( vId );
    }

}

//------------------------------------------------------------------------------
/// Takes a map: tag -> vertexIDs and stores each map entry in a node-list stream,
/// i.e. each tag gets a stream containing its vertexIDs
///
/// \param[in]  tagToIds   Map: tag -> vertexIDs
/// \param[out] tagLists   Node-list streams for each tag of map
template< typename MAPTAGTOVERTEXIDS, typename OSTREAM >
void tools::input::partitioners::tagfile::writeJoinedTagNodeLists( const MAPTAGTOVERTEXIDS & tagToIds, 
                                                                   std::vector<OSTREAM *> & tagLists )
{
    typedef typename MAPTAGTOVERTEXIDS::key_type    Tag;
    typedef typename MAPTAGTOVERTEXIDS::mapped_type IdVec;
    typedef typename IdVec::value_type              Id;

    // write files names for node list containing vertex tags
    typename MAPTAGTOVERTEXIDS::const_iterator tIt = tagToIds.begin(); 
    for ( ; tIt != tagToIds.end(); ++tIt ) {
        //const Tag     vTag           = tIt->first;
        const IdVec & taggedVertices = tIt->second;
        OSTREAM & tgl = *tagLists[ std::distance( tagToIds.begin(), tIt ) ];
        
        // write header, ie number of tagged vertices with tag #vTag
        tgl << taggedVertices.size() << std::endl;

        // write the vertex IDs
        typename IdVec::const_iterator vIt = taggedVertices.begin();
        for ( ; vIt != taggedVertices.end(); ++vIt ) {
            tgl << *vIt << std::endl;
        }
    }
}

//------------------------------------------------------------------------------
/// Read tags and their node-lists and store them in a map: vertexID -> tags
/// 
/// \param[in]  tags       The tags
/// \param[in]  tagLists   Node-list streams of the tags (these may be partitioned)
/// \param[out] idToTag    Map: global vertexID -> tags
template< typename ISTREAM, typename MAPVERTEXIDTOTAG >
void tools::input::partitioners::tagfile::readPartitionedTagNodeLists( const std::vector<unsigned> & tags,
                                                                       const std::vector<ISTREAM *> & tagLists,
                                                                       MAPVERTEXIDTOTAG & idToTag )
{
    FTL_VERIFY_DESCRIPTIVE( tags.size() == tagLists.size(),
                            "with tags.size()=%d, and tagLists.size()=%d\n", 
                            tags.size(), tagLists.size() );

    for ( unsigned t = 0; t < tagLists.size(); ++t ) {
        const unsigned vTag = tags[ t ];

        ISTREAM & tgl = *(tagLists[ t ]);
        unsigned numTaggedVertices;
        tgl >> numTaggedVertices;
        
        for ( unsigned v = 0; v < numTaggedVertices; ++v ) {
            unsigned vId;
            tgl >> vId;

            FTL_VERIFY( idToTag.find( vId ) == idToTag.end() );
            idToTag.insert( std::make_pair( vId, vTag ) );
        }
    }
}

//------------------------------------------------------------------------------
/// Take map: vertexID -> tags and write them to shared stream (PTG file)
/// containing pairs of _global_ vertexID and its tag after each other
///
/// \param[in]   idToTag    Map: global vertexID -> tags
/// \param[out]  ptg        Tag file stream (PTG file, partitioned due to _global_ IDs)
template< typename MAPVERTEXIDTOTAG, typename OSTREAM >
void tools::input::partitioners::tagfile::writePartitionedTagFile( const MAPVERTEXIDTOTAG & idToTag,
                                                                   OSTREAM & ptg )
{
    // write header
    ptg << "TG" << std::endl;
    ptg << idToTag.size() << " 0" << std::endl;

    // write vertex IDs with tag
    typename MAPVERTEXIDTOTAG::const_iterator it = idToTag.begin();
    for ( ; it != idToTag.end(); ++it ) {
        ptg << it->first << " " << it->second << std::endl;  // vId vTag
    }
}

#endif
