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

#ifndef subdiv_surf_io_h
#define subdiv_surf_io_h

#include <algorithm>

#include <subdiv/surf/Vertex.hpp>
#include <subdiv/surf/Facet.hpp>
#include <subdiv/surf/FacetTree.hpp>

//==============================================================================
namespace subdiv {
    namespace surf {

        std::ostream & write( std::ostream & os,
                              const Vertex * vertex );

        template< corlib::shape SHAPE, typename VERTEX >
        std::ostream & write( std::ostream & os,
                              Facet< SHAPE, VERTEX > * facet );

        template< corlib::shape SHAPE, typename VERTEX >
        void printConnectivityData( std::ostream & os,
                                    Facet< SHAPE, VERTEX > * facet );
        
        template< corlib::shape SHAPE, typename VERTEX >
        std::ostream & write( std::ostream & os,
                              FacetTree< SHAPE, VERTEX > * facet );


    }
}

//==============================================================================

//------------------------------------------------------------------------------
std::ostream & subdiv::surf::write( std::ostream & os,
                                    const Vertex * vertex )
{
    os << "Vertex " << vertex->index( ) << std::endl;
    os << "    tag : " << vertex->tag( ) << std::endl;
    return os;
}

//------------------------------------------------------------------------------
template< corlib::shape SHAPE, typename VERTEX >
std::ostream & subdiv::surf::write( std::ostream & os,
                                    Facet< SHAPE, VERTEX > * facet )
{
    typedef VERTEX                 Vertex;
    typedef Facet< SHAPE, Vertex > Facet_;

    os << "Facet " << facet->index( ) << std::endl;

    os << "    vertex IDs : ";
    std::transform( facet->verticesBegin(), facet->verticesEnd(),
                    std::ostream_iterator< unsigned >( os, " " ),
                    std::bind( &Vertex::index, std::placeholders::_1 ) );
    os << std::endl;

    os << "    neighbour IDs : ";
    corlib::transform_if( facet->neighborsBegin( ), facet->neighborsEnd( ),
                          std::ostream_iterator< unsigned >( os, " " ),
                          std::bind( &Facet_::index, std::placeholders::_1 ),
                          std::bind( std::not_equal_to< Facet_ * >( ), std::placeholders::_1,
                                       static_cast< Facet_ * >( NULL ) ) );
    os << std::endl;

    os << "    all (non-)manifold neighbour IDs : ";
    for ( unsigned n = 0; n < Facet_::numNeighbors; n++ ) {
        std::vector< Facet_ * > all;
        facet->neighborsAtEdge( n, all );
        std::transform( all.begin(), all.end(),
                        std::ostream_iterator< unsigned >( os, " " ),
                        std::bind( &Facet_::index, std::placeholders::_1 ) );
    }
    os << std::endl;

    os << "    non-manifold edges local IDs : ";
    for ( unsigned n = 0; n < Facet_::numNeighbors; ++n ) {
        if ( not facet->isManifold( n ) )  os << n << " ";
    }
    os << std::endl;

    os << "    manifold neighbor IDs : ";
#if 1
    for ( unsigned n = 0; n < Facet_::numNeighbors; ++n )
        if ( facet->manifoldNeighbor( n ) )
            os << facet->manifoldNeighbor( n )->index( ) << " ";
        else 
            os << "-1 ";
#else
    corlib::transform_if( facet->manNeighborsBegin( ), facet->manNeighborsEnd( ),
                          std::ostream_iterator< unsigned >( os, " " ),
                          std::bind( &Facet_::index, std::placeholders::_1 ),
                          std::bind( std::not_equal_to< Facet_ * >( ), std::placeholders::_1,
                                       static_cast< Facet_ * >( NULL ) ) );
#endif
    os << std::endl;

    //os << "    local edge tags : ";
    //std::copy( edgeTags_.begin( ), edgeTags_.end( ),
    //           std::ostream_iterator< EdgeTag >( os, " " ) );
    //os << std::endl;

    os << "    global edge tags : ";
    for ( unsigned n = 0; n < Facet_::numNeighbors; ++n )
        os << facet->getEdgeTag( n ) << " ";
    os << std::endl;

    return os;
}

//------------------------------------------------------------------------------
template< corlib::shape SHAPE, typename VERTEX >
void subdiv::surf::printConnectivityData( std::ostream & os,
                                          Facet< SHAPE, VERTEX > * facet )
{
    std::transform( facet->verticesBegin( ), facet->verticesEnd( ),
                    std::ostream_iterator< unsigned >( os, " " ),
                    std::bind( &VERTEX::index, std::placeholders::_1 ) );
}

//------------------------------------------------------------------------------
template< corlib::shape SHAPE, typename VERTEX >
std::ostream & subdiv::surf::write( std::ostream & os,
                                    FacetTree< SHAPE, VERTEX > * ftree )
{
    // a few types
    typedef VERTEX                                    Vertex;
    typedef FacetTree< SHAPE, Vertex >                FTree;
    typedef std::vector< FTree * >                    VecFTPtr;
    typedef typename VecFTPtr::iterator               VecFTPtrIter;

    // collect leafs
    VecFTPtr leafs;
    ftree->addLeafs( std::back_inserter( leafs ) );

    // loop leafs
    const VecFTPtrIter lBegin = leafs.begin();
    const VecFTPtrIter lEnd = leafs.end();
    for ( VecFTPtrIter l = lBegin; l != lEnd; ++l ) {

        os << "Facet tree " << (*l)->index( ) << std::endl;

        os << "    vertex IDs : ";
        std::transform( (*l)->verticesBegin(), (*l)->verticesEnd(),
                        std::ostream_iterator< unsigned >( os, " " ),
                        std::bind( &Vertex::index, std::placeholders::_1 ) );
        os << std::endl;

        os << "    neighbour IDs : ";
        corlib::transform_if( (*l)->neighborsBegin( ), (*l)->neighborsEnd( ),
                              std::ostream_iterator< unsigned >( os, " " ),
                              std::bind( &FTree::index, std::placeholders::_1 ),
                              std::bind( std::not_equal_to< FTree * >( ), std::placeholders::_1,
                                           static_cast< FTree * >( NULL ) ) );
        os << std::endl;

        os << "    all (non-)manifold neighbour IDs : ";
        for ( unsigned n = 0; n < FTree::numNeighbors; n++ ) {
            std::vector< FTree * > all;
            (*l)->neighborsAtEdge( n, all );
            std::transform( all.begin(), all.end(),
                            std::ostream_iterator< unsigned >( os, " " ),
                            std::bind( &FTree::index, std::placeholders::_1 ) );
        }
        os << std::endl;

        os << "    non-manifold edges local IDs : ";
        for ( unsigned n = 0; n < FTree::numNeighbors; ++n ) {
            if ( not (*l)->isManifold( n ) )  os << n << " ";
        }
        os << std::endl;

        os << "    manifold neighbor IDs : ";
#if 1
        for ( unsigned n = 0; n < FTree::numNeighbors; ++n )
            if ( (*l)->manifoldNeighbor( n ) )
                os << (*l)->manifoldNeighbor( n )->index( ) << " ";
            else
                os << "-1 ";
#else
        corlib::transform_if( (*l)->manNeighborsBegin( ), (*l)->manNeighborsEnd( ),
                              std::ostream_iterator< unsigned >( os, " " ),
                              std::bind( &FTree::index, std::placeholders::_1 ),
                              std::bind( std::not_equal_to< FTree * >( ), std::placeholders::_1,
                                           static_cast< FTree * >( NULL ) ) );
#endif
        os << std::endl;

        // os << "    local edge tags : ";
        // std::copy( edgeTags_.begin( ), edgeTags_.end( ),
        //            std::ostream_iterator< EdgeTag >( os, " " ) );
        // os << std::endl;

        os << "    global edge tags : ";
        for ( unsigned n = 0; n < FTree::numNeighbors; ++n )
            os << (*l)->getEdgeTag( n ) << " ";
        os << std::endl;

        os << "    parent ID : ";
        if ( ftree->parent() ) os << ftree->parent()->index();
        os << std::endl;

    }

    return os;
}

#endif
