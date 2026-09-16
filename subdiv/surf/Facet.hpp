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

#ifndef subdiv_surf_facet_h
#define subdiv_surf_facet_h


#include <cstdlib>
#include <iostream>
#include <set>
#include <vector>
#include <iterator>

#include <Eigen/Core>

#include <corlib/Shape.hpp>
#include <corlib/verify.hpp>
#include <corlib/misc.hpp>
#include <subdiv/surf/ShapeDim.hpp>
#include <subdiv/surf/Walker.hpp>

//------------------------------------------------------------------------------
namespace subdiv{
    namespace surf{

        template< corlib::shape SHAPE, typename V >
        class Facet;

    }
}

//------------------------------------------------------------------------------
template < corlib::shape SHAPE, typename V >
class subdiv::surf::Facet
{    
public:
    static const corlib::shape myShape           = SHAPE;
    static const unsigned      numVertices       = corlib::ShapeTraits< myShape >::numVertices; 
    static const unsigned      numNeighbors      = corlib::ShapeTraits< myShape >::numFaces;
    static const unsigned      numEdges          = corlib::ShapeTraits< myShape >::numVertices;

private:
    typedef Facet< SHAPE, V >                                Facet_;

public:
    typedef V                                                Vertex;

    typedef std::array< Vertex *, numVertices >      VecVPtrNV;
    typedef typename VecVPtrNV::const_iterator       VecVPtrNVConstIter;

    typedef std::vector< Vertex * >                  VecVPtr;
    typedef typename VecVPtr::iterator               VecVPtrIter;

    typedef std::vector< Facet_ * >                  VecFPtr;
    typedef typename VecFPtr::iterator               VecFPtrIter;

    typedef std::set< Facet_ * >                     SetFPtr;
    typedef typename SetFPtr::iterator               SetFPtrIter;

    typedef std::array< Facet_ *, numNeighbors >     VecFPtrNN;
    typedef typename VecFPtrNN::const_iterator       VecFPtrNNConstIter;

    // Edge tags
    typedef enum subdiv::surf::edgeTag               EdgeTag;
    typedef std::array< EdgeTag, numEdges >          VecETagNE;

public:
    /// Empty constructor
    Facet( ) { }

    /// Constructor
    ///
    /// \param[in]     vec    Vector with pointers to facet vertices
    Facet( VecVPtrNV vec ) : vertices_( vec )
    {
        neighbors_.fill( nullptr );
        manifoldNeighbors_.fill( nullptr );
        edgeTags_.fill( subdiv::surf::EDGE_ASKNEIGHBOR );
    }

    /// Destructor
    ~Facet( ) { }

    /// Reverse orientation of vertices, neighbors, etc
    void reverseOrientation( );
    
    /// Return i-th vertex of face
    ///
    /// \param[in]    i     Local index of vertex
    Vertex * vertex( const unsigned i ) { return vertices_[ i ]; }

    /// Get local vertex index by vertex pointer
    ///
    /// \param[in]   v   Pointer to vertex
    /// \return          Local index of vertex
    unsigned localVertexIndex( const Vertex * v ) const;

    /// Return iterator to beginning of vertices
    VecVPtrNVConstIter verticesBegin( ) const { return vertices_.begin( ); }

    /// Return iterator to end of vertices
    VecVPtrNVConstIter verticesEnd( ) const { return vertices_.end( ); }

    /// Get local edge index
    ///
    /// \param[in]   v0  Pointer to a vertex on edge
    /// \param[in]   v1  Pointer to a vertex on edge
    /// \return          Local index of edge
    int localEdgeIndex( const Vertex * v0, const Vertex * v1 );

    /// Return pair of vertices on (i)-th edge
    ///
    ///\verbatim
    ///                edge(2)
    ///               
    ///            3-----------2                          2
    ///            |           |                         / `
    ///            |           |                        /   `
    ///    edge(3) |  [this]   |  edge(1)     edge(2)  /     ` edge(1)
    ///            |           |                      /[this] `
    ///            |           |                     /         `
    ///            0-----------1                    0-----------1
    ///
    ///                edge(0)                         edge(0)
    ///\endverbatim
    ///
    /// \param[in]     i     Local edge index
    /// \return              Pair of vertices on edge
    std::pair< Vertex *, Vertex * > vertexPairOnEdge( const int i )
    {
        return std::make_pair( vertices_[ i ], vertices_[ next_[ i ] ] );
    }

    /// Set edge tag (only in top level elements)
    ///
    /// \param[in]   e      Local edge index
    /// \param[in]   eTag   The edge tag
    void setEdgeTag( const unsigned e, const EdgeTag eTag )
    {
        edgeTags_[ e ] = eTag;
        return;
    }

    /// Returns global edge tag of parent, neighbor or itself,
    /// or local edge tag of itself
    ///
    /// \param[in]   e         Local edge index
    /// \param[in]   local     If true, get local edge tag of facet, rather than globally looked-up
    /// \return                (Global or local) edge tag
    EdgeTag getEdgeTag( const unsigned e, const bool local = false );

    /// Check if manifold over edge iedge
    ///
    /// \param[in]    iedge   Edge index
    /// \return               True if manifold edge
    bool isManifold( const int iedge );

    /// Return [i]-th neighbour (see sketch for quadrilaterals) of [this] facet.
    ///
    /// <i>Please note:</i> In a non-manifold setting, the neighbours need not
    ///                     be in one "plane". Moreover, the neighbours are
    ///                     not (necessarily) unique.
    /// 
    ///\verbatim
    ///             *-----*
    ///             |     |
    ///             | [2] |
    ///             |     |
    ///             *-----*
    ///     *-----* 3-----2 *-----*            *-------* 2 *-------*
    ///     |     | |     | |     |             ` [2] / / ` ` [1] /
    ///     | [3] | | this| | [1] |              `   / /   ` `   /
    ///     |     | |     | |     |               ` / /this ` ` /
    ///     *-----* 0-----1 *-----*                * 0-------1 *
    ///             *-----*         <                 *-------*
    ///             |     |                           ` [0] /
    ///             | [0] |                            `   /
    ///             |     |                             ` /
    ///             *-----*                              *
    ///\endverbatim
    inline Facet< SHAPE, V > * & neighbor( const int i ) { return neighbors_[ i ]; }

    /// Return iterator to beginning of neighbours
    VecFPtrNNConstIter neighborsBegin( ) const { return neighbors_.begin( ); }

    /// Return iterator to end of neighbours
    VecFPtrNNConstIter neighborsEnd( ) const { return neighbors_.end( ); }

    /// Return/Set "opposite" face provided to (unique) vertex pointers pair
    /// to identify edge. The pair of vertex pointers needs not to be
    /// oriented.
    inline Facet< SHAPE, V > * & opposite( const Vertex * v1, const Vertex * v2 );

    /// Get all facets at a (non-)manifold edge
    ///
    /// \param[in]   iedge           Edge index
    /// \param[out]  eNeighbors      (Non-)manifold neigbours at edge
    void neighborsAtEdge( const int iedge, VecFPtr & eNeighbors );

    /// Get all neighbouring facet including non-manifold
    void allNeighbors( SetFPtr & allNeighbors );

    /// Give manifold neighbour
    inline Facet< SHAPE, V > * & manifoldNeighbor( const int i ) { return manifoldNeighbors_[ i ]; }

    /// Return iterator to beginning of neighbours
    VecFPtrNNConstIter manNeighborsBegin( ) const { return manifoldNeighbors_.begin( ); }

    /// Return iterator to end of neighbours
    VecFPtrNNConstIter manNeighborsEnd( ) const { return manifoldNeighbors_.end( ); }

    /// Return index
    unsigned index( ) const { return index_; } 
    /// Set index
    void setIndex( const unsigned i ) { index_ = i; return; } 

private:
    /// Permutation: Gives local index of next node depending on current node index
    static const std::array< int, corlib::ShapeTraits< SHAPE >::numVertices >  next_;
    /// Permutation: Gives local index of previous node depending on current node index
    static const std::array< int, corlib::ShapeTraits< SHAPE >::numVertices >  previous_;

    unsigned          index_;                ///< vertex index

    VecVPtrNV         vertices_;             ///< vertices
    VecFPtrNN         neighbors_;            ///< non-manifold neighboring facets
    VecFPtrNN         manifoldNeighbors_;    ///< manifold neighboring facets

    VecETagNE         edgeTags_;             ///< edge tags
};

//------------------------------------------------------------------------------
#include "Facet.ipp"
    
#endif
