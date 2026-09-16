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

#ifndef subdiv_surf_facettree_h
#define subdiv_surf_facettree_h

#include <cstdlib>
#include <iterator>
#include <functional>
#include <algorithm>

#include <corlib/verify.hpp>
#include <corlib/Shape.hpp>
#include <corlib/misc.hpp>

#include <subdiv/surf/ShapeDim.hpp>
#include <subdiv/surf/MeshTags.hpp>

//------------------------------------------------------------------------------
namespace subdiv{
    namespace surf{

        template< corlib::shape SHAPE, typename V >
        class FacetTree;

    }
}

//------------------------------------------------------------------------------
/// ???
///
/// \tparam  SHAPE     Simplex shape type
/// \tparam  V         Vertex type
template< corlib::shape SHAPE, typename V >
class subdiv::surf::FacetTree
{
public:
    static const corlib::shape myShape      = SHAPE;
    static const unsigned      numVertices  = corlib::ShapeTraits< myShape >::numVertices;
    static const unsigned      numChildren  = 4;
    static const unsigned      numNeighbors = corlib::ShapeTraits< myShape >::numFaces;
    static const unsigned      numEdges     = corlib::ShapeTraits< myShape >::numVertices;

    typedef V                                       Vertex;
    typedef FacetTree< myShape, V >                 FTree;

    /// Vector of pointers to vertices of facet
    typedef std::array< Vertex *, numVertices >     VecVPtrNV;
    typedef typename VecVPtrNV::const_iterator      VecVPtrNVConstIter;
    /// Set of pointers of vertices
    typedef std::vector< Vertex * >                 VecVPtr;
    typedef typename VecVPtr::iterator              VecVPtrIter;
    /// Vector of pointer to children vertices
    typedef std::array< FTree *, numChildren >      VecFTPtrNC;
    typedef typename VecFTPtrNC::iterator           VecFTPtrNCItr;
    typedef std::array< FTree *, numNeighbors >     VecFTPtrNN;
    typedef typename VecFTPtrNN::const_iterator     VecFTPtrNNConstIter;
    // Edge tags
    typedef enum subdiv::surf::edgeTag              EdgeTag;
    typedef std::array< EdgeTag, numEdges >         VecETagNE;

private:

public:
    /// Constructor
    ///
    /// \param[in]    vertices  Vector of pointers to vertices
    /// \param[in]    depth     Refinement level
    /// \param[in]    parent    Parent facet, NULL in case of top level
    FacetTree( const VecVPtrNV & vertices, const int depth, FTree * parent );

    /// Destructor
    ~FacetTree( );

    //! Get vertex pointer by local vertex index
    ///
    ///\verbatim
    ///                                 2
    ///       3--------2               / `
    ///       |        |              /   `
    ///       | [this] |             /     `
    ///       |        |            / [this]`
    ///       0--------1           0---------1
    ///
    ///\endverbatim
    ///
    /// \param[in]   i   Local vertex index
    /// \return          Vertex there
    Vertex * & vertex( const int i ) { return vertices_[ i ]; }

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
    int localEdgeIndex( const Vertex * v0, const Vertex * v1 ) const;

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

    /// Determine if edge is creased
    bool isCreased( const int iedge );

    /// Check if manifold _or_ boundary over edge iedge
    ///
    /// \param[in]    iedge   Edge index
    /// \return               True if manifold edge, _or_ boundary edge
    bool isManifold( const int iedge );

    /// Get [i]-th neighbor of [this] facet
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
    ///             *-----*                          *-------*
    ///             |     |                           ` [0] /
    ///             | [0] |                            `   /
    ///             |     |                             ` /
    ///             *-----*                              *
    ///\endverbatim
    ///
    /// \param[in]     i     i-th neighbour to get
    /// \return              Pointer to i-th neighbour, may be NULL
    FTree * & neighbor( const int i ) { return neighbors_[ i ]; }

    /// Return iterator to beginning of neighbours
    VecFTPtrNNConstIter neighborsBegin( ) const { return neighbors_.begin( ); }

    /// Return iterator to end of neighbours
    VecFTPtrNNConstIter neighborsEnd( ) const { return neighbors_.end( ); }

    /// Get neighbor opposite to edge v1-v2
    ///
    /// \param[in]  v1    A vertex on edge
    /// \param[in]  v2    Another vertex on edge
    /// \return           Opposite neighbour at facet edge
    FacetTree< SHAPE, V > * & opposite( const Vertex * v1, const Vertex * v2 );

    /// Get all facet trees at a (non-)manifold edge
    ///
    /// \param[in]   iedge           Edge index
    /// \param[out]  allNeighbors    (Non-)manifold neigbours at edge
    void neighborsAtEdge( const int iedge, std::vector< FTree * > & eNeighbors );

    /// Get [i]-th manifold neighbor of [this] facet
    ///
    /// \param[in]     i     i-th neighbour to get
    /// \return              Pointer to i-th neighbour, may be NULL
    FTree * & manifoldNeighbor( const int i ) { return manifoldNeighbors_[ i ]; }

    /// Return iterator to beginning of neighbours
    VecFTPtrNNConstIter manNeighborsBegin( ) const { return manifoldNeighbors_.begin( ); }

    /// Return iterator to end of neighbours
    VecFTPtrNNConstIter manNeighborsEnd( ) const { return manifoldNeighbors_.end( ); }

    /// Get parent
    ///
    /// \return           Pointer to parent element (or NULL)
    FTree * & parent( ) { return parent_; }

    /// Get [i]-th child of [this] facet
    ///
    ///\verbatim
    ///                                              2
    ///      3-----------------2                    / `
    ///      | 3-----2 3-----2 |                   / 2 `
    ///      | |     | |     | |                  / / ` `
    ///      | | [3] | | [2] | |                 / /[2]` `
    ///      | |     | |     | |                / /     ` `
    ///      | 0-----1 0-----1 |               / 0-------1 `
    ///      |     [this]      |              /   [this]    `
    ///      | 3-----2 3-----2 |             / 2 2-------1 2 `
    ///      | |     | |     | |            / / ` ` [3] / / ` `
    ///      | | [0] | | [1] | |           / /[0]` `   / /[1]` `
    ///      | |     | |     | |          / /     ` ` / /     ` `
    ///      | 0-----1 0-----1 |         / 0-------1 0 0-------1 `
    ///      0-----------------1        0-------------------------1
    ///
    ///\endverbatim
    FTree * & child( const int i ) { return children_[ i ]; }

    //! Check if it is a leaf (has no children)
    bool isLeaf( ) const;

    //! Get depth of refinement
    int depth( ) const { return depth_; }

    //! TwoDide this facet tree
    template< typename SUBDIV >
    void subdivide( SUBDIV * subdiv );

    //@}

    /// @name Output
    //@{
    /// Return Index
    unsigned index( ) const { return index_; }

    /// Set index
    void setIndex( const unsigned i ) { index_ = i; }

    //! Collect all leafs to a vector
    template< typename OUTITER >
    void addLeafs( OUTITER iter );

    /// Collect all new vertices
    template< typename OUTITER >
    void addNewVertices( OUTITER vtxs );
    //@}

private:
    /// Refine the facet vertex
    template< typename SUBDIV >
    void newCoordFacetVertex_( Vertex * & vf, SUBDIV * subdiv );

    //! Refine the edges
    template< typename SUBDIV >
    void newCoordEdgeVertices_( VecVPtrNV vtx, SUBDIV * subdiv );

    //! Refine the existing vertices
    template< typename SUBDIV >
    void newCoordExistingVertices_( SUBDIV * subdiv );

    /// Find edge vertex from neighbor facet attached to same edge
    ///
    /// \param[in]   iedge     Edge index
    Vertex * findNeighborEdgeMidVertex_( const int iedge );

    //! Create new children during facet tree leaf subdivision
    template< typename SUBDIV >
    void createNewElements_( SUBDIV * subdiv );

    /// Subdivide facet tree leaf after triggering refinement of neighbors
    /// to allow maximal one level difference in refinement at neighbours
    ///
    /// \tparam       SUBDIV    Subdivision scheme
    ///
    /// \param[in]    subdiv    Subdivisor
    template< typename SUBDIV >
    void subdivideLeaf_( SUBDIV * subdiv );

    /// Subdivide all neighbour facet trees connected to edge
    ///
    /// \tparam       SUBDIV    Subdivision scheme
    ///
    /// \param[in]    n         Local index of edge
    /// \param[in]    subdiv    Subdivisor
    template< typename SUBDIV >
    void subdivideNeighbors_( const int n, SUBDIV * subdiv );

protected:
    /// Return iterator to beginning of children
    VecFTPtrNCItr cBegin( ) { return children_.begin( ); }

    /// Return iterator to end of children
    VecFTPtrNCItr cEnd( ) { return children_.end( ); }

    /// Clear new vertices
    void clearNewVertices( ) { newVertices_.clear( ); }

    /// Clear children
    void clearChildren( ) { children_.clear( ); }

    FTree* returnSelf() { return this; }

    /// Make or find edge vertex
    ///
    /// \param[in]       iedge     Edge index
    /// \param[in,out]   vNew      Newly allocated vertices
    Vertex * makeOrFindEdgeMidVertex_( const int iedge,
                                       VecVPtr & vNew );

private:
    /// Permutation: Gives local index of next node depending on current node index
    static const std::array< int, corlib::ShapeTraits< SHAPE >::numVertices >       next_;
    /// Permutation: Gives local index of previous node depending on current node index
    static const std::array< int, corlib::ShapeTraits< SHAPE >::numVertices >       previous_;

    VecVPtrNV                  vertices_;            ///< vertices

    FTree *                    parent_;              ///< parent pointer
    VecFTPtrNC                 children_;            ///< child pointers
    VecFTPtrNN                 neighbors_;           ///< neighbors
    VecFTPtrNN                 manifoldNeighbors_;   ///< manifold neighbors
    VecETagNE                  edgeTags_;            ///< edge tags
    unsigned                   index_;               ///< index
    int                        depth_;               ///< depth of refinement

protected:
    VecVPtr                    newVertices_;         ///< newly allocated vertices used by children

};


#include "FacetTree.ipp"

#endif
