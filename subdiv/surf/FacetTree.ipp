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

#include <algorithm>
#include <boost/numeric/ublas/io.hpp>

#include <subdiv/surf/Subdivision.hpp>
#include <subdiv/surf/BuildChildren.hpp>

//------------------------------------------------------------------------------
/// Local index of next vertex
template< corlib::shape SHAPE, typename V >
const std::array< int, corlib::ShapeTraits< SHAPE >::numVertices >
subdiv::surf::FacetTree< SHAPE, V >::next_ = subdiv::surf::ShapeProp< SHAPE >::next;

/// Local index of previous vertex
template< corlib::shape SHAPE, typename V >
const std::array< int, corlib::ShapeTraits< SHAPE >::numVertices >
subdiv::surf::FacetTree< SHAPE, V >::previous_ = subdiv::surf::ShapeProp< SHAPE >::previous;

//------------------------------------------------------------------------------
/// Constructor
template< corlib::shape SHAPE, typename V >
subdiv::surf::FacetTree< SHAPE, V >::FacetTree( const VecVPtrNV & vertices,
                                              const int depth,
                                              FTree * parent )
    : vertices_( vertices ),
      parent_( parent ),  // this is NULL if on top level
      children_( std::array< FTree *, numChildren > { } ),  // detection of leaf status relies on NULLs
      neighbors_( std::array< FTree *, numNeighbors > { } ),  // Mesh relies on NULL if unset
      manifoldNeighbors_( std::array< FTree *, numNeighbors > { } ),  // Mesh relies on NULL if unset
      edgeTags_( std::array< EdgeTag, numEdges > { subdiv::surf::EDGE_ASKNEIGHBOR } ),
      depth_( depth )
{ }

//------------------------------------------------------------------------------
/// Destructor
template< corlib::shape SHAPE, typename V >
subdiv::surf::FacetTree< SHAPE, V >::~FacetTree( )
{
    // delete children
    if ( not this->isLeaf( ) ) {
        std::for_each( children_.begin( ), children_.end( ),
                       corlib::deleteFunctor( ) );
    }

    // delete vertices
    std::for_each( newVertices_.begin( ), newVertices_.end( ),
                   corlib::deleteFunctor( ) );
    newVertices_.clear( );

    return;
}

//------------------------------------------------------------------------------
template< corlib::shape SHAPE, typename V >
unsigned subdiv::surf::FacetTree< SHAPE, V >::localVertexIndex( const Vertex * v ) const
{
    typename VecVPtrNV::const_iterator vIter =
        std::find( vertices_.begin(), vertices_.end(), v );
    FTL_VERIFY_DESCRIPTIVE( vIter != vertices_.end( ),
                            "Vertex is not in facet of facettree\n" );
    return std::distance( vertices_.begin(), vIter );
}

//------------------------------------------------------------------------------
template< corlib::shape SHAPE, typename V >
int subdiv::surf::FacetTree< SHAPE, V >::localEdgeIndex( const Vertex * v0,
                                                         const Vertex * v1 ) const
{
    assert( v0 );  assert( v1 );  assert( v1 != v0 );
    const int iv0 = this->localVertexIndex( v0 );
    const int iv1 = this->localVertexIndex( v1 );
    assert( iv1 == next_[ iv0 ] or iv1 == previous_[ iv0 ] );
    return ( iv1 == next_[ iv0 ] ) ? iv0 : iv1;
}

//------------------------------------------------------------------------------
template< corlib::shape SHAPE, typename V >
enum subdiv::surf::edgeTag subdiv::surf::FacetTree< SHAPE, V >::getEdgeTag(
    const unsigned e,
    const bool local
    )
{
    EdgeTag eTag = edgeTags_[ e ];
    if ( local ) {
        return eTag;
    }
    else {
//        if ( eTag == subdiv::surf::EDGE_ASKMANIFOLDNEIGHBOR ) {
//            // vertices on edge
//            const Vertex * v0 = this->vertex( e );
//            const Vertex * v1 = this->vertex( next_[ e ] );
//            assert( this->manifoldNeighbor( e ) );
//            FTree * fm = this->manifoldNeighbor( e );
//            const unsigned em = fm->localEdgeIndex( v0, v1 );
//            return fm->getEdgeTag( em );
//        }
//        else
        if ( eTag == subdiv::surf::EDGE_ASKNEIGHBOR ) {
            // vertices on edge
            const Vertex * v0 = this->vertex( e );
            const Vertex * v1 = this->vertex( next_[ e ] );
            // loop all neighbours at edge, which are stored in ring
            FTree * fn = this->neighbor( e );
            int en = fn->localEdgeIndex( v0, v1 );
            while ( fn->opposite( v0, v1 ) != this and
                    fn->getEdgeTag( en ) == subdiv::surf::EDGE_ASKNEIGHBOR ) {
                fn = fn->opposite( v0, v1 );
                en = fn->localEdgeIndex( v0, v1 );
            }
            assert( fn->getEdgeTag( en ) != subdiv::surf::EDGE_ASKNEIGHBOR );
            return fn->getEdgeTag( en );
        }
        else if ( eTag == subdiv::surf::EDGE_ASKPARENT ) {
            return parent_->getEdgeTag( e );
        }
        else {
            return eTag;
        }
    }
}

//------------------------------------------------------------------------------
template< corlib::shape SHAPE, typename V >
bool subdiv::surf::FacetTree< SHAPE, V >::isCreased( const int iedge )
{
    return ( this->getEdgeTag( iedge ) == subdiv::surf::EDGE_CREASE or
             this->getEdgeTag( iedge ) == subdiv::surf::EDGE_SMOOTHCREASE );
}

//------------------------------------------------------------------------------
/// Check if manifold over edge i
template< corlib::shape SHAPE, typename V >
bool subdiv::surf::FacetTree< SHAPE, V >::isManifold( const int iedge )
{
    bool manifold = true;
    FTree * fn = this->neighbor( iedge );
    if ( fn ) {
        const int in = fn->localVertexIndex( this->vertex( iedge ) );
        if ( fn->neighbor( previous_[ in ] ) != this )
            manifold = false;
    }
    else if ( this->parent() ) {
        FacetTree * fp = this->parent();
#if 0
        manifold = fp->isManifold( iedge );
#else
        fn = fp->neighbor( iedge );  // neighbor of parent
        if ( fn ) {
            const int in = fn->localVertexIndex( fp->vertex( iedge ) );
            if ( fn->neighbor( previous_[ in ] ) != fp )
                manifold = false;
        }
#endif
    }
    return manifold;
}

//------------------------------------------------------------------------------
/// Get neighbor opposite to edge v1-v2
template< corlib::shape SHAPE, typename V >
subdiv::surf::FacetTree< SHAPE, V > * &
subdiv::surf::FacetTree< SHAPE, V >::opposite( const Vertex * v1, const Vertex * v2 )
{
    FTL_VERIFY( v1 != NULL );
    FTL_VERIFY( v2 != NULL );

    unsigned indexOpposite = numNeighbors;
    for ( unsigned i = 0; i < numVertices; i++ ) {
        if ( v1 == vertices_[ i ] and v2 == vertices_[ next_[ i ] ] ) {
            indexOpposite = i;
            break;
        }
        else if ( v1 == vertices_[ i ] and v2 == vertices_[ previous_[ i ] ]  ) {
            indexOpposite = previous_[ i ];
            break;
        }
        else {
            ; // DO NOTHING
        }
    }
    FTL_VERIFY( indexOpposite != numNeighbors );
    return neighbors_[ indexOpposite ];
}

//------------------------------------------------------------------------------
/// Find all facet trees at a non-manifold edge
template< corlib::shape SHAPE, typename V >
void subdiv::surf::FacetTree< SHAPE, V >::neighborsAtEdge( const int iedge,
                                                         std::vector< FTree * > & eNeighbors )
{
    FTree * fn = this->neighbor( iedge );  // my neighbour
    if ( fn ) {
        eNeighbors.push_back( fn );
        // vertices on edge
        Vertex * edge1stVertex = this->vertex( iedge );
        Vertex * edge2ndVertex = this->vertex( next_[ iedge ] );
        // loop all neighbours at edge, which are stored in ring
        while ( fn->opposite( edge1stVertex, edge2ndVertex ) != this ) {
            fn = fn->opposite( edge1stVertex, edge2ndVertex );
            eNeighbors.push_back( fn );
        }
    }
    return;
}



//------------------------------------------------------------------------------
/// Check if is a leaf (has no children)
template< corlib::shape SHAPE, typename V >
bool subdiv::surf::FacetTree< SHAPE, V >::isLeaf( ) const
{
#if 0
    typename VecFTPtrNC::const_iterator child =
        std::find_if( children_.begin(), children_.end(),
                      std::not1( std::bind2nd( std::equal_to< FTree * >( ),
                                               static_cast< FTree * >( NULL ) ) ) );
    return ( child == children_.end() );  // every child is NULL
#else
    bool leaf = true;
    for ( unsigned i = 0; i < numChildren; i++ ) {
        if ( children_[ i ] ) {
            leaf = false;
            break;
        }
    }
    return leaf;
#endif
}

//------------------------------------------------------------------------------
/// Subdivide facet tree
template< corlib::shape SHAPE, typename V >
template< typename SUBDIV >
void subdiv::surf::FacetTree< SHAPE, V >::subdivide( SUBDIV * subdiv )
{
    // check if facet tree leaf
    if ( this->isLeaf() ) {
        // subdivide the leaf (no children)
        this->subdivideLeaf_( subdiv );
    }
    else {
        // subdivide each children
        for ( unsigned i = 0; i < numChildren; i++ )
            children_[ i ]->subdivide( subdiv );
    }

    return;
}

//------------------------------------------------------------------------------
/// Collect all leafs to a vector
template< corlib::shape SHAPE, typename V >
template< typename OUTITER >
void subdiv::surf::FacetTree< SHAPE, V >::addLeafs( OUTITER iter )
{
    if ( this->isLeaf() ) {
        *(iter++) = this;
    }
    else {
        for ( unsigned i = 0; i < numChildren; i++ )
            children_[ i ]->addLeafs( iter );
    }

    return;
}

//------------------------------------------------------------------------------
template< corlib::shape SHAPE, typename V >
template< typename OUTITER >
void subdiv::surf::FacetTree< SHAPE, V >::addNewVertices( OUTITER vtxs )
{
    std::copy( newVertices_.begin( ), newVertices_.end( ), vtxs );

    if ( not this->isLeaf( ) )
        for ( unsigned i = 0; i < numChildren; i++ )
            children_[ i ]->addNewVertices( vtxs );

    return;
}

//------------------------------------------------------------------------------
/// Refine the facet vertex
template< corlib::shape SHAPE, typename V >
template< typename SUBDIV >
void subdiv::surf::FacetTree< SHAPE, V >::newCoordFacetVertex_( Vertex * & vf,
                                                                SUBDIV * subdiv )
{
    vf = new Vertex( depth_ + 1 );
    typename Vertex::Point newPos = Vertex::makePoint( 0. ); // must be initialised to zero
    subdiv->refineFacet( this, depth_, newPos );
    vf->addNewLevel( newPos, depth_ + 1 );
    return;
}

//------------------------------------------------------------------------------
/// Refine the edges
template< corlib::shape SHAPE, typename V >
template< typename SUBDIV >
void subdiv::surf::FacetTree< SHAPE, V >::newCoordEdgeVertices_( VecVPtrNV vtx,
                                                                 SUBDIV * subdiv )
{
    for ( unsigned i= 0; i<numVertices; ++i ) {
        typename Vertex::Point newPos = Vertex::makePoint( 0. ); // must be initialised to zero
        subdiv->refineEdge( this, i, depth_, newPos );
        vtx[ i ]->addNewLevel( newPos, depth_ + 1 );
    }
    return;
}

//------------------------------------------------------------------------------
/// Refine the existing vertices
template< corlib::shape SHAPE, typename V >
template< typename SUBDIV >
void subdiv::surf::FacetTree< SHAPE, V >::newCoordExistingVertices_( SUBDIV * subdiv )
{
    for ( unsigned i=0; i<numVertices; ++i ) {

        const int maxLev = vertices_[ i ]->maxActiveLevel();

        if ( ( depth_ + 1 ) != maxLev ) {

            typename Vertex::Point newPos = Vertex::makePoint( 0. ); // must be initialised to zero
            subdiv->refineVertex( this, i, depth_, newPos );

            vertices_[ i ]->addNewLevel( newPos, depth_ + 1 );
        }
    }

    return;
}

//------------------------------------------------------------------------------
/// Make or find edge vertex
template< corlib::shape SHAPE, typename V >
V * subdiv::surf::FacetTree< SHAPE, V >::makeOrFindEdgeMidVertex_( const int iedge,
                                                                   VecVPtr & vNew )
{
    Vertex * midVertex = NULL;
    bool found = false;
    FTree * fn = neighbors_[ iedge ];

    if ( this->isManifold( iedge ) ) {
        // manifold edge

        if ( fn and not fn->isLeaf( ) ) {
            // the neighbour is already refined and created the edge vertex
            found = true;
            // link vertex
            midVertex = this->findNeighborEdgeMidVertex_( iedge );
        }
    }
    else {
        // non-manifold edge

        // get all facet trees connected at the edge
        std::vector< FTree * > ftVec;
        this->neighborsAtEdge( iedge, ftVec );

        // check if any of the connected facet trees has the vertex
        for ( unsigned i = 0; i < ftVec.size(); i++ ) {
            if ( not ftVec[ i ]->isLeaf() ) {
                found = true;

                // find neighbour edge id
                const int iVtx0 = ftVec[ i ]->localVertexIndex( vertices_[ iedge ] );
                const int iVtx1 = ftVec[ i ]->localVertexIndex( vertices_[ next_[ iedge ] ] );
                const int nEdge = ( iVtx1 == next_[ iVtx0 ] ) ? iVtx0 : previous_[ iVtx0 ];

                midVertex = ftVec[ i ]->child( nEdge )->vertex( next_[ nEdge ] );
                break;
            }
        }
    }

    // check if vertex has not been found
    if ( not found ) {
        midVertex = new Vertex( depth_ + 1 );
        assert( midVertex != NULL );
        newVertices_.push_back( midVertex );

        // properly mark vertex on creased edge
        if ( this->getEdgeTag( iedge ) == subdiv::surf::EDGE_CREASE )
            midVertex->setTag( subdiv::surf::VERTEX_CREASE );
        else if ( this->getEdgeTag( iedge ) == subdiv::surf::EDGE_SMOOTHCREASE )
            midVertex->setTag( subdiv::surf::VERTEX_SMOOTHCREASE );

    }
    return midVertex;
}

//------------------------------------------------------------------------------
/// Find the edge vertex from neighbour
template< corlib::shape SHAPE, typename V >
V * subdiv::surf::FacetTree< SHAPE, V >::findNeighborEdgeMidVertex_( const int iedge )
{
    Vertex * v0 = vertices_[ iedge ];
    FTree * ftN = neighbors_[ iedge ];
    const int iedgeN = ftN->localVertexIndex( v0 );
    Vertex * midVertex = ftN->child( iedgeN )->vertex( previous_[ iedgeN ] );
    return midVertex;
}

//------------------------------------------------------------------------------
/// Create new children during facet tree leaf subdivision
template< corlib::shape SHAPE, typename V >
template< typename SUBDIV >
void subdiv::surf::FacetTree< SHAPE, V >::createNewElements_( SUBDIV * subdiv )
{
    // a) compute and create new facet vertex coordinate
    Vertex * vFace = NULL;
    this->newCoordFacetVertex_( vFace, subdiv );
    newVertices_.push_back( vFace );

    // b) create new or find existing mid-point vertices on edges
    VecVPtrNV vEdge;
    for ( unsigned i= 0 ; i < numVertices; i++ )
        vEdge[ i ] = this->makeOrFindEdgeMidVertex_( i, newVertices_ );
    // compute new edge vertex coordinates
    this->newCoordEdgeVertices_( vEdge, subdiv );

    // c) recompute existing coordinates
    this->newCoordExistingVertices_( subdiv );

    // build the children and their relationships
    subdiv::surf::BuildChildren< myShape, FTree >()( this, vFace, vEdge, newVertices_ );
    return;
}

//------------------------------------------------------------------------------
/// Subdivide facet tree leaf after checking neighbour refinement levels
template< corlib::shape SHAPE, typename V >
template< typename SUBDIV >
void subdiv::surf::FacetTree< SHAPE, V >::subdivideLeaf_( SUBDIV * subdiv )
{
    // check if level difference is not to high
    for ( unsigned n = 0; n < numVertices; ++n ) {

        int edgeId = previous_[ n ];
        bool edgeFound = false;
        const Vertex * v = this->vertex( n );
        FacetTree * fn = neighbors_[ previous_[ n ] ]; // neighbour
        FacetTree * fp = this; // parent

        // circulate anti-clockwise
        while ( ( not edgeFound ) and ( fn != this ) ) {
            if ( not fp->isManifold( edgeId ) ) {
                // non-manifold edge found
                if ( fp->parent() )
                    fp->subdivideNeighbors_( edgeId, subdiv );
                edgeFound = true;
            }
            else if ( ( not fn ) and fp->parent() ) {
                // no neighbour but child, i.e. neighbour insufficiently refined
                fp->subdivideNeighbors_( edgeId, subdiv );  // subdivide parent neighbour at edge
                fn = fp->neighbor( edgeId );
                if ( not fn ) {
                    // border found
                    edgeFound = true;
                }
                else {
                    fp = fn;
                    edgeId = previous_[ fp->localVertexIndex( v ) ];
                    fn = fp->neighbor( edgeId );
                }
            }
            else if ( fn ) {
                // Neighbour in same level
                fp = fn;
                edgeId = previous_[ fp->localVertexIndex( v ) ];
                fn = fp->neighbor( edgeId );
            }
            else {
                // border
                edgeFound = true;
            }
        }

        // rotate clockwise if edge is found
        if ( edgeFound ) {
            int edgeId = n;
            edgeFound = false;
            fn = neighbors_[ n ];
            fp = this;
            while ( not edgeFound ) {
                if ( not fp->isManifold( edgeId ) ) {
                    // non-manifold edge found
                    if ( fp->parent() )
                        fp->subdivideNeighbors_( edgeId, subdiv );
                    edgeFound = true;
                }
                else if ( ( not fn ) and fp->parent() ) {
                    // no neighbour but child
                    fp->subdivideNeighbors_( edgeId, subdiv );
                    fn = fp->neighbor( edgeId );
                    if ( not fn ) {
                        //border found
                        edgeFound = true;
                    }
                    else {
                        fp = fn;
                        edgeId = fp->localVertexIndex( v );
                        fn = fp->neighbor( edgeId );
                    }
                }
                else if (fn) {
                    // has neighbour in same level, go ahead
                    fp = fn;
                    edgeId = fp->localVertexIndex( v );
                    fn = fp->neighbor( edgeId );
                }
                else {
                    // border
                    edgeFound = true;
                }
            }
        }
    }  // loop over edges

    this->createNewElements_( subdiv );

    return;
}

//------------------------------------------------------------------------------
/// Subdivide all facet trees connected to edge
template< corlib::shape SHAPE, typename V >
template< typename SUBDIV >
void subdiv::surf::FacetTree< SHAPE, V >::subdivideNeighbors_( const int n,
                                                               SUBDIV * subdiv )
{
    assert( parent_ );

    // check if manifold
    if ( parent_->isManifold( n ) ) {
        // check if parent has neighbour
        FTree * fn = parent_->neighbor( n );
        if ( fn and fn->isLeaf() ) {
            fn->subdivideLeaf_( subdiv );
        }
    }
    else {
        // get vector of all facet trees connected at the edge
        std::vector< FTree * > ftVec;
        parent_->neighborsAtEdge( n, ftVec );

        // subdivide all facet trees connected to the edge
        for ( unsigned i = 0; i < ftVec.size(); i++ ) {
            if ( ftVec[ i ]->isLeaf() ) {
                ftVec[ i ]->subdivide( subdiv );
            }
        }
    }

    return;
}
