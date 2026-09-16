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

//------------------------------------------------------------------------------
//! Local index of next vertex
template< corlib::shape SHAPE, typename V >
const std::array< int, corlib::ShapeTraits< SHAPE >::numVertices > 
subdiv::surf::Facet< SHAPE, V >::next_ = subdiv::surf::ShapeProp< SHAPE >::next;

//------------------------------------------------------------------------------
//! Local index of previous vertex
template< corlib::shape SHAPE, typename V >
const std::array< int, corlib::ShapeTraits< SHAPE >::numVertices >
subdiv::surf::Facet< SHAPE, V >::previous_ = subdiv::surf::ShapeProp< SHAPE >::previous;

//------------------------------------------------------------------------------
template< corlib::shape SHAPE, typename V >
void subdiv::surf::Facet< SHAPE, V >::reverseOrientation( )
{
    // reverse vertices
    std::reverse( vertices_.begin(), vertices_.end() );
    // reverse neighbours
    std::reverse( neighbors_.begin(), neighbors_.end() );
    std::rotate( neighbors_.begin(), (neighbors_.begin())+1, neighbors_.end() );
    // reverse manifold neighbours
    std::reverse( manifoldNeighbors_.begin(), manifoldNeighbors_.end() );
    std::rotate( manifoldNeighbors_.begin(), (manifoldNeighbors_.begin())+1,
                 manifoldNeighbors_.end() );
    return;
}

//------------------------------------------------------------------------------
template< corlib::shape SHAPE, typename V >
unsigned subdiv::surf::Facet< SHAPE, V >::localVertexIndex( const Vertex * v ) const
{
    typename VecVPtrNV::const_iterator vIter =
        std::find( vertices_.begin(), vertices_.end(), v );
    FTL_VERIFY_DESCRIPTIVE( vIter != vertices_.end( ),
                            "Vertex is not in facet\n" );
    return std::distance( vertices_.begin(), vIter );
}

//------------------------------------------------------------------------------
template< corlib::shape SHAPE, typename V >
int subdiv::surf::Facet< SHAPE, V >::localEdgeIndex( const Vertex * v0,
                                                     const Vertex * v1 )
{
    assert( v0 );  assert( v1 );  assert( v1 != v0 );
    const int iv0 = this->localVertexIndex( v0 );
    const int iv1 = this->localVertexIndex( v1 );
    assert( iv1 == next_[ iv0 ] or iv1 == previous_[ iv0 ] );
    return ( iv1 == next_[ iv0 ] ) ? iv0 : iv1;
}

//------------------------------------------------------------------------------
template< corlib::shape SHAPE, typename V >
enum subdiv::surf::edgeTag subdiv::surf::Facet< SHAPE, V >::getEdgeTag(
    const unsigned e,
    const bool local
    )
{
    EdgeTag eTag = edgeTags_[ e ];
    if ( local ) {
        return eTag;
    }
    else {
        if ( eTag == subdiv::surf::EDGE_ASKNEIGHBOR ) {
            Vertex * v0 = this->vertex( e );
            Vertex * v1 = this->vertex( next_[ e ] );
            // loop all neighbours at edge, which are stored in ring
            Facet< SHAPE, V > * fn = this->neighbor( e );
            int en = fn->localEdgeIndex( v0, v1 );
            while ( fn->opposite( v0, v1 ) != this and
                    fn->getEdgeTag( en ) == subdiv::surf::EDGE_ASKNEIGHBOR ) {
                fn = fn->opposite( v0, v1 );
                en = fn->localEdgeIndex( v0, v1 );
            }
            assert( fn->getEdgeTag( en ) != subdiv::surf::EDGE_ASKNEIGHBOR );
            return fn->getEdgeTag( en );
        }
        else {
            return eTag;
        }
    }
}

//------------------------------------------------------------------------------
/// Check if manifold over edge i
template< corlib::shape SHAPE, typename V >
bool subdiv::surf::Facet< SHAPE, V >::isManifold( const int iedge )
{
    bool manifold = true;
    Facet< SHAPE, V > * fn = this->neighbor( iedge );
    if ( fn ) {
#if 1
        // works for non-aligned manifolds too
        if ( std::find( fn->neighborsBegin(), fn->neighborsEnd(), this ) ==
             fn->neighborsEnd() ) {
            manifold = false;
        }
#else
        // works only on aligned manifolds
        const int in = fn->localVertexIndex( this->vertex( iedge ) );
        if ( fn->neighbor( previous_[ in ] ) != this ) {
            manifold = false;
        }
#endif
    }
    return manifold;
}

//------------------------------------------------------------------------------
//! Give opposite facet
template< corlib::shape SHAPE, typename V >
inline subdiv::surf::Facet< SHAPE, V > * &
subdiv::surf::Facet< SHAPE, V >::opposite( const Vertex * v1, const Vertex * v2 )
{    
    FTL_VERIFY( v1 != NULL );
    FTL_VERIFY( v2 != NULL );
    
    unsigned indexOpposite = numNeighbors;
    for ( unsigned i = 0; i < numVertices; i++ ) {
        // counter-clockwise orientation
        if ( v1 == vertices_[ i ] and v2 == vertices_[ next_[ i ] ] ) {
            indexOpposite = i;
            break;
        }
        // clockwise orientation
        else if ( v1 == vertices_[ i ] and v2 == vertices_[ previous_[ i ] ] ) {
            indexOpposite = previous_[ i ];
            break;
        }
    }
    FTL_VERIFY( indexOpposite != numNeighbors );
    return neighbors_[ indexOpposite ];
}

//------------------------------------------------------------------------------
/// Find all facet trees at a non-manifold edge
template< corlib::shape SHAPE, typename V >
void subdiv::surf::Facet< SHAPE, V >::neighborsAtEdge( const int iedge,
                                                       VecFPtr & eNeighbors )
{
    Facet< SHAPE, V > * fn = this->neighbor( iedge );  // my neighbour
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
/// Find all facet trees at a non-manifold edge
template< corlib::shape SHAPE, typename V >
void subdiv::surf::Facet< SHAPE, V >::allNeighbors( SetFPtr & allNeighbors )
{
    for ( unsigned n = 0; n < numNeighbors; ++n ) {
        VecFPtr eNeighbors;
        this->neighborsAtEdge( n, eNeighbors );
        allNeighbors.insert( eNeighbors.begin( ), eNeighbors.end( ) );
    }
    return;
}

