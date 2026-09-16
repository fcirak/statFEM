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
/// Local index of next vertex
template< typename FT >
const std::array< int, corlib::ShapeTraits< corlib::TRIANGLE >::numVertices >
subdiv::surf::BuildChildren< corlib::TRIANGLE, FT >::next_ =
    subdiv::surf::ShapeProp< corlib::TRIANGLE >::next;

/// Local index of previous vertex
template< typename FT >
const std::array< int,  corlib::ShapeTraits< corlib::TRIANGLE >::numVertices >
subdiv::surf::BuildChildren< corlib::TRIANGLE, FT >::previous_ =
    subdiv::surf::ShapeProp< corlib::TRIANGLE >::previous;

//------------------------------------------------------------------------------
template< typename FT >
void subdiv::surf::BuildChildren< corlib::TRIANGLE, FT >::operator()(
    FTree * ft,  
    Vertex * fv, 
    const VecVPtrNV & ev, 
    VecVPtr & nv
    )
{
    static_assert( myShape == FT::myShape );

    // no need for facet vertex
    delete fv;
    nv.erase( std::find( nv.begin(), nv.end(), fv ) );

    // make children 0, 1, 2
    for ( unsigned i = 0; i < FTree::numChildren - 1; ++i ) {
        VecVPtrNV vtxs;
        vtxs[ i ] = ft->vertex( i );
        vtxs[ next_[ i ] ] = ev[ i ];
        vtxs[ previous_[ i ] ] = ev[ previous_[ i ] ];

        // create child
        ft->child( i ) = new FTree( vtxs, ft->depth()+1, ft );
        assert( ft->child( i ) != NULL );

        // relate edge tags
        ft->child( i )->setEdgeTag( i, subdiv::surf::EDGE_ASKPARENT );
        ft->child( i )->setEdgeTag( previous_[ i ], subdiv::surf::EDGE_ASKPARENT );
    }

    // center child 3
    ft->child( 3 ) = new FTree( ev, ft->depth()+1, ft );
    assert( ft->child( 3 ) != NULL );

    // neighbor relations with the center element
    for ( unsigned i = 0; i < FTree::numNeighbors; ++i ) {
        // non-manifolf neighbours
        ft->child( i )->neighbor( next_[ i ] ) = ft->child( 3 );
        ft->child( 3 )->neighbor( i ) = ft->child( next_[ i ] );
        // manifold neighbours
        ft->child( i )->manifoldNeighbor( next_[ i ] ) = ft->child( 3 );
        ft->child( 3 )->manifoldNeighbor( i ) = ft->child( next_[ i ] );
        // edge tag
        ft->child( i )->setEdgeTag( next_[ i ], subdiv::surf::EDGE_NOTAG );
        ft->child( 3 )->setEdgeTag( i, subdiv::surf::EDGE_ASKNEIGHBOR );
    }

    // non-manifold neighbours among other elements
    for ( unsigned i = 0; i < FTree::numChildren-1; ++i ) {

        // first edge "i"
        FTree * fn = ft->neighbor( i );
        FTree * mn = ft->manifoldNeighbor( i );
        if ( ft->isManifold( i ) ) {
            if ( ( fn ) and ( not fn->isLeaf() ) ) {
                // neighbor exists and has children
                const int j = fn->localVertexIndex( ft->vertex( i ) );
                ft->child( i )->neighbor( i ) = fn->child( j );
                fn->child( j )->neighbor( previous_[ j ] ) = ft->child( i );

                assert( fn == mn );
                ft->child( i )->manifoldNeighbor( i ) = mn->child( j );
                mn->child( j )->manifoldNeighbor( previous_[ j ] ) = ft->child( i );
            }
        }
        else {
            Vertex * v0 = ft->vertex( i );
            Vertex * v1 = ft->vertex( next_[ i ] );
            FTree * fp = fn;
            std::vector< FTree * > mps;
            while ( fp->neighbor( fp->localEdgeIndex( v0, v1 ) ) != ft ) {
                if ( fp->manifoldNeighbor( fp->localEdgeIndex( v0, v1 ) ) == ft )
                    mps.push_back( fp );
                fp = fp->neighbor( fp->localEdgeIndex( v0, v1 ) );
            }
            if ( fp->manifoldNeighbor( fp->localEdgeIndex( v0, v1 ) ) == ft )
                mps.push_back( fp );
            if ( not fn->isLeaf() ) {
                // forward neighbor has children
                const int j = fn->localVertexIndex( ft->vertex( i ) );
                if ( fn->localEdgeIndex( v0, v1 ) == j )  // negatively oriented
                    ft->child( i )->neighbor( i ) = fn->child( j );
                else                                      // positively oriented
                    ft->child( i )->neighbor( i ) = fn->child( j );
            }
            if ( not fp->isLeaf() ) {
                // backward neighbor has children
                const int j = fp->localVertexIndex( ft->vertex( i ) );
                if ( fp->localEdgeIndex( v0, v1 ) == j )  // negatively oriented
                    fp->child( j )->neighbor( j ) = ft->child( i );
                else                                      // positively oriented
                    fp->child( j )->neighbor( previous_[ j ] ) = ft->child( i );
            }
            if ( mn and not mn->isLeaf() ) {
                // forward neighbor has children
                const int j = mn->localVertexIndex( ft->vertex( i ) );
                if ( mn->localEdgeIndex( v0, v1 ) == j )  // negatively oriented
                    ft->child( i )->manifoldNeighbor( i ) = mn->child( j );
                else                                      // positively oriented
                    ft->child( i )->manifoldNeighbor( i ) = mn->child( j );
            }
            for ( unsigned k = 0; k < mps.size( ); ++k ) {
                FTree * mp = mps[ k ];
                if ( mp and not mp->isLeaf() ) {
                    // backward neighbor has children
                    const int j = mp->localVertexIndex( ft->vertex( i ) );
                    if ( mp->localEdgeIndex( v0, v1 ) == j )  // negatively oriented
                        mp->child( j )->manifoldNeighbor( j ) = ft->child( i );
                    else                                      // positively oriented
                        mp->child( j )->manifoldNeighbor( previous_[ j ] ) = ft->child( i );
                }
            }
        }

        // second edge "--i"
        fn = ft->neighbor( previous_[ i ] );
        mn = ft->manifoldNeighbor( previous_[ i ] );
        if ( ft->isManifold( previous_[ i ] ) ) {
            if ( ( fn ) and ( not fn->isLeaf() ) ) {
                // neighbor exists and has children
                const int j = fn->localVertexIndex( ft->vertex( i ) );
                ft->child( i )->neighbor( previous_[ i ] ) = fn->child( j );
                fn->child( j )->neighbor( j ) = ft->child( i );

                assert( fn == mn );  // implies mn != NULL
                ft->child( i )->manifoldNeighbor( previous_[ i ] ) = mn->child( j );
                mn->child( j )->manifoldNeighbor( j ) = ft->child( i );
            }
        }
        else {
            Vertex * v0 = ft->vertex( previous_[ i ] );
            Vertex * v1 = ft->vertex( i );
            FTree * fp = fn;
            std::vector< FTree * > mps;
            while ( fp->neighbor( fp->localEdgeIndex( v0, v1 ) ) != ft ) {
                if ( fp->manifoldNeighbor( fp->localEdgeIndex( v0, v1 ) ) == ft )
                    mps.push_back( fp );
                fp = fp->neighbor( fp->localEdgeIndex( v0, v1 ) );
            }
            if ( fp->manifoldNeighbor( fp->localEdgeIndex( v0, v1 ) ) == ft )
                mps.push_back( fp );
            if ( not fn->isLeaf() ) {
                // forward neighbor has children
                const int j = fn->localVertexIndex( ft->vertex( i ) );
                if ( fn->localEdgeIndex( v0, v1 ) == j )  // fn negatively oriented
                    ft->child( i )->neighbor( previous_[ i ] ) = fn->child( j );
                else                                      // fn positively oriented
                    ft->child( i )->neighbor( previous_[ i ] ) = fn->child( j );
            }
            if ( not fp->isLeaf() ) {
                // backward neighbor has children
                const int j = fp->localVertexIndex( ft->vertex( i ) );
                if ( fp->localEdgeIndex( v0, v1 ) != j )  // fp negatively oriented
                    fp->child( j )->neighbor( previous_[ j ] ) = ft->child( i );
                else                                      // fp positively oriented
                    fp->child( j )->neighbor( j ) = ft->child( i );
            }
            if ( mn and not mn->isLeaf() ) {
                // forward neighbor has children
                const int j = mn->localVertexIndex( ft->vertex( i ) );
                if ( mn->localEdgeIndex( v0, v1 ) == j )  // mn negatively oriented
                    ft->child( i )->manifoldNeighbor( previous_[ i ] ) = mn->child( j );
                else                                      // mn positively oriented
                    ft->child( i )->manifoldNeighbor( previous_[ i ] ) = mn->child( j );
            }
            for ( unsigned k = 0; k < mps.size( ); ++k ) {
                FTree * mp = mps[ k ];
                if ( mp and not mp->isLeaf() ) {
                    // backward neighbor has children
                    const int j = mp->localVertexIndex( ft->vertex( i ) );
                    if ( mp->localEdgeIndex( v0, v1 ) != j )  // mp negatively oriented
                        mp->child( j )->manifoldNeighbor( previous_[ j ] ) = ft->child( i );
                    else                                      // mp positively oriented
                        mp->child( j )->manifoldNeighbor( j ) = ft->child( i );
                }
            }
        }
    }

    return;
}

//------------------------------------------------------------------------------
/// Local index of next vertex
template< typename FT >
const std::array< int, corlib::ShapeTraits< corlib::QUADRILATERAL >::numVertices >
subdiv::surf::BuildChildren< corlib::QUADRILATERAL, FT >::next_ =
    subdiv::surf::ShapeProp< corlib::QUADRILATERAL >::next;

/// Local index of previous vertex
template< typename FT >
const std::array< int, corlib::ShapeTraits< corlib::QUADRILATERAL >::numVertices >
subdiv::surf::BuildChildren< corlib::QUADRILATERAL, FT >::previous_ =
    subdiv::surf::ShapeProp< corlib::QUADRILATERAL >::previous;

//------------------------------------------------------------------------------
template< typename FT >
void subdiv::surf::BuildChildren< corlib::QUADRILATERAL, FT >::operator()(
    FTree * ft, 
    Vertex * fv, 
    const VecVPtrNV & ev, 
    VecVPtr & nv
    )
{
    static_assert( myShape == FT::myShape );

    // make new children
    for ( unsigned i = 0; i < FTree::numChildren; ++i ) {
        // collect the vertices for new children
        VecVPtrNV vtxs;
        vtxs[ i ] = ft->vertex( i );
        vtxs[ next_[ i ] ] = ev[ i ];
        vtxs[ next_[ next_[ i ] ] ] = fv;
        vtxs[ previous_[ i ] ] = ev[ previous_[ i ] ];

        // create child
        ft->child( i ) = new FTree( vtxs, ft->depth()+1, ft );
        assert( ft->child( i ) != NULL );
    }

    // neighbor relations within children
    for ( unsigned i = 0; i < FTree::numNeighbors; ++i ) {
        ft->child( i )->neighbor( next_[ i ] )          = ft->child( next_[ i ] );
        ft->child( i )->neighbor( next_[ next_[ i ] ] ) = ft->child( previous_[ i ] );

        ft->child( i )->manifoldNeighbor( next_[ i ] )          = ft->child( next_[ i ] );
        ft->child( i )->manifoldNeighbor( next_[ next_[ i ] ] ) = ft->child( previous_[ i ] );
    }

    // relate edges
    for ( unsigned i = 0; i < FTree::numChildren; ++i ) {
        ft->child( i )->setEdgeTag( previous_[ i ], subdiv::surf::EDGE_ASKPARENT );
        ft->child( i )->setEdgeTag( i, subdiv::surf::EDGE_ASKPARENT );
        ft->child( i )->setEdgeTag( next_[ i ], subdiv::surf::EDGE_NOTAG );
        ft->child( i )->setEdgeTag( next_[ next_[ i ] ], subdiv::surf::EDGE_ASKNEIGHBOR );
    }

    // (non-manifold) neighbor relations to other facets : loop over edges
    for ( unsigned i = 0; i < FTree::numVertices; ++i ) {
        // get vector of all facet trees connected at the edge
        std::vector< FTree * > ftVec;
        ft->neighborsAtEdge( i, ftVec );

        // loop neighbours at edge
        for ( unsigned e = 0; e < ftVec.size(); e++ ) {
            // look only for neighbours in same refinement level
            if ( not ftVec[ e ]->isLeaf() ) {
                // find edge id
                const int iVtx0 = ftVec[ e ]->localVertexIndex( ft->vertex( i ) );
                const int iVtx1 = ftVec[ e ]->localVertexIndex( ft->vertex( next_[ i ] ) );
                const int nEdge = ( iVtx1 == next_[ iVtx0 ] ) ? iVtx0 : previous_[ iVtx0 ];

                // non-manifold neighbor relations
                if ( ftVec[ e ]->neighbor( nEdge ) == ft ) {
                    ftVec[ e ]->child( iVtx0 )->neighbor( nEdge ) = ft->child( i );
                    ftVec[ e ]->child( iVtx1 )->neighbor( nEdge ) = ft->child( next_[ i ] );
                }
                if ( ftVec[ e ] == ft->neighbor( i ) ) {
                    ft->child( i )->neighbor( i )          = ftVec[ e ]->child( iVtx0 );
                    ft->child( next_[ i ] )->neighbor( i ) = ftVec[ e ]->child( iVtx1 );
                }

                // manifold neighbor relations
                if ( ftVec[ e ]->manifoldNeighbor( nEdge ) == ft ) {
                    ftVec[ e ]->child( iVtx0 )->manifoldNeighbor( nEdge ) = ft->child( i );
                    ftVec[ e ]->child( iVtx1 )->manifoldNeighbor( nEdge ) = ft->child( next_[ i ] );
                }
                if ( ftVec[ e ] == ft->manifoldNeighbor( i ) ) {
                    ft->child( i )->manifoldNeighbor( i )          = ftVec[ e ]->child( iVtx0 );
                    ft->child( next_[ i ] )->manifoldNeighbor( i ) = ftVec[ e ]->child( iVtx1 );
                }

            }
        }
    }

    return;
}

