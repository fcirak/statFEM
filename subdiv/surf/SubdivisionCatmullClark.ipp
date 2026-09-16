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
template< typename FT >
const std::array< int, corlib::ShapeTraits< corlib::QUADRILATERAL >::numVertices > 
subdiv::surf::Subdivision< subdiv::surf::CATMULL_CLARK, FT >::next_ =
    subdiv::surf::ShapeProp< corlib::QUADRILATERAL >::next;

//! Local index of previous vertex
template< typename FT >
const std::array< int, corlib::ShapeTraits< corlib::QUADRILATERAL >::numVertices >
subdiv::surf::Subdivision< subdiv::surf::CATMULL_CLARK, FT >::previous_ =
    subdiv::surf::ShapeProp< corlib::QUADRILATERAL >::previous;

//------------------------------------------------------------------------------
// refine facet in given mesh
template< typename FT > 
void subdiv::surf::Subdivision< subdiv::surf::CATMULL_CLARK, FT >::refineFacet( FT * ft,
                                                                                const int level, 
                                                                                Point & posNew )
{
    // compute average of facet vertex co-ordinates
    //
    //       1/4.........1/4
    //        :           :
    //        :     *     :
    //        :           :
    //       1/4.........1/4
    //
    for (unsigned i = 0; i < numVertices; i++)
        posNew += ft->vertex( i )->point( level );
    posNew /= static_cast< double >( numVertices );

    return;
}

//------------------------------------------------------------------------------
template< typename FT >
void subdiv::surf::Subdivision< subdiv::surf::CATMULL_CLARK, FT >::refineEdge( FT * ft,
                                                                               const int iedge,
                                                                               const int level, 
                                                                               Point & posNew )
{
    Vertex * v0 = ft->vertex( iedge );
    Vertex * v1 = ft->vertex( next_[ iedge ] );
    const Point & pos0 = v0->point( level );
    const Point & pos1 = v1->point( level );

    // neigbour facet
    FT * fn = ft->neighbor( iedge );

    // check if border vertex or creased
    if ( ft->getEdgeTag( iedge ) == subdiv::surf::EDGE_CREASE ) {
        // border edge or creased edge
        //
        //       1/2----*----1/2
        //
        posNew += pos0;
        posNew += pos1;
        posNew /= 2.;
    }
    else if ( ft->getEdgeTag( iedge ) == subdiv::surf::EDGE_SMOOTHCREASE ) {
        // "smooth crease" averaged crease and non-singular edge
        //
        //        0...........0
        //        :       (2) :
        //        :    1/4    :
        //    (0) :           : (1)        (0)             (1)
        //       1/4----*----1/4      +      1/2----*----1/2
        //        :           :
        //        :    1/4    :
        //        :     (3)   :
        //        0...........0
        //
        assert( not ft->isManifold( iedge ) );

        // collect all neighbour facets at the edge
        std::vector< FT * > eNeighbors;
        ft->neighborsAtEdge( iedge, eNeighbors );
        assert( eNeighbors.size( ) > 0 );
        eNeighbors.insert( eNeighbors.begin(), ft );

        // loop all common neighbors at edge
        // to find average new co-ordinate
        unsigned numEnvironFacets = 0;
        for ( unsigned n = 0; n < eNeighbors.size( ); ++n ) {
            FT * nFacet = eNeighbors[ n ];
            const unsigned nEdgeId = nFacet->localEdgeIndex( v0, v1 );
            FT * mFacet = nFacet->manifoldNeighbor( nEdgeId );
            if ( mFacet ) {             // like non-singular edge
                Point pos2 = Vertex::makePoint( 0. );
                this->refineFacet( nFacet, level, pos2 );
                Point pos3 = Vertex::makePoint( 0. );
                this->refineFacet( mFacet, level, pos3 );
                posNew += .25*pos0 + .25*pos1 + .25*pos2 + .25*pos3;
                numEnvironFacets += 1;
            }
            else {                     // like crease edge
                posNew += .5*pos0 + .5*pos1;
                numEnvironFacets += 1;
            }
        }
        // average
        FTL_VERIFY( numEnvironFacets > 0 );
        posNew /= static_cast< double >( numEnvironFacets );
    }
    else {
        // ordinary edge
        //
        //        0...........0
        //        :           :
        //        :    1/4    :
        //        :           :
        //       1/4----*----1/4
        //        :           :
        //        :    1/4    :
        //        :           :
        //        0...........0
        //
        Point posFace = Vertex::makePoint( 0. );
        Point posFaceN = Vertex::makePoint( 0. );
        this->refineFacet( ft, level, posFace );
        this->refineFacet( fn, level, posFaceN );
        
        posNew += pos0 + pos1 + posFace + posFaceN;
        posNew /= 4.;
    }

    return;
}

//------------------------------------------------------------------------------
template< typename FT >
void subdiv::surf::Subdivision< subdiv::surf::CATMULL_CLARK, FT >::refineVertex( FT * ft,
                                                                                 const int ivtx,
                                                                                 const int level, 
                                                                                 Point & posNew )
{
    // pointer to vertex
    Vertex * vc = ft->vertex( ivtx );
    std::vector< Vertex * > vtxs;

    if ( vc->tag() == subdiv::surf::VERTEX_CORNER or
         vc->tag() == subdiv::surf::VERTEX_DART ) {
        // corner vertex: no change in position
        posNew = ft->vertex( ivtx )->point( level );
    }
    else if ( vc->tag() == subdiv::surf::VERTEX_CREASE ) {
        // vertex on crease
        //
        //         1/8---------3/4---------1/8
        //

        // retrieve one-ring of vertex with edge information
        typedef std::pair< enum subdiv::surf::edgeTag,
                           enum subdiv::surf::edgeCategory > PairE;
        typedef std::map< Vertex *, PairE >                MapVE;
        typedef typename MapVE::iterator                   MapVEIter;
        typedef std::vector< FT * >                        VecFTPtr;

        MapVE oneRingV;
        VecFTPtr oneRingF;
        subdiv::surf::collectOneRingVertices< FT >( ft, ivtx, oneRingV, oneRingF );

        // check if "smooth crease"
        vtxs.clear( );
        subdiv::surf::extractByEdgeTag( oneRingV, subdiv::surf::EDGE_SMOOTHCREASE, vtxs );
        if ( vtxs.size( ) > 0 ) {
            FTL_VERIFY( false );
        }
        
        // act acc to edge type
        vtxs.clear( );
        subdiv::surf::extractByEdgeCategory( oneRingV, subdiv::surf::ECAT_NONMANIFOLD, vtxs );
        if ( vtxs.size( ) != 2 ) {
            vtxs.clear( );
            subdiv::surf::extractByEdgeCategory( oneRingV, subdiv::surf::ECAT_BOUNDARY, vtxs );
            if ( vtxs.size( ) != 2 ) {
                vtxs.clear( );
                subdiv::surf::extractByEdgeTag( oneRingV, subdiv::surf::EDGE_CREASE, vtxs );
                if ( vtxs.size( ) != 2 ) {
                    vtxs.clear( );
                    // this last resort makes the vertex a fixed/no-subdiv/corner vertex
                    vtxs.push_back( vc );
                    vtxs.push_back( vc );
                }
            }
        }

        assert( vtxs.size() == 2 );
        
        posNew += vtxs[ 0 ]->point( level );
        posNew += vtxs[ 1 ]->point( level );    
        posNew *= 0.125;
        
        Point p0( vc->point( level ) );
        p0 *= 0.75;   
        posNew += p0;        
    }
    else if ( vc->tag() == subdiv::surf::VERTEX_SMOOTHCORNER ) {
        // "smooth corner" vertex: no change in position
        //
        //         1/8----3/4----1/8    +    1
        //

        // retrieve one-ring of vertex with edge information
        typedef std::pair< enum subdiv::surf::edgeTag,
                           enum subdiv::surf::edgeCategory >       PairE;
        typedef std::map< Vertex *, PairE >                        MapVE;
        typedef typename MapVE::iterator                           MapVEIter;
        typedef std::vector< FT * >                                VecFTPtr;

        MapVE oneRingV;
        VecFTPtr oneRingF;
        subdiv::surf::collectOneRingVertices< FT >( ft, ivtx, oneRingV, oneRingF );

#if DEBUG
        // check if vertex has a "smooth crease" attached
        vtxs.clear( );
        subdiv::surf::extractByEdgeTag( oneRingV, subdiv::surf::EDGE_SMOOTHCREASE, vtxs );
        assert( vtxs.size( ) > 0 );
#endif

        // loop one-ring of facets and sum up new position
        unsigned numEnvironFacets = 0;
        for ( unsigned i = 0; i < oneRingF.size( ); ++i ) {
            FT * fn = oneRingF[ i ];
            const unsigned vId = fn->localVertexIndex( vc );
            const typename FT::EdgeTag eTag0 = fn->getEdgeTag( previous_[ vId ] );
            const typename FT::EdgeTag eTag1 = fn->getEdgeTag( vId );
            if ( ( ( eTag0 == subdiv::surf::EDGE_SMOOTHCREASE and
                     eTag1 == subdiv::surf::EDGE_SMOOTHCREASE ) or
                   ( eTag0 == subdiv::surf::EDGE_SMOOTHCREASE and
                     eTag1 == subdiv::surf::EDGE_CREASE ) or
                   ( eTag0 == subdiv::surf::EDGE_CREASE and
                     eTag1 == subdiv::surf::EDGE_SMOOTHCREASE ) ) and
                 ( ( fn->manifoldNeighbor( previous_[ vId ] ) or
                     fn->manifoldNeighbor( vId ) ) ) ) {            // "smooth corner"
                vtxs.clear( );
                VecFTPtr fcts;
                subdiv::surf::collectManifoldOneRingVertices< FT >( fn, vId,
                                                                    &FT::manifoldNeighbor,
                                                                    vtxs, fcts );
                assert( vtxs.size() >= 2 );
                posNew += .125 * (*vtxs.begin( ))->point( level );
                posNew += .125 * (*(--vtxs.end( )))->point( level );
                posNew += .75 * vc->point( level );
                numEnvironFacets += 1;
            }
            else {                                             // true corner
                // corner vertex: no change in position
                posNew += ft->vertex( ivtx )->point( level );
                numEnvironFacets += 1;
            }
        }
            
        // average new position
        FTL_VERIFY( numEnvironFacets > 0 );
        posNew /= static_cast< double >( numEnvironFacets );

    }
    else if ( vc->tag() == subdiv::surf::VERTEX_SMOOTHCREASE ) {
        // "smooth crease" vertex
        //
        // mask with n = number of valence (shown is n=4)
        //
        //       0............0............0
        //       :            |            :
        //       :  1/n^2   2/n^2   1/n^2  :
        //       :            |            :
        //       0--2/n^2--(n-3)/n--2/n^2--0       +       1/8------------3/4----------1/8
        //       :            |            :
        //       :  1/n^2   2/n^2   1/n^2  :
        //       :            |            :
        //       0............0............0
        //
        // retrieve one-ring of vertex with edge information
        typedef std::pair< enum subdiv::surf::edgeTag,
                           enum subdiv::surf::edgeCategory >       PairE;
        typedef std::map< Vertex *, PairE >                        MapVE;
        typedef typename MapVE::iterator                           MapVEIter;
        typedef std::vector< FT * >                                VecFTPtr;

        MapVE oneRingV;
        VecFTPtr oneRingF;
        subdiv::surf::collectOneRingVertices< FT >( ft, ivtx, oneRingV, oneRingF );

        // check if vertex surrounded by (smooth) crease edges
        vtxs.clear( );
        subdiv::surf::extractByEdgeTag( oneRingV, subdiv::surf::EDGE_SMOOTHCREASE, vtxs );
        assert( vtxs.size( ) >= 2 );

        // loop one-ring of facets and sum up new position
        unsigned numEnvironFacets = 0;
        for ( unsigned i = 0; i < oneRingF.size( ); ++i ) {
            FT * fn = oneRingF[ i ];
            const unsigned vId = fn->localVertexIndex( vc );
            if ( fn->getEdgeTag( previous_[ vId ] ) == subdiv::surf::EDGE_SMOOTHCREASE or
                 fn->getEdgeTag( vId ) == subdiv::surf::EDGE_SMOOTHCREASE ) {
                if ( fn->manifoldNeighbor( previous_[ vId ] ) and
                     fn->manifoldNeighbor( vId ) ) {            // "smooth crease"

                    // co-ordinate at vertex
                    Point pc = vc->point( level );
                    // sum of edge points attached to vertex
                    Point pe = Vertex::makePoint( 0. );
                    // sum of face points in 1-neighbourhood of vertex
                    Point pf = Vertex::makePoint( 0. );

                    // collect 1-ring
                    vtxs.clear( );
                    VecFTPtr fcts;
                    subdiv::surf::collectManifoldOneRingVertices< FT >( fn, vId,
                                                                        &FT::manifoldNeighbor,
                                                                        vtxs, fcts );
                    FTL_VERIFY( vtxs.size() == fcts.size() );

                    // compute new position
                    for ( unsigned i = 0; i < vtxs.size(); ++i ) {
                        // edge points
                        pe += .5 * ( pc + vtxs[ i ]->point( level ) );
                        // face points
                        Point pfi = Vertex::makePoint( 0. );
                        this->refineFacet( fcts[ i ], level, pfi );
                        pf += pfi;
                    }
                    const double valence = static_cast< double >( vtxs.size() );
                    posNew += ( pf/valence + 2.*pe/valence + (valence-3.)*pc ) / valence;
                    numEnvironFacets += 1;
                }
                else if ( fn->manifoldNeighbor( previous_[ vId ] ) or
                          fn->manifoldNeighbor( vId ) ) {       // true crease
                    vtxs.clear( );
                    VecFTPtr fcts;
                    subdiv::surf::collectManifoldOneRingVertices< FT >( fn, vId,
                                                                        &FT::manifoldNeighbor,
                                                                        vtxs, fcts );
                    assert( vtxs.size() >= 2 );
                    posNew += .125 * (*vtxs.begin( ))->point( level );
                    posNew += .125 * (*(--vtxs.end( )))->point( level );
                    posNew += .75 * vc->point( level );
                    numEnvironFacets += 1;
                }
                else {
                    // there should be at least one manifold neigbour
                    FTL_VERIFY_DESCRIPTIVE( false,
                                            "a vertex on a smooth crease must have"
                                            " two crease edges attached.\n" );
                }
            }
        }
            
        // average new position
        FTL_VERIFY( numEnvironFacets > 0 );
        posNew /= static_cast< double >( numEnvironFacets );
        
    }
    else { 
        // field vertex
        //
        // mask with n = number of valence (shown is n=4)
        //
        //       0............0............0
        //       :            |            :
        //       :  1/n^2   2/n^2   1/n^2  :
        //       :            |            :
        //       0--2/n^2--(n-3)/n--2/n^2--0
        //       :            |            :
        //       :  1/n^2   2/n^2   1/n^2  :
        //       :            |            :
        //       0............0............0
        //
        
        // collect 1-ring
        std::vector< FT * > fcts;
        subdiv::surf::collectManifoldOneRingVertices< FT >( ft, ivtx, &FT::neighbor,
                                                            vtxs, fcts );
        FTL_VERIFY( vtxs.size() == fcts.size() );  // valence

        // co-ordinate at vertex
        Point pc = vc->point( level );
        // sum of edge points attached to vertex
        Point pe = Vertex::makePoint( 0. );
        // sum of face points in 1-neighbourhood of vertex
        Point pf = Vertex::makePoint( 0. );

        // compute new position
        for ( unsigned i = 0; i < vtxs.size(); ++i ) {
            // edge points
            pe += .5 * ( pc + vtxs[ i ]->point( level ) );
            // face points
            Point pfi = Vertex::makePoint( 0. );
            this->refineFacet( fcts[ i ], level, pfi );
            pf += pfi;
        }
        const double valence = static_cast< double >( vtxs.size() );
        posNew += ( pf/valence + 2.*pe/valence + (valence-3.)*pc ) / valence;

    }   
    return;
}

//------------------------------------------------------------------------------
template< typename FT >
void subdiv::surf::Subdivision< subdiv::surf::CATMULL_CLARK, FT >::computeTangents( FT * ft,
                                                                                    const int ivtx,
                                                                                    const int level,
                                                                                    MapVPtrD & coeff0,
                                                                                    MapVPtrD & coeff1 )
{
    // centre vertex
    Vertex * vc = ft->vertex( ivtx );

    // retrieve one-ring of vertex with edge information
    typedef std::vector< Vertex * >                    VecVPtr;
    typedef typename VecVPtr::iterator                 VecVPtrIter;
    typedef std::vector< FT * >                        VecFTPtr;

    VecVPtr vtxsE;  // edge vertices on manifold
    VecFTPtr fcts;  // (open/closed) ring of facets on manifold
    subdiv::surf::collectManifoldOneRingVertices< FT >( ft, ivtx,
                                                        &FT::manifoldNeighbor,
                                                        vtxsE, fcts, true );

    // get face vertices on manifold --- maintain ordering
    VecVPtr vtxsF;
    vtxsF.reserve( vtxsF.size( ) );
    for ( unsigned f = 0; f < fcts.size( ); ++f ) {
        for ( unsigned v = 0; v < FT::numVertices; ++v ) {
            if ( not corlib::contains( vtxsE.begin(), vtxsE.end(), fcts[ f ]->vertex( v ) ) and
                 not ( fcts[ f ]->vertex( v ) == vc ) ) {
                vtxsF.push_back( fcts[ f ]->vertex( v ) );
            }
        }
    }
    FTL_VERIFY( vtxsF.size( ) == fcts.size( ) );

    // number of adjacent polygons
    const double val = static_cast< double >( fcts.size( ) );

    // compute tangents
    if ( vc->tag() == subdiv::surf::VERTEX_CORNER ) {
        FTL_VERIFY( vtxsE.size( ) == ( fcts.size( ) + 1 ) );
        FTL_VERIFY( vtxsF.size( ) == fcts.size( ) );

        // check if first and last vertex are on one line
        // to prevent collinear tangents
        Vertex * vFirst = *(vtxsE.begin( ));
        Vertex * vLast = *(--vtxsE.end( ));
        typedef typename Vertex::VecDim VecDim;
        const VecDim pCentre = vc->giveCoordinates( );
        const VecDim pFirst = vFirst->giveCoordinates( );
        const VecDim pLast = vLast->giveCoordinates( );
        const double dirSin = ( corlib::cross_prod( pFirst - pCentre,
                                                    pLast - pCentre ) ).norm( );
        if ( corlib::fuzzyEqual( dirSin, 0. ) ) {
            vLast = ft->vertex( next_[ ivtx ] );
            vFirst = ft->vertex( previous_[ ivtx ] );
        }

        // centre vertex
        coeff1[ vc ] = -1.;
        // edge vertices
        coeff1[ vLast ] = 1.;
        // face vertex
        ;  // do nothing

        // centre vertex
        coeff0[ vc ] = -1.;
        // edge vertices
        coeff0[ vFirst ] = 1.;
        // face vertex
        ;  // do nothing
    }
    else if ( vc->tag() == subdiv::surf::VERTEX_DART ) {
        FTL_VERIFY( vtxsF.size( ) == fcts.size( ) );
        if ( vtxsE.size( ) == ( fcts.size( ) + 1 ) and
             *(vtxsE.begin( )) == *(--vtxsE.end( )) ) {  // same as VERTEX_NOTAG

            // PLEASE NOTE: One of tangents at the dart vertex is not unique.
            //              Here, we simple deliver the tangents of an ordinary vertex
            //              ignoring the dart status of the vertex. This is not perfect.
            //              However, simply using the corner rule leads to two identical
            //              tangent vectors, which is not good either.

            // factors
            const double th = 2.*M_PI/val;
            const double lam = 5./16.+
                (std::cos(th)+std::cos(th/2.)*std::sqrt(9.+std::cos(2.*th)))/16.;  // ????

            // centre vertex
            ;  // do nothing

            // edge vertices
            for ( unsigned i = 0; i < fcts.size( ); ++i ) {
                Vertex * v = vtxsE[ i ];
                coeff1[ v ] = 4.*std::sin(i*th);
                coeff0[ v ] = 4.*std::cos(i*th);
            }
            
            // face vertices
            for ( unsigned i = 0; i < vtxsF.size( ); ++i ) {
                Vertex * v = vtxsF[ i ];
                coeff1[ v ] = (std::sin(i*th)+std::sin((i+1.)*th))/(4.*lam-1.);
                coeff0[ v ] = (std::cos(i*th)+std::cos((i+1.)*th))/(4.*lam-1.);
            }
        }
        else {  // same as VERTEX_CORNER

            // centre vertex
            coeff0[ vc ] = -1.;
            // edge vertices
            coeff0[ *(--vtxsE.end( )) ] = 1.;
            // face vertex
            ;  // do nothing

            // centre vertex
            coeff1[ vc ] = -1.;
            // edge vertices
            coeff1[ *(vtxsE.begin( )) ] = 1.;
            // face vertex
            ;  // do nothing
        }
    }
    else if ( vc->tag() == subdiv::surf::VERTEX_CREASE ) {
        if ( val == 1 ) {
            FTL_VERIFY( vtxsE.size( ) == 2 );
            FTL_VERIFY( vtxsF.size( ) == 1 );
            
            // centre vertex
            coeff1[ vc ] = 6.;
            // edge vertices
            coeff1[ *(vtxsE.begin()) ] = -3.;
            coeff1[ *(--vtxsE.end()) ] = -3.;
            // face vertex
            ;  // do nothing

            // centre vertex
            ;  // do nothing 
            // edge vertices
            coeff0[ *(vtxsE.begin()) ] = -1.;
            coeff0[ *(--vtxsE.end()) ] = 1.;
            // face vertex
            ;  // do nothing 
        }
        else {
            FTL_VERIFY_DESCRIPTIVE( vtxsE.size( ) == ( fcts.size( ) + 1 ),
                                    "with vtxsE.size( )=%d, fcts.size( )=%d\n",
                                    vtxsE.size( ), fcts.size( ) );
            FTL_VERIFY( vtxsF.size( ) == fcts.size( ) );

            // factors
            const double th = M_PI/val;
            const double r = (std::cos(th)+1.)/(val*std::sin(th)*(3.+std::cos(th)));

            // centre vertex
            coeff1[ vc ] = 4.*r*(std::cos(th)-1.);
            // edge vertices
            coeff1[ *(vtxsE.begin( )) ] = -r*(1.+2.*std::cos(th));
            for ( VecVPtrIter v = ++vtxsE.begin( ); v != --vtxsE.end( ); ++v ) {
                const unsigned i = std::distance( vtxsE.begin( ), v );
                coeff1[ *v ] = 4.*std::sin(i*th)/(3.+std::cos(th))/val;
            }
            coeff1[ *(--vtxsE.end( )) ] = -r*(1.+2.*std::cos(th));
            // face vertices
            for ( VecVPtrIter v = vtxsF.begin( ); v != vtxsF.end( ); ++v ) {
                const unsigned i = std::distance( vtxsF.begin( ), v );
                coeff1[ *v ] = (std::sin(i*th)+std::sin((i+1.)*th))/(3.+std::cos(th))/val;  // corrected reference, factor 4 too much
            }
            
            // centre vertex
            ;  // do nothing
            // edge vertices
            coeff0[ *(vtxsE.begin( )) ] = .5;
            coeff0[ *(--vtxsE.end( )) ] = -.5;
            // face vertices
            ;  // do nothing
        }
    }
    else if ( vc->tag() == subdiv::surf::VERTEX_NOTAG ) {
        FTL_VERIFY( vtxsE.size( ) == fcts.size( ) );
        FTL_VERIFY( vtxsF.size( ) == fcts.size( ) );

        // factors
        const double th = 2.*M_PI/val;
        const double lam = 5./16.+
            (std::cos(th)+std::cos(th/2.)*std::sqrt(9.+std::cos(2.*th)))/16.;  // ????

        // centre vertex
        ;  // do nothing

        // edge vertices
        for ( unsigned i = 0; i < vtxsE.size( ); ++i ) {
            Vertex * v = vtxsE[ i ];
            coeff1[ v ] = 4.*std::sin(i*th);
            coeff0[ v ] = 4.*std::cos(i*th);
        }

        // face vertices
        for ( unsigned i = 0; i < vtxsF.size( ); ++i ) {
            Vertex * v = vtxsF[ i ];
            coeff1[ v ] = (std::sin(i*th)+std::sin((i+1.)*th))/(4.*lam-1.);
            coeff0[ v ] = (std::cos(i*th)+std::cos((i+1.)*th))/(4.*lam-1.);
        }
    }
    else {
        FTL_VERIFY_DESCRIPTIVE( false,
                                "Not implemented vertexTag=%d\n",
                                vc->tag() );
    }

    return;
}

//------------------------------------------------------------------------------
template< typename FT >
void subdiv::surf::Subdivision< subdiv::surf::CATMULL_CLARK, FT >::computeLimit( FT * ft,
                                                                                 const int ivtx,
                                                                                 const int level,
                                                                                 MapVPtrD & coeff )
{
    // centre vertex
    Vertex * vc = ft->vertex( ivtx );

    // retrieve one-ring of vertex with edge information
    typedef std::vector< Vertex * >                    VecVPtr;
    typedef typename VecVPtr::iterator                 VecVPtrIter;
    typedef std::vector< FT * >                        VecFTPtr;

    VecVPtr vtxsE;  // edge vertices on manifold
    VecFTPtr fcts;  // (open/closed) ring of facets on manifold
    subdiv::surf::collectManifoldOneRingVertices< FT >( ft, ivtx,
                                                        &FT::manifoldNeighbor,
                                                        vtxsE, fcts, true );

    // get face vertices on manifold --- maintain ordering
    VecVPtr vtxsF;
    vtxsF.reserve( vtxsF.size( ) );
    for ( unsigned f = 0; f < fcts.size( ); ++f ) {
        for ( unsigned v = 0; v < FT::numVertices; ++v ) {
            if ( not corlib::contains( vtxsE.begin(), vtxsE.end(), fcts[ f ]->vertex( v ) ) and
                 not ( fcts[ f ]->vertex( v ) == vc ) ) {
                vtxsF.push_back( fcts[ f ]->vertex( v ) );
            }
        }
    }
    FTL_VERIFY( vtxsF.size( ) == fcts.size( ) );

    // number of adjacent polygons
    const double val = static_cast< double >( fcts.size( ) );

    // compute tangents
    if ( vc->tag() == subdiv::surf::VERTEX_CORNER ) {
        FTL_VERIFY( vtxsE.size( ) == ( fcts.size( ) + 1 ) );
        FTL_VERIFY( vtxsF.size( ) == fcts.size( ) );

        // centre vertex
        coeff[ vc ] = 1.;
        // edge vertices
        ;  // do nothing
        // face vertex
        ;  // do nothing
    }
    else if ( vc->tag() == subdiv::surf::VERTEX_DART ) {
        FTL_VERIFY( vtxsF.size( ) == fcts.size( ) );
        if ( vtxsE.size( ) == ( fcts.size( ) + 1 ) and
             *(vtxsE.begin( )) == *(--vtxsE.end( )) ) {  // same as VERTEX_NOTAG

            // PLEASE NOTE: One of tangents at the dart vertex is not unique.
            //              Here, we simple deliver the tangents of an ordinary vertex
            //              ignoring the dart status of the vertex. This is not perfect.
            //              However, simply using the corner rule leads to two identical
            //              tangent vectors, which is not good either.

            // factors
            const double factC = val/(val+5.);
            const double factE = 4./(val*(val+5.));
            const double factF = 1./(val*(val+5.));

            // centre vertex
            coeff[ vc ] = factC;
            // edge vertices
            for ( unsigned i = 0; i < vtxsE.size( ); ++i )
                coeff[ vtxsE[ i ] ] = factE;
            // face vertices
            for ( unsigned i = 0; i < vtxsF.size( ); ++i )
                coeff[ vtxsF[ i ] ] = factF;
        }
        else {  // same as VERTEX_CORNER

            // centre vertex
            coeff[ vc ] = 1.;
            // edge vertices
            ;  // do nothing
            // face vertex
            ;  // do nothing
        }
    }
    else if ( vc->tag() == subdiv::surf::VERTEX_CREASE ) {
        FTL_VERIFY_DESCRIPTIVE( vtxsE.size( ) == ( fcts.size( ) + 1 ),
                                "with vtxsE.size()=%d, fcts.size()=%d\n",
                                vtxsE.size( ), fcts.size( ) );
        FTL_VERIFY( vtxsF.size( ) == fcts.size( ) );

        // centre vertex
        coeff[ vc ] = 2./3.;
        // edge vertices
        coeff[ *(vtxsE.begin()) ] = 1./6.;
        coeff[ *(--vtxsE.end()) ] = 1./6.;
        // face vertex
        ;  // do nothing
    }
    else if ( vc->tag() == subdiv::surf::VERTEX_NOTAG ) {
        FTL_VERIFY( vtxsE.size( ) == fcts.size( ) );
        FTL_VERIFY( vtxsF.size( ) == fcts.size( ) );

        // factors
        const double factC = val/(val+5.);
        const double factE = 4./(val*(val+5.));
        const double factF = 1./(val*(val+5.));

        // centre vertex
        coeff[ vc ] = factC;
        // edge vertices
        for ( unsigned i = 0; i < vtxsE.size( ); ++i )
            coeff[ vtxsE[ i ] ] = factE;
        // face vertices
        for ( unsigned i = 0; i < vtxsF.size( ); ++i )
            coeff[ vtxsF[ i ] ] = factF;
    }
    else {
        FTL_VERIFY_DESCRIPTIVE( false,
                                "Not implemented vertexTag=%d\n",
                                vc->tag() );
    }

    return;
}
