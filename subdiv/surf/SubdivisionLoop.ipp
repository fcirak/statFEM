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

#include <subdiv/surf/collect.hpp>

//------------------------------------------------------------------------------
//! Local index of next vertex
template< typename FT >
const std::array< int, 3 >
subdiv::surf::Subdivision< subdiv::surf::LOOP, FT >::next_ =
    ShapeProp< corlib::TRIANGLE >::next;

//! Local index of previous vertex
template< typename FT >
const std::array< int, 3 >
subdiv::surf::Subdivision< subdiv::surf::LOOP, FT >::previous_ =
    ShapeProp< corlib::TRIANGLE >::previous;

//------------------------------------------------------------------------------
/// refine facet in given mesh
template< typename FT >
void subdiv::surf::Subdivision< subdiv::surf::LOOP,FT >::refineFacet( FT * ft,
                                                                      const int level,
                                                                      Point & posNew )
{
    // do nothing
    return;
}

//------------------------------------------------------------------------------
/// refine edge in given mesh
template< typename FT >
void subdiv::surf::Subdivision< subdiv::surf::LOOP, FT >::refineEdge( FT * ft,
                                                                      const int iedge,
                                                                      const int level,
                                                                      Point & posNew )
{
    Vertex * v0 = ft->vertex( iedge );
    Vertex * v1 = ft->vertex( next_[ iedge ] );
    const Point & pos0 = v0->point( level );
    const Point & pos1 = v1->point( level );

    if ( ft->getEdgeTag( iedge ) == subdiv::surf::EDGE_CREASE ) {
        // border edge or creased edge
        //
        // pos0  1/2----*----1/2  pos1
        //
        posNew += pos0;
        posNew += pos1;
        posNew /= 2.;
    }
    else if ( ft->getEdgeTag( iedge ) == subdiv::surf::EDGE_SMOOTHCREASE ) {
        // "smooth crease" averaged crease and non-singular edge
        //
        //              1/8 (2)
        //              / `
        //             /   `
        //            /     `
        //       (0) /       ` (1)        (0)           (1)
        //        3/8----*----3/8    +     1/2----*----1/2
        //           `       /
        //            `     /
        //             `   /
        //              ` /
        //              1/8 (3)
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
            if ( mFacet ) {
                const unsigned mEdgeId = mFacet->localEdgeIndex( v0, v1 );
                // like ordinary edge
                const Point & pos2 = nFacet->vertex( previous_[ nEdgeId ] )->point( level );
                const Point & pos3 = mFacet->vertex( previous_[ mEdgeId ] )->point( level );
                posNew += .375*pos0 + .375*pos1 + .125*pos2 + .125*pos3;
                numEnvironFacets += 1;
            }
            else {
                // like crease edge
                posNew += .5*pos0 + .5*pos1;
                numEnvironFacets += 1;
            }
        }
        // average
        FTL_VERIFY( numEnvironFacets > 0 );
        posNew /= static_cast< double >( numEnvironFacets );
    }
    else {
        // non-singular edge
        //
        //              1/8   pos2
        //              / `
        //             /   `
        //            /     `
        //           /       `
        // pos0   3/8----*----3/8   pos1
        //           `       /
        //            `     /
        //             `   /
        //              ` /
        //              1/8   pos3
        //
        const Point & pos2 = ft->vertex( previous_[ iedge ] )->point( level );

        FT * fn = ft->neighbor( iedge );
        const int iv3 = next_[ fn->localVertexIndex( v0 ) ];
        const Point & pos3 = fn->vertex( iv3 )->point( level );

#if 1
        posNew += pos0 + pos1;
        posNew *= 0.375;
        posNew += 0.125*pos2 + 0.125*pos3;
#else
        posNew += pos0;
        posNew += pos1;
        posNew *= 0.375;

        Point tmp( pos2 );
        tmp += pos3;
        tmp *= 0.125;

        posNew += tmp;
#endif
    }
    return;
}


//------------------------------------------------------------------------------
template< typename FT >
void subdiv::surf::Subdivision< subdiv::surf::LOOP, FT >::refineVertex( FT * ft,
                                                                        const int ivtx,
                                                                        const int level,
                                                                        Point & posNew )
{
    Vertex * vc = ft->vertex( ivtx );
    std::vector< Vertex * > vtxs;

    if ( vc->tag() == subdiv::surf::VERTEX_CORNER or 
         vc->tag() == subdiv::surf::VERTEX_DART ) {
        // corner vertex: no change in position
        //
        //                 1
        //
        posNew = ft->vertex( ivtx )->point( level );
    }
    else if ( vc->tag() == subdiv::surf::VERTEX_CREASE ) {
        // vertex on crease
        //
        //         1/8----3/4----1/8
        //

        // retrieve one-ring of vertex with edge information
        typedef std::pair< enum subdiv::surf::edgeTag,
                           enum subdiv::surf::edgeCategory > PairE;
        typedef std::map< Vertex *, PairE >                  MapVE;
        typedef typename MapVE::iterator                     MapVEIter;
        typedef std::vector< FT * >                          VecFTPtr;

        MapVE oneRingV;
        VecFTPtr oneRingF;
        subdiv::surf::collectOneRingVertices< FT >( ft, ivtx, oneRingV, oneRingF );

        // check no "smooth crease"
#ifdef DEBUG
        vtxs.clear( );
        subdiv::surf::extractByEdgeTag( oneRingV, subdiv::surf::EDGE_SMOOTHCREASE, vtxs );
        assert( vtxs.size( ) == 0 );
#endif

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

        Point tmp( vc->point( level ) );
        tmp *= 0.75;
        posNew += tmp;

    }
    else if ( vc->tag() == subdiv::surf::VERTEX_SMOOTHCORNER ) {
        // "smooth corner" vertex: no change in position
        //
        //         1/8----3/4----1/8       +        1
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

        // check if vertex has a "smooth crease" attached
        vtxs.clear( );
        subdiv::surf::extractByEdgeTag( oneRingV, subdiv::surf::EDGE_SMOOTHCREASE, vtxs );
        FTL_VERIFY( vtxs.size( ) > 0 );
        
        // loop one-ring of facets and sum up new position
        unsigned numEnvironFacets = 0;
        for ( unsigned i = 0; i < oneRingF.size( ); ++i ) {
            FT * fn = oneRingF[ i ];
            const unsigned vId = fn->localVertexIndex( vc );
            const typename FT::EdgeTag eTag0 = fn->getEdgeTag( previous_[ vId ] );
            const typename FT::EdgeTag eTag1 = fn->getEdgeTag( vId );
            if ( ( eTag0 == subdiv::surf::EDGE_SMOOTHCREASE and
                   fn->manifoldNeighbor( previous_[ vId ] ) ) or
                 ( eTag1 == subdiv::surf::EDGE_SMOOTHCREASE  and
                   fn->manifoldNeighbor( vId ) ) ) {            // "smooth corner"

                vtxs.clear( );
                VecFTPtr fcts;
                subdiv::surf::collectManifoldOneRingVertices< FT >( fn, vId,
                                                                    &FT::manifoldNeighbor,
                                                                    vtxs, fcts );
                FTL_VERIFY( vtxs.size() >= 2 );

                posNew += .125 * (*vtxs.begin( ))->point( level );
                posNew += .125 * (*(--vtxs.end( )))->point( level );
                posNew += .75 * vc->point( level );
                numEnvironFacets += 1;
            } 
            else if ( eTag0 == subdiv::surf::EDGE_SMOOTHCREASE or
                      eTag1 == subdiv::surf::EDGE_SMOOTHCREASE ) { // true corner
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
        // mask with n = number of valence (shown is n=6)
        //
        //               3/8n----3/8n
        //               / `     / `
        //              /   `   /   `
        //             /     ` /     `
        //           3/8n----5/8-----3/8n    +     1/8----3/4----1/8
        //             `     / `     /
        //              `   /   `   /
        //               ` /     ` /
        //               3/8n----3/8n
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
        FTL_VERIFY( vtxs.size( ) == 2 );

        // loop one-ring of facets and sum up new position
        unsigned numEnvironFacets = 0;
        for ( unsigned i = 0; i < oneRingF.size( ); ++i ) {
            FT * fn = oneRingF[ i ];
            const unsigned vId = fn->localVertexIndex( vc );
            if ( fn->getEdgeTag( previous_[ vId ] ) == subdiv::surf::EDGE_SMOOTHCREASE or
                 fn->getEdgeTag( vId ) == subdiv::surf::EDGE_SMOOTHCREASE ) {
                if ( fn->manifoldNeighbor( previous_[ vId ] ) and
                     fn->manifoldNeighbor( vId ) ) {            // "smooth crease"
                    vtxs.clear( );
                    VecFTPtr fcts;
                    subdiv::surf::collectManifoldOneRingVertices< FT >( fn, vId,
                                                                        &FT::manifoldNeighbor,
                                                                        vtxs, fcts );
                    FTL_VERIFY( vtxs.size() == fcts.size() );
                    const double fact = 0.375 / static_cast< double >( vtxs.size( ) );
                    for ( unsigned j = 0; j < vtxs.size( ); ++j )
                        posNew += fact * vtxs[ j ]->point( level );
                    posNew += 0.625 * vc->point( level );
                    numEnvironFacets += 1;
                }
                else if ( fn->manifoldNeighbor( previous_[ vId ] ) or
                          fn->manifoldNeighbor( vId ) ) {       // true crease
                    vtxs.clear( );
                    VecFTPtr fcts;
                    subdiv::surf::collectManifoldOneRingVertices< FT >( fn, vId,
                                                                        &FT::manifoldNeighbor,
                                                                        vtxs, fcts );
                    FTL_VERIFY( vtxs.size() >= 2 );
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
        // mask with n = number of valence (shown is n=6)
        //
        //               3/8n----3/8n
        //               / `     / `
        //              /   `   /   `
        //             /     ` /     `
        //           3/8n----5/8-----3/8n
        //             `     / `     /
        //              `   /   `   /
        //               ` /     ` /
        //               3/8n----3/8n
        //
        std::vector< FT * > fcts;
        subdiv::surf::collectManifoldOneRingVertices< FT >( ft, ivtx, &FT::neighbor,
                                                            vtxs, fcts );
        const unsigned valence = vtxs.size( );
        for ( unsigned i=0; i<valence; ++i ) {
            posNew += vtxs[ i ]->point( level );
        }
//#define SUBDIVISIONLOOP_USE_LOOPORIGINALFACTORS
#ifdef SUBDIVISIONLOOP_USE_LOOPORIGINALFACTORS
        // (original) Loop's factor
        const double n = static_cast< double >( valence );
        const double a = (40.-std::pow((3.+2.*std::cos(2.*M_PI/n)),2.))/64.;
#else
        // Warren's factor
        const double a = 0.375;  // equal to 3./8.
#endif
        posNew *= ( a / static_cast< double >( valence ) );
        posNew += ( 1. - a ) * vc->point( level );
    }
    return;
}

//------------------------------------------------------------------------------
template< typename FT >
void subdiv::surf::Subdivision< subdiv::surf::LOOP, FT >::computeTangents( FT * ft,
                                                                           const int ivtx,
                                                                           const int level,
                                                                           MapVPtrD & coeff0,
                                                                           MapVPtrD & coeff1 )
{
    // centre vertex
    Vertex * vc = ft->vertex( ivtx );

    // retrieve one-ring of vertex with edge information
    typedef std::vector< Vertex * >                    VecVPtr;
    typedef std::vector< FT * >                        VecFTPtr;

    VecVPtr vtxs;
    VecFTPtr fcts;
    subdiv::surf::collectManifoldOneRingVertices< FT >( ft, ivtx,
                                                        &FT::manifoldNeighbor,
                                                        vtxs, fcts, true );

    // number of adjacent polygons
    const double val = static_cast< double >( fcts.size( ) );

    // compute tangents
    if ( vc->tag() == subdiv::surf::VERTEX_CORNER ) {
        FTL_VERIFY( vtxs.size() == ( fcts.size() + 1 ) );

        // check if first and last vertex are on one line
        // to prevent collinear tangents
        Vertex * vFirst = *(vtxs.begin( ));
        Vertex * vLast = *(--vtxs.end( ));
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

        // centre vertex
        coeff0[ vc ] = -1.;
        // edge vertices
        coeff0[ vFirst ] = 1.;
    }
    else if ( vc->tag() == subdiv::surf::VERTEX_DART ) {
        if ( vtxs.size() == ( fcts.size() + 1 ) and
             *(vtxs.begin()) == *(--vtxs.end()) ) {  // same as VERTEX_NOTAG

            // PLEASE NOTE: One of tangents at the dart vertex is not unique.
            //              Here, we simple deliver the tangents of an ordinary vertex
            //              ignoring the dart status of the vertex. This is not perfect.
            //              However, simply using the corner rule leads to two identical
            //              tangent vectors, which is not good either.

            // centre vertex
            ;  // do nothing

            // edge vertices
            for ( unsigned i = 0; i < fcts.size( ); ++i ) {
                const double theta = 2. * M_PI / val;
                coeff1[ vtxs[ i ] ] = 2.*std::sin(i*theta)/val;
                coeff0[ vtxs[ i ] ] = 2.*std::cos(i*theta)/val;
            }
        }
        else {  // same as VERTEX_CORNER

            // centre vertex
            coeff1[ vc ] = -1.;
            // edge vertices
            coeff1[ *(--vtxs.end()) ] = 1.;
            
            // centre vertex
            coeff0[ vc ] = -1.;
            // edge vertices
            coeff0[ *(vtxs.begin()) ] = 1.;
        }
    }
    else if ( vc->tag() == subdiv::surf::VERTEX_CREASE ) {
        if ( val == 1 ) {
            FTL_VERIFY( vtxs.size( ) == 2 );
            
            // centre vertex
            ;  // do nothing
            // edge vertices
            coeff0[ *(vtxs.begin()) ] = .5;
            coeff0[ *(--vtxs.end()) ] = -.5;

            // centre vertex
            coeff1[ vc ] = 1.;
            // edge vertices
            coeff1[ *(vtxs.begin()) ] = -.5;
            coeff1[ *(--vtxs.end()) ] = -.5;            
        }
        else {
            FTL_VERIFY( vtxs.size( ) == ( fcts.size( ) + 1 ) );

            // factors
            const double th = M_PI/val;
            const double zeta = std::acos(std::cos(th)-1.);
            const double a = .25*(1.+std::cos(th))/(3.*(.5-.25*std::cos(th)));
            const double b = (2./3.-a)/std::cos(val*zeta/2.);
            const double sig1 = std::sin(th)/(1.-std::cos(th));
            const double sig3 = std::cos(val*zeta/2.)*std::sin(th)/(std::cos(zeta)-std::cos(th));

            // centre vertex
            ;  // do nothing
            // edge vertices
            coeff0[ *(vtxs.begin()) ] = .5;
            coeff0[ *(--vtxs.end()) ] = -.5;

            // centre vertex
            coeff1[ vc ] = -2./val*((2./3.-a)*sig1-b*sig3);
            // edge vertices
            coeff1[ *(vtxs.begin()) ] = -2./val*((.5*a+1./6.)*sig1+.5*b*sig3);
            for ( unsigned i = 1; i < fcts.size( ); ++i ) {
                coeff1[ vtxs[ i ] ] = 2./val*std::sin(i*th);
            }
            coeff1[ *(--vtxs.end()) ] = -2./val*((.5*a+1./6.)*sig1+.5*b*sig3);

        }
    }
    else if ( vc->tag() == subdiv::surf::VERTEX_NOTAG ) {
        FTL_VERIFY( vtxs.size( ) == fcts.size( ) );

        // centre vertex
        ;  // do nothing

        // edge vertices
        for ( unsigned i = 0; i < vtxs.size( ); ++i ) {
            const double theta = 2. * M_PI / val;
            coeff1[ vtxs[ i ] ] = 2.*std::sin(i*theta)/val;
            coeff0[ vtxs[ i ] ] = 2.*std::cos(i*theta)/val;
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
void subdiv::surf::Subdivision< subdiv::surf::LOOP, FT >::computeLimit( FT * ft,
                                                                        const int ivtx,
                                                                        const int level,
                                                                        MapVPtrD & coeff )
{
    // centre vertex
    Vertex * vc = ft->vertex( ivtx );

    // retrieve one-ring of vertex with edge information
    typedef std::vector< Vertex * >                    VecVPtr;
    typedef std::vector< FT * >                        VecFTPtr;

    VecVPtr vtxs;
    VecFTPtr fcts;
    subdiv::surf::collectManifoldOneRingVertices< FT >( ft, ivtx,
                                                        &FT::manifoldNeighbor,
                                                        vtxs, fcts, true );

    // number of adjacent polygons
    const double val = static_cast< double >( fcts.size( ) );

    // compute tangents
    if ( vc->tag() == subdiv::surf::VERTEX_CORNER ) {
        //FTL_VERIFY( vtxs.size() == ( fcts.size() + 1 ) );

        // // check if first and last vertex are on one line
        // // to prevent collinear tangents
        // Vertex * vFirst = *(vtxs.begin( ));
        // Vertex * vLast = *(--vtxs.end( ));
        // typedef typename Vertex::VecDim VecDim;
        // const VecDim pCentre = vc->giveCoordinates( );
        // const VecDim pFirst = vFirst->giveCoordinates( );
        // const VecDim pLast = vLast->giveCoordinates( );
        // const double dirSin = ( corlib::cross_prod( pFirst - pCentre,
        //                                             pLast - pCentre ) ).norm( );
        // if ( corlib::fuzzyEqual( dirSin, 0. ) ) {
        //     vLast = ft->vertex( next_[ ivtx ] );
        //     vFirst = ft->vertex( previous_[ ivtx ] );
        // }

        // centre vertex
        coeff[ vc ] = 1.;
        // edge vertices
        ;  // do nothing

    }
    else if ( vc->tag() == subdiv::surf::VERTEX_DART ) {
        if ( vtxs.size() == ( fcts.size() + 1 ) and
             *(vtxs.begin()) == *(--vtxs.end()) ) {  // same as VERTEX_NOTAG

            // PLEASE NOTE: One of tangents at the dart vertex is not unique.
            //              Here, we simple deliver the tangents of an ordinary vertex
            //              ignoring the dart status of the vertex. This is not perfect.
            //              However, simply using the corner rule leads to two identical
            //              tangent vectors, which is not good either.

#ifdef SUBDIVISIONLOOP_USE_LOOPORIGINALFACTORS
            const double beta = 3./(11.-(24.+std::pow((3.+2.*std::cos(2.*M_PI/n)),2.))/8.);
#else
            const double beta = 1./2.;  // Warren's factor
#endif
            // centre vertex
            coeff[ vc ] = beta;
            // edge vertices
            for ( unsigned i = 0; i < vtxs.size( ); ++i ) {
                coeff[ vtxs[ i ] ] = (1.-beta)/val;
            }

        }
        else {  // same as VERTEX_CORNER

            // centre vertex
            coeff[ vc ] = 1.;
            // edge vertices
            ;  // do nothing

        }
    }
    else if ( vc->tag() == subdiv::surf::VERTEX_CREASE ) {
        // centre vertex
        coeff[ vc ] = 2./3.;
        // edge vertices
        coeff[ *(vtxs.begin()) ] = 1./6.;
        coeff[ *(--vtxs.end()) ] = 1./6.;            
    }
    else if ( vc->tag() == subdiv::surf::VERTEX_NOTAG ) {
        FTL_VERIFY( vtxs.size( ) == fcts.size( ) );

        // constant
#ifdef SUBDIVISIONLOOP_USE_LOOPORIGINALFACTORS
        const double beta = 3./(11.-(24.+std::pow((3.+2.*std::cos(2.*M_PI/n)),2.))/8.);
#else
        const double beta = 1./2.;   // Warren's factor
#endif
        // centre vertex
        coeff[ vc ] = beta;

        // edge vertices
        for ( unsigned i = 0; i < vtxs.size( ); ++i ) {
            coeff[ vtxs[ i ] ] = (1.-beta)/val;
        }
    }
    else {
        FTL_VERIFY_DESCRIPTIVE( false,
                                "Not implemented vertexTag=%d\n",
                                vc->tag() );
    }

    return;
}
