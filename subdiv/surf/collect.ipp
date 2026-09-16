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

#include<functional>

#ifndef subdiv_surf_onering_h
#define subdiv_surf_onering_h

//==============================================================================
template< typename FACET >
void subdiv::surf::collectManifoldOneRingFacets( FACET * facet,
                                               std::vector< FACET * > & oneRingF,
                                               std::set< typename FACET::Vertex * > & oneRingV,
                                               bool & closedRing )
{

    typedef FACET                                      Facet;
    typedef std::vector< Facet * >                     VecFPtr;
    typedef typename VecFPtr::iterator                 VecFPtrIter;
    typedef std::set< Facet * >                        SetFPtr;
    typedef typename SetFPtr::iterator                 SetFPtrIter;
    typedef typename Facet::VecVPtrNVConstIter         VecVPtrNVConstIter;

    const unsigned numNeighbors = Facet::numNeighbors;
    const corlib::shape myShape = Facet::myShape;

    oneRingF.clear( );
    closedRing = true;

    for ( unsigned n = 0; n < numNeighbors; n++ ) {
        Facet * neighFacet = facet->neighbor( n );
        if ( neighFacet ) {
            oneRingF.push_back( neighFacet );
            VecFPtrIter nIt = --(oneRingF.end());

            // go clock-wise
            Walker< Facet > walker( facet, n );
            walker = walker.opposite().next().opposite();
            while ( ( walker.getFacet() != NULL ) and
                    ( walker.getFacet() != facet->neighbor( ShapeProp< Facet::myShape >::previous[ n ] ) ) ) {
                VecFPtrIter fIt = std::find( oneRingF.begin(), oneRingF.end(), walker.getFacet() );
                if ( fIt == oneRingF.end() )
                    nIt = oneRingF.insert( nIt, walker.getFacet() );
                walker = walker.next().opposite();
            }
            if ( walker.getFacet() == NULL ) {
                closedRing = false;
            }

            // go counter-clockwise
            walker.set( facet, n );
            walker = walker.opposite().previous().opposite();
            while ( ( walker.getFacet() != NULL ) and
                    ( walker.getFacet() != facet->neighbor( ShapeProp< Facet::myShape >::next[ n ] ) ) ) {
                VecFPtrIter fIt = std::find( oneRingF.begin(), oneRingF.end(), walker.getFacet() );
                if ( fIt == oneRingF.end() )
                    oneRingF.push_back( walker.getFacet() );
                walker = walker.previous().opposite();
            }
            if ( walker.getFacet() == NULL ) {
                closedRing = false;
            }
        }
    }

    // collect vertices
    for ( VecFPtrIter f = oneRingF.begin(); f != oneRingF.end(); ++f ) {
        oneRingV.insert( (*f)->verticesBegin( ), (*f)->verticesEnd( ) );
    }
    for ( VecVPtrNVConstIter v = facet->verticesBegin( ); v != facet->verticesEnd( ); ++v ) {
        oneRingV.erase( *v );
    }

    return;
}

//==============================================================================
template< typename FACET >
void subdiv::surf::collectOneRingFacets( FACET * facet,
                                       std::vector< FACET * > & oneRingF,
                                       std::set< typename FACET::Vertex * > & oneRingV )
{

    typedef FACET                                      Facet;
    typedef std::vector< Facet * >                     VecFPtr;
    typedef typename VecFPtr::iterator                 VecFPtrIter;
    typedef std::set< Facet * >                        SetFPtr;
    typedef typename SetFPtr::iterator                 SetFPtrIter;
    typedef typename Facet::VecVPtrNVConstIter         VecVPtrNVConstIter;

    const unsigned numNeighbors = Facet::numNeighbors;
    const unsigned numVertices = Facet::numVertices;

    // immediate neigbours at edges
    SetFPtr ordNeighF;   ordNeighF.insert( facet );
    for ( unsigned n = 0; n < numNeighbors; n++ ) {
        VecFPtr edgeNeighF;
        facet->neighborsAtEdge( n, edgeNeighF );
        ordNeighF.insert( edgeNeighF.begin(), edgeNeighF.end() );
    }

    // find indirect neighbours
    unsigned newNeigh;
    do {
        newNeigh = 0;

        // find new neighbours
        for ( SetFPtrIter f = ordNeighF.begin(); f != ordNeighF.end(); ++f ) {
            SetFPtr fNeighF;   (*f)->allNeighbors( fNeighF );
            for ( SetFPtrIter fn = fNeighF.begin(); fn != fNeighF.end(); ++fn ) {
                // check if candidate is already in one-ring
                if ( ordNeighF.find( *fn ) == ordNeighF.end() ) {
                    // verify that at least one vertex is shared
                    for ( unsigned i = 0; i < numVertices; ++i ) {
                        if ( std::find( facet->verticesBegin( ), facet->verticesEnd( ), (*fn)->vertex( i ) )
                             != facet->verticesEnd( ) ) {
                            ordNeighF.insert( *fn );
                            newNeigh++;
                            break;
                        }
                    }
                }
            }
        }

    } while ( newNeigh > 0 );

    // copy to vector
    ordNeighF.erase( facet );
    std::copy( ordNeighF.begin(), ordNeighF.end(), std::back_inserter( oneRingF ) );

    // collect vertices
    for ( VecFPtrIter f = oneRingF.begin(); f != oneRingF.end(); ++f ) {
        oneRingV.insert( (*f)->verticesBegin( ), (*f)->verticesEnd( ) );
    }
    for ( VecVPtrNVConstIter v = facet->verticesBegin( ); v != facet->verticesEnd( ); ++v ) {
        oneRingV.erase( *v );
    }

    return;
}

//==============================================================================
template< typename FACET >
void subdiv::surf::collectNRingFacets( const unsigned order,
                                       FACET * facet,
                                       std::vector< FACET * > & nRingF,
                                       std::set< typename FACET::Vertex * > & nRingV )
{
    // nothing to collect if 0-ring is requested
    if ( order == 0 )
        return;

    // types
    typedef std::vector< FACET * >                 VecFPtr;
    typedef typename VecFPtr::iterator             VecFPtrIter;
    typedef std::set< FACET * >                    SetFPtr;
    typedef typename std::set< FACET * >::iterator SetFPtrIter;

    // collect 1-ring
    subdiv::surf::collectOneRingFacets( facet, nRingF, nRingV );

    // collect N-ring
    if ( order > 1 ) {
        // insert the facet at _front_
        nRingF.insert( nRingF.begin(), facet );

        unsigned begin = 1;  // at 1, because at 0 is #facet, whose 1-ring we already have
        unsigned end = nRingF.size();

        for ( unsigned i = 1; i < order; ++i ) {
            for ( unsigned f = begin; f < end; ++f ) {
                VecFPtr oneRingF;
                subdiv::surf::collectOneRingFacets( nRingF[ f ], oneRingF, nRingV );
                for ( VecFPtrIter s = oneRingF.begin(); s != oneRingF.end(); ++s ) {
                    if ( std::find( nRingF.begin(), nRingF.end(), *s ) == nRingF.end() ) {
                        nRingF.push_back( *s );
                        nRingV.insert( (*s)->verticesBegin(), (*s)->verticesEnd() );
                    }
                }
            }
            begin = end;
            end = nRingF.size();

        }

        // remove centre facet from _front_
        nRingF.erase( nRingF.begin() );

        // remove direct vertices of centre facet
        typedef typename FACET::VecVPtrNVConstIter         VecVPtrNVConstIter;
        for ( VecVPtrNVConstIter v = facet->verticesBegin( ); v != facet->verticesEnd( ); ++v ) {
            nRingV.erase( *v );
        }
    }

    // done
    return;
}

//------------------------------------------------------------------------------
template< typename FACET >
bool subdiv::surf::collectManifoldOneRingVertices( FACET * ft,
                                                 const int ivtx,
                                                 FACET * & (FACET::*getNeighborFun)( const int ),
                                                 std::vector< typename FACET::Vertex * > & oneRingV,
                                                 std::vector< FACET * > & oneRingF,
                                                 const bool stopOnCrease )
{
    const enum corlib::shape myShape = FACET::myShape;
    typedef typename FACET::Vertex Vertex;
    const Vertex * vc = ft->vertex( ivtx );
    oneRingV.clear( );
    oneRingF.clear( );
    bool isClosed = true;  // closed 1-ring made of 1 petal (starting assumption)

    // go counter-clockwise collect on ft's petal
    oneRingF.push_back( ft );
    int cId = ivtx;
    int pId = ShapeProp< myShape >::previous[ cId ];
    typename FACET::Vertex * vp = ft->vertex( pId );
    oneRingV.push_back( vp );
    int nId = ShapeProp< myShape >::previous[ cId ];
    FACET * fn = (ft->*getNeighborFun)( nId );
    std::vector< FACET * > eNeighbors;
    ft->neighborsAtEdge( nId, eNeighbors );
    bool isCrease = ( ft->getEdgeTag( nId ) == subdiv::surf::EDGE_CREASE );
    while ( fn and ( fn != ft ) and ( eNeighbors.size() == 1 ) and 
            not ( stopOnCrease and isCrease ) ) {
        cId = fn->localVertexIndex( vc );
        if ( fn->localEdgeIndex( vc, vp ) == cId ) {     // (+)
            pId = ShapeProp< myShape >::previous[ cId ];
            nId = ShapeProp< myShape >::previous[ cId ];
        }
        else {                                           // (-)
            pId = ShapeProp< myShape >::next[ cId ];
            nId = cId;
        }
        oneRingF.push_back( fn );
        vp = fn->vertex( pId );
        oneRingV.push_back( vp );
        eNeighbors.clear( );
        fn->neighborsAtEdge( nId, eNeighbors );
        isCrease = ( fn->getEdgeTag( nId ) == subdiv::surf::EDGE_CREASE );
        fn = (fn->*getNeighborFun)( nId );
    }

    // ring on #ft's petal is open.
    if ( fn != ft ) {  // this includes fn == NULL
        isClosed = false;  // open 1-ring made of 1 petals

        // go clockwise and continue to collect on #ft's petal
        pId = ShapeProp< myShape >::next[ ivtx ];
        vp = ft->vertex( pId );
        oneRingV.insert( oneRingV.begin( ), vp );  // maintain orientation/order of vertices
        nId = ivtx;
        fn = (ft->*getNeighborFun)( nId );
        eNeighbors.clear( );
        ft->neighborsAtEdge( nId, eNeighbors );
        isCrease = ( ft->getEdgeTag( nId ) == subdiv::surf::EDGE_CREASE );
        while ( fn and ( eNeighbors.size() == 1 ) and
                not ( stopOnCrease and isCrease ) ) {  // this must end at open edge
            cId = fn->localVertexIndex( vc );
            if ( fn->localEdgeIndex( vc, vp ) != cId ) {    // (+)
                pId = ShapeProp< myShape >::next[ cId ];
                nId = cId;
            }
            else {                                          // (-)
                pId = ShapeProp< myShape >::previous[ cId ];
                nId = ShapeProp< myShape >::previous[ cId ];
            }
            oneRingF.insert( oneRingF.begin( ), fn );
            vp = fn->vertex( pId );
            oneRingV.insert( oneRingV.begin( ), vp );  // maintain orientation of vertices
            eNeighbors.clear( );
            fn->neighborsAtEdge( nId, eNeighbors );
            isCrease = ( fn->getEdgeTag( nId ) == subdiv::surf::EDGE_CREASE );
            fn = (fn->*getNeighborFun)( nId );
        }

        // set pointers to end of 1-rings on #ft's petal
        vp = *(--oneRingV.end( ));
        fn = *(--oneRingF.end( ));

        // go counter-clockwise on different petal of "same" manifold
        cId = fn->localVertexIndex( vc );
        pId = fn->localVertexIndex( vp );
        nId = fn->localEdgeIndex( vc, vp );
        eNeighbors.clear( );
        isCrease = ( fn->getEdgeTag( nId ) == subdiv::surf::EDGE_CREASE );
        fn->neighborsAtEdge( nId, eNeighbors );
        fn = (fn->*getNeighborFun)( nId );
        bool isOpenPetal = true;
        while ( fn and isOpenPetal and not ( stopOnCrease and isCrease ) ) {
            cId = fn->localVertexIndex( vc );
            if ( fn->localEdgeIndex( vc, vp ) == cId ) {     // (+)
                pId = ShapeProp< myShape >::previous[ cId ];
                nId = ShapeProp< myShape >::previous[ cId ];
            }
            else {                                           // (-)
                pId = ShapeProp< myShape >::next[ cId ];
                nId = cId;
            }
            oneRingF.push_back( fn );
            vp = fn->vertex( pId );
            eNeighbors.clear( );
            fn->neighborsAtEdge( nId, eNeighbors );
            if ( eNeighbors.size() > 1 )
                isOpenPetal = false;
            else
                oneRingV.push_back( vp );
            isCrease = ( fn->getEdgeTag( nId ) == subdiv::surf::EDGE_CREASE );
            fn = (fn->*getNeighborFun)( nId );
        }

        if ( not isOpenPetal ) {
            FTL_VERIFY( corlib::contains( eNeighbors.begin(), eNeighbors.end(),
                                          *(oneRingF.begin( )) ) );
            isClosed = true;  // closed 1-ring made of 2 petals
        }
        else {  // equivalent if ( not fn )

            isClosed = false;  // definitely open 1-ring (made of 3 petals)

            // set pointers to begin of 1-rings on #ft's petal
            vp = *(oneRingV.begin( ));
            fn = *(oneRingF.begin( ));

            // go clockwise on attached petal of "same" manifold
            cId = fn->localVertexIndex( vc );
            pId = fn->localVertexIndex( vp );
            nId = fn->localEdgeIndex( vc, vp );
            eNeighbors.clear( );
            fn->neighborsAtEdge( nId, eNeighbors );
            isCrease = ( fn->getEdgeTag( nId ) == subdiv::surf::EDGE_CREASE );
            fn = (fn->*getNeighborFun)( nId );
            bool isOpenPetal = true;
            while ( fn and isOpenPetal and not ( stopOnCrease and isCrease ) ) {
                cId = fn->localVertexIndex( vc );
                if ( fn->localEdgeIndex( vc, vp ) != cId ) {    // (+)
                    pId = ShapeProp< myShape >::next[ cId ];
                    nId = cId;
                }
                else {                                          // (-)
                    pId = ShapeProp< myShape >::previous[ cId ];
                    nId = ShapeProp< myShape >::previous[ cId ];
                }
                oneRingF.insert( oneRingF.begin( ), fn );
                vp = fn->vertex( pId );
                eNeighbors.clear( );
                fn->neighborsAtEdge( nId, eNeighbors );
                if ( eNeighbors.size() > 1 )
                    isOpenPetal = false;
                else
                    oneRingV.insert( oneRingV.begin( ), vp );  // maintain orientation of vertices
                isCrease = ( fn->getEdgeTag( nId ) == subdiv::surf::EDGE_CREASE );
                fn = (fn->*getNeighborFun)( nId );
            }
            FTL_VERIFY( isOpenPetal );
        }
    }

    if ( isClosed )
        FTL_VERIFY( oneRingV.size( ) == oneRingF.size( ) );

    return isClosed;
}

//==============================================================================
template< typename FACET >
void subdiv::surf::collectOneRingVertices(
    FACET * ft,
    const int ivtx,
    std::map< typename FACET::Vertex *,
              std::pair< enum subdiv::surf::edgeTag,
                         enum subdiv::surf::edgeCategory > > & oneRingV,
    std::vector< FACET * > & oneRingF
    )
{
    // a few types
    typedef FACET                                           Facet;
    typedef typename Facet::Vertex                          Vertex;
    typedef std::tuple< Facet *, int, int >                 SomeTuple;
    typedef std::vector< SomeTuple >                        SomeVec;
    typedef typename SomeVec::iterator                      SomeVecIter;
    typedef enum subdiv::surf::edgeTag                      EdgeTag;
    typedef enum subdiv::surf::edgeCategory                 EdgeCategory;

    const corlib::shape myShape = Facet::myShape;

    // vertex to pointer whose one-ring is collected
    Vertex * vc = ft->vertex( ivtx );

    // go counter-clockwise
    SomeVec eNEV;
    eNEV.push_back( std::make_tuple( ft,
                                       ShapeProp< myShape >::previous[ ivtx ],
                                       ShapeProp< myShape >::previous[ ivtx ] ) );
    
    unsigned start = 0;
    unsigned end = eNEV.size();  // equals 1
    while ( start < end ) {
        // only loop over newly collected facets
        for ( unsigned j = start; j < end; ++j ) {
            Facet * fn = std::get< 0 >( eNEV[ j ] );
            int nEdge = std::get< 1 >( eNEV[ j ] );
            Vertex * vp = fn->vertex( std::get< 2 >( eNEV[ j ] ) );
            assert( vp != vc );

            // collect neighbours at edge
            std::vector< Facet * > eNeighbors;
            fn->neighborsAtEdge( nEdge, eNeighbors );
            
            // loop neighbors at edge
            for ( unsigned n = 0; n < eNeighbors.size(); ++n ) {
                if ( eNeighbors[ n ] != ft ) {
                    const int nVertex = eNeighbors[ n ]->localVertexIndex( vc );
                    const int nEdge = eNeighbors[ n ]->localEdgeIndex( vc, vp );
                    SomeTuple newNEV;
                    if ( nVertex == nEdge ) {
                        newNEV = std::make_tuple( eNeighbors[ n ],
                                                    ShapeProp< myShape >::previous[ nEdge ],
                                                    ShapeProp< myShape >::previous[ nVertex ] ); // (+)
                    }
                    else {
                        newNEV = std::make_tuple( eNeighbors[ n ],
                                                    ShapeProp< myShape >::next[ nEdge ],
                                                    ShapeProp< myShape >::next[ nVertex ] );  // (-)
                    }
                    if ( std::find( eNEV.begin(), eNEV.end(), newNEV ) == eNEV.end() ) {
                        eNEV.push_back( newNEV );
                    }
                }
            }
        }
        start = end;
        end = eNEV.size();
    }

    // go clockwise
    eNEV.push_back( std::make_tuple( ft, ivtx, ShapeProp< myShape >::next[ ivtx ] ) );

    start = end;
    end = eNEV.size();
    while ( start < end ) {
        // only loop over newly collected facets
        for ( unsigned j = start; j < end; ++j ) {
            Facet * fn = std::get< 0 >( eNEV[ j ] );
            int nEdge = std::get< 1 >( eNEV[ j ] );
            Vertex * vp = fn->vertex( std::get< 2 >( eNEV[ j ] ) );
            assert( vp != vc );

            // collect neighbours at edge
            std::vector< Facet * > eNeighbors;
            fn->neighborsAtEdge( nEdge, eNeighbors );

            // loop neighbors at edge
            for ( unsigned n = 0; n < eNeighbors.size(); ++n ) {
                if ( eNeighbors[ n ] != ft ) {
                    const int nVertex = eNeighbors[ n ]->localVertexIndex( vc );
                    const int nEdge = eNeighbors[ n ]->localEdgeIndex( vc, vp );
                    SomeTuple newNEV;
                    if ( nVertex != nEdge ) {
                        newNEV = std::make_tuple( eNeighbors[ n ],
                                                    ShapeProp< myShape >::next[ nEdge ],
                                                    ShapeProp< myShape >::next[ nVertex ] ); // (+)
                    }
                    else {
                        newNEV = std::make_tuple( eNeighbors[ n ],
                                                    ShapeProp< myShape >::previous[ nEdge ],
                                                    ShapeProp< myShape >::previous[ nVertex ] );  // (-)
                    }
                    if ( std::find( eNEV.begin(), eNEV.end(), newNEV ) == eNEV.end() ) {
                        eNEV.push_back( newNEV );
                    }
                }
            }
        }
        start = end;
        end = eNEV.size();
    }

    // loop collected neighboring facets at vertex
    for ( SomeVecIter n = eNEV.begin(); n != eNEV.end(); ++n ) {
        Facet * fn = std::get< 0 >( *n );
        int nEdge = std::get< 1 >( *n );
        Vertex * vp = fn->vertex( std::get< 2 >( *n ) );

        // get edge tag
        EdgeTag eTag = fn->getEdgeTag( nEdge );

        // collect neighbours at edge
        std::vector< Facet * > eNeighbors;
        fn->neighborsAtEdge( nEdge, eNeighbors );

        // deterine edge category
        EdgeCategory eCat = subdiv::surf::ECAT_UNKNOWN;
        if ( eNeighbors.size() == 0 ) {
            eCat = subdiv::surf::ECAT_BOUNDARY;
        }
        else if ( eNeighbors.size() == 1 ) {
            eCat = subdiv::surf::ECAT_MANIFOLD;
        }
        else {
            eCat = subdiv::surf::ECAT_NONMANIFOLD;
        }
        
        oneRingV.insert( std::make_pair( vp, std::make_pair( eTag, eCat ) ) );

        // store facets
        if ( not corlib::contains( oneRingF.begin(), oneRingF.end(), fn ) )
             oneRingF.push_back( fn );
    }

    // done
    return;
}

//------------------------------------------------------------------------------
template< typename VERTEX >
void subdiv::surf::extractByEdgeTag(
    const std::map< VERTEX *,
                    std::pair< enum subdiv::surf::edgeTag,
                               enum subdiv::surf::edgeCategory > > & oneRingV,
    const enum subdiv::surf::edgeTag eTag,
    std::vector< VERTEX * > & vertices
    )
{
    typedef std::pair< enum subdiv::surf::edgeTag,
                       enum subdiv::surf::edgeCategory > Pair;
    typedef std::map< VERTEX *, Pair > Map;

    corlib::transform_if( oneRingV.begin( ), oneRingV.end( ),
                          std::back_inserter( vertices ),
                          std::bind( &Map::value_type::first, std::placeholders::_1 ),
                          std::bind( std::equal_to< enum subdiv::surf::edgeTag >(),
                                       std::bind( &Pair::pair::first,
                                                    std::bind( &Map::value_type::second, std::placeholders::_1 ) ),
                                       eTag ) );

    return;
}

//------------------------------------------------------------------------------
template< typename VERTEX >
void subdiv::surf::extractByEdgeCategory(
    const std::map< VERTEX *,
                    std::pair< enum subdiv::surf::edgeTag,
                               enum subdiv::surf::edgeCategory > > & oneRingV,
    const enum subdiv::surf::edgeCategory eCat,
    std::vector< VERTEX * > & vertices
    )
{
    typedef std::pair< enum subdiv::surf::edgeTag,
                       enum subdiv::surf::edgeCategory > Pair;
    typedef std::map< VERTEX *, Pair > Map;

    corlib::transform_if( oneRingV.begin( ), oneRingV.end( ),
                          std::back_inserter( vertices ),
                          std::bind( &Map::value_type::first, std::placeholders::_1 ),
                          std::bind( std::equal_to< enum subdiv::surf::edgeCategory >(),
                                       std::bind( &Pair::pair::second,
                                                    std::bind( &Map::value_type::second, std::placeholders::_1 ) ),
                                       eCat ) );

    return;
}

//==============================================================================

//------------------------------------------------------------------------------
template< typename FACET >
void subdiv::surf::collectManifoldRegularVertices(
    FACET * centreFacet,
    std::array< typename FACET::Vertex *,
                           Subdivision< subdiv::surf::CATMULL_CLARK,
                                        FACET >::numFunctions > & regularVertices
    )
{
    typedef FACET FacetTree;
    const unsigned numVertices = FacetTree::numVertices;
    const subdiv::surf::method myMethod = subdiv::surf::CATMULL_CLARK;
    const unsigned numFunctions = Subdivision< myMethod, FACET >::numFunctions;

    // index translation of collection route
    const std::array< unsigned, numFunctions > route = { {
            5, 6, 10, 9, 0, 1, 2, 3, 7, 11, 15, 14, 13, 12, 8, 4
        } };
    // initialise walker
    Walker< FacetTree > walker( centreFacet, 0 );   
    // inner ring
    for ( unsigned i = 0; i < numVertices; ++i ) {
        regularVertices[ route[ i ] ] = walker.getVertex( );
        walker.next( );
    }
    // outer ring
    walker = walker.opposite( ).next( ).opposite( ).previous( );
    unsigned i = numVertices;
    for ( unsigned l = 0; l < numVertices; ++l ) {
        regularVertices[ route[ i++ ] ] = walker.getVertex( );
        walker.next( );
        regularVertices[ route[ i++ ] ] = walker.getVertex( );
        walker = walker.opposite( ).next( ).next( );
        regularVertices[ route[ i++ ] ] = walker.getVertex( );
        walker = walker.opposite( ).next( ).next( );
    }

    return;
}

//------------------------------------------------------------------------------
template< typename FACET >
void subdiv::surf::collectManifoldRegularVertices(
    FACET * centreFacet,
    std::array< typename FACET::Vertex *,
                           Subdivision< subdiv::surf::LOOP,
                                        FACET >::numFunctions > & regularVertices
    )
{
    typedef FACET FacetTree;
    const unsigned numVertices = FacetTree::numVertices;
    const subdiv::surf::method myMethod = subdiv::surf::LOOP;
    const unsigned numFunctions = Subdivision< myMethod, FACET >::numFunctions;

    // index translation of collection route
    const std::array< unsigned, numFunctions > route = { {
            3, 6, 7,  0, 2, 5,  9, 10, 11,  8, 4, 1
        } };
    // initialise walker
    Walker< FacetTree > walker( centreFacet, 0 );
    // inner ring
    for ( unsigned i = 0; i < numVertices; ++i ) {
        regularVertices[ route[ i ] ] = walker.getVertex( );
        walker.next( );
    }
    // outer ring
    walker = walker.opposite( ).next( ).opposite( ).next( ).next( );
    unsigned i = numVertices;
    for ( unsigned l = 0; l < numVertices; ++l ) {
        regularVertices[ route[ i++ ] ] = walker.getVertex( );
        walker.next( );
        regularVertices[ route[ i++ ] ] = walker.getVertex( );
        walker = walker.opposite( ).next( ).opposite( ).next( ).next( );
        regularVertices[ route[ i++ ] ] = walker.getVertex( );
        walker = walker.opposite( ).next( ).next( ).opposite( ).next( );
    }

    return;
}

#endif
        
