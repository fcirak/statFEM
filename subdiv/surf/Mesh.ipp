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

#include <cassert>
#include <cstdlib>
#include <algorithm>
#include <utility>
#include <functional>
#include <iterator>

#include <corlib/SmfHead.hpp>
#include <corlib/verify.hpp>

#include <subdiv/surf/MeshTags.hpp>
#include <subdiv/surf/ShapeDim.hpp>

//------------------------------------------------------------------------------
/// Local index of next vertex
template< typename V, typename E, typename F, typename FT >
const std::array< int, corlib::ShapeTraits< FT::myShape >::numVertices >
subdiv::surf::Mesh< V, E, F, FT >::next_ = 
    subdiv::surf::ShapeProp< FT::myShape >::next;

/// Local index of previous vertex
template< typename V, typename E, typename F, typename FT > 
const std::array< int, corlib::ShapeTraits< FT::myShape >::numVertices >
subdiv::surf::Mesh< V, E, F, FT >::previous_ = 
    subdiv::surf::ShapeProp< FT::myShape >::previous;

//------------------------------------------------------------------------------
template< typename V, typename E, typename F, typename FT > 
subdiv::surf::Mesh< V, E, F, FT >::Mesh( std::istream & smf ) :
    verbose_( false )
{    
    static_assert( Facet::myShape == FTree::myShape );

    this->readSmf( smf );

    return;
}

//------------------------------------------------------------------------------
/// Destructor
template< typename V, typename E, typename F, typename FT >
subdiv::surf::Mesh< V, E, F, FT >::~Mesh( )
{    
    // remove facet tree
    this->iterateOverFacetTrees( corlib::deleteFunctor( ) );
    facetTrees_.clear( );

    // remove facets (level 0)
    this->iterateOverFacets( corlib::deleteFunctor( ) );
    facets_.clear( );

    // remove nodes
    this->iterateOverVertices( corlib::deleteFunctor( ) );
    vertices_.clear( );

    // remove edges
    std::for_each( edges_.begin( ), edges_.end( ), corlib::deleteFunctor( ) );
    edges_.clear( );
}

//------------------------------------------------------------------------------
template< typename V, typename E, typename F, typename FT >
template< typename OP >
OP subdiv::surf::Mesh< V, E, F, FT >::iterateOverVertices( OP op )
{
    VertexIterator begin = vertices_.begin( );
    VertexIterator end = vertices_.end( );
    return std::for_each( begin, end, op );
}

//------------------------------------------------------------------------------
template< typename V, typename E, typename F, typename FT>
template< typename OP >
OP subdiv::surf::Mesh< V, E, F, FT >::iterateOverEdges( OP op )
{
    const EdgeIterator begin = edges_.begin( );
    const EdgeIterator end = edges_.end( );
    return std::for_each( begin, end, op );
}

//------------------------------------------------------------------------------
template< typename V, typename E, typename F, typename FT>
template< typename OP >
OP subdiv::surf::Mesh< V, E, F, FT >::iterateOverFacets( OP op )
{
    FacetIterator begin = facets_.begin( );
    FacetIterator end = facets_.end( );
    return std::for_each( begin, end, op );
}

//------------------------------------------------------------------------------
template< typename V, typename E, typename F, typename FT>
template< typename OP >
OP subdiv::surf::Mesh< V, E, F, FT >::iterateOverFacetTrees( OP op )
{
    FacetIterator begin = this->fBegin();
    FacetIterator end   = this->fEnd();
    for ( ; begin!=end; ++begin )
        if ( facetTrees_[ *begin ] )
            op( facetTrees_[ *begin ] );
    return op;
}

//------------------------------------------------------------------------------
/// Read SMF file
template< typename V, typename E, typename F, typename FT > 
void subdiv::surf::Mesh< V, E, F, FT >::readSmf( std::istream & smf )
{
    // validate input file
    corlib::SmfHead smfHead;
    smfHead.readValidated( smf, myShape, numVertices );

    // number of vertices and faces
    unsigned numVectors = 0;
    unsigned numFaces = 0;
    smf >> numVectors >> numFaces;
    
    // reserve space
    vertices_.reserve( numVectors );
    facets_.reserve( numFaces );

    // add vertices
    for ( unsigned i = 0; i < numVectors; ++i ) {
        Vertex * v = new Vertex( );
        smf >> v;
        this->addVertex_( v );
        v->setIndex( i );   //TEMP
    }

    // add elements
    for ( unsigned i = 0; i < numFaces; ++i ) {
        VecVPtrNV nodesOfFace;
        // collect vertices by node IDs
        for ( unsigned j = 0; j < numVertices; ++j ) {
            unsigned id;
            smf >> id;
            nodesOfFace[ j ] = vertices_[ id ];  // set pointer to node
        }

        Facet * f = new Facet( nodesOfFace );
        this->addFacet_( f );
        f->setIndex( i );   //TEMP
    }

    return;
}

//------------------------------------------------------------------------------
/// Read tag file
template< typename V, typename E, typename F, typename FT > 
void subdiv::surf::Mesh< V, E, F, FT >::readTags( std::istream & tg )
{
    typedef typename Vertex::VertexTag     VertexTag;
    typedef typename Edge::EdgeTag         EdgeTag;

    std::string name;    
    tg >> name;
    FTL_VERIFY( name == "TG" );

    unsigned numVertexTags, numEdgeTags;
    tg >> numVertexTags >> numEdgeTags;

    // set node tags
    for ( unsigned vt = 0; vt < numVertexTags; vt ++ ) {
        unsigned v;
        int tagInt;
        tg >> v >> tagInt;

        const VertexTag tag = static_cast< VertexTag >( tagInt );
        vertices_[ v ]->setTag( tag );
    }

    // store edge tags in edges_ container
    for ( unsigned et = 0; et < numEdgeTags; et ++ ) {
        unsigned v0;
        unsigned v1;
        int tagInt;
        tg >> v0 >> v1 >> tagInt;

        // create edge based on start and end node
        Edge * e = new Edge( vertices_[ v0 ], vertices_[ v1 ] );
        // insert into edge set
        const std::pair< EdgeSetConstIterator, bool > ei = edges_.insert( e );
        // if edge already appeared in the EdgeSet, then ei's bool is false
        if ( not ei.second )  delete e;
        // set tag
        const EdgeTag tag = static_cast< EdgeTag >( tagInt );
        // check if additional information must be read
        if ( tag == subdiv::surf::EDGE_SMOOTHCREASE ) {
            unsigned f0, f1;
            tg >> f0 >> f1;
            (*ei.first)->setTag( tag );
            (*ei.first)->insertFacet( facets_[ f0 ], facets_[ f1 ] );
        }
        else {
            if ( (*ei.first)->getTag( ) == subdiv::surf::EDGE_NOTAG )
                (*ei.first)->setTag( tag );
        }
    }

    return;
}

//------------------------------------------------------------------------------
/// Build the tree structure
template< typename V, typename E, typename F, typename FT > 
void subdiv::surf::Mesh< V, E, F, FT >::buildFacetTopology( ) 
{
    // inform curious user
    if ( verbose_ ) {
        std::cout << "Mesh contains "
                  << vertices_.size( ) << " vertices and "
                  << facets_.size() << " facets"
                  << std::endl;
    }

    // facet range
    const FacetIterator fBegin = this->fBegin();
    const FacetIterator fEnd = this->fEnd();

    // loop facets and create edge as well establish neighbouring relations
    for ( FacetIterator fIter = fBegin; fIter != fEnd; ++fIter ) {
        Facet * f = *fIter;
        for ( unsigned i = 0; i < numVertices; i++ ) {
            this->addOrMatch_( f->vertex( i ), f->vertex( next_[ i ] ), f );
        }
    }

    // inform nosey user
    if ( verbose_ ) {
        unsigned numEdgeManifold = 0;
        unsigned numEdgeBoundary = 0;
        unsigned numEdgeNonManifold = 0;

        for ( EdgeSetConstIterator e = edges_.begin(); e != edges_.end( ); ++e ) {
            if ( (*e)->getCategory( ) == subdiv::surf::ECAT_MANIFOLD )
                numEdgeManifold += 1;
            else if ( (*e)->getCategory( ) == subdiv::surf::ECAT_BOUNDARY )
                numEdgeBoundary += 1;
            else if ( (*e)->getCategory( ) == subdiv::surf::ECAT_NONMANIFOLD )
                numEdgeNonManifold += 1;
        }
        assert( numEdgeManifold+numEdgeBoundary+numEdgeNonManifold == edges_.size() );

        std::cout << "Mesh contains "
                  << numEdgeManifold << " manifold, "
                  << numEdgeBoundary << " boundary, and "
                  << numEdgeNonManifold << " non-manifold edges"
                  << std::endl;
    }
    
    // make sure every manifold is consistently oriented
    // set-up of (true) manifold neigbourships
    this->orientateManifolds_( );

    // set edge tags
    for ( EdgeIterator eIter = edges_.begin(); eIter != edges_.end(); ++eIter ) {
        // convenience
        typedef typename Edge::EdgeTag                  EdgeTag;
        typedef typename Edge::FacetPair                FacetPair;
        typedef typename Edge::FacetPairIter            FacetPairIter;

        // crease boundary and non-manifold edges
        if ( (*eIter)->getCategory( ) == subdiv::surf::ECAT_BOUNDARY or
             (*eIter)->getCategory( ) == subdiv::surf::ECAT_NONMANIFOLD ) {
            if ( (*eIter)->getTag( ) == subdiv::surf::EDGE_NOTAG ) {
                (*eIter)->setTag( subdiv::surf::EDGE_CREASE );
            }
        }

        // set edge tags in facets
        for ( FacetPairIter fIter = (*eIter)->facetsBegin( );
              fIter != (*eIter)->facetsEnd( ); ++fIter ) {
            const EdgeTag eTag = (*eIter)->getTag( );
            Facet * f = fIter->first;
            const int eIndex = f->localEdgeIndex( (*eIter)->first( ),
                                                  (*eIter)->second( ) );
            // let the master facet store the tag
            if ( f == (*eIter)->getFacet( ) )
                f->setEdgeTag( eIndex, eTag );
            // associate facet IDs for "smooth creases"
            if ( eTag == subdiv::surf::EDGE_SMOOTHCREASE ) {
                Facet * fn = fIter->second;
                if ( fn ) f->manifoldNeighbor( eIndex ) = fn;
            }
        }

        // propagate crease information to vertices
        if ( (*eIter)->getTag( ) == subdiv::surf::EDGE_CREASE ) {
            Vertex * v0 = (*eIter)->first( );
            if ( v0->tag( ) == subdiv::surf::VERTEX_NOTAG ) {
                v0->setTag( subdiv::surf::VERTEX_CREASE );
            }
            Vertex * v1 = (*eIter)->second( );
            if ( v1->tag( ) == subdiv::surf::VERTEX_NOTAG ) {
                v1->setTag( subdiv::surf::VERTEX_CREASE );
            }
        }
        if ( (*eIter)->getTag( ) == subdiv::surf::EDGE_SMOOTHCREASE ) {
            Vertex * v0 = (*eIter)->first( );
            if ( v0->tag( ) == subdiv::surf::VERTEX_NOTAG ) {
                v0->setTag( subdiv::surf::VERTEX_SMOOTHCREASE );
            }
            Vertex * v1 = (*eIter)->second( );
            if ( v1->tag( ) == subdiv::surf::VERTEX_NOTAG ) {
                v1->setTag( subdiv::surf::VERTEX_SMOOTHCREASE );
            }
        }
    }

    // set corner/dart vertex tags, if vertex has _not_ exactly 2 creased edges
    // note that non-manifold and boundary edges are creases as well
#if 1
    const VertexIterator vBegin = this->vBegin( );
    const VertexIterator vEnd = this->vEnd( );
    for ( VertexIterator vIter = vBegin; vIter != vEnd; ++vIter ) {
        if ( (*vIter)->tag( ) == subdiv::surf::VERTEX_CREASE or
             (*vIter)->tag( ) == subdiv::surf::VERTEX_SMOOTHCREASE ) {
            Vertex * vPtr = *vIter;
            // find an edge which contains the vertex
            EdgeIterator eIter = std::find_if( 
                edges_.begin( ), edges_.end( ),
                std::bind( std::logical_or< bool >(),
                             std::bind( std::equal_to< Vertex * >(),
                                          std::bind( &Edge::first, std::placeholders::_1 ), vPtr ),
                             std::bind( std::equal_to< Vertex * >(),
                                          std::bind( &Edge::second, std::placeholders::_1 ), vPtr ) )
                );
            assert( eIter != edges_.end() );
            Facet * fPtr = (*eIter)->getFacet( );
            const unsigned vIndex = fPtr->localVertexIndex( vPtr );
            std::map< Vertex *, std::pair< enum subdiv::surf::edgeTag,
                                           enum subdiv::surf::edgeCategory > > oneRingV;
            std::vector< Facet * > oneRingF;
            subdiv::surf::collectOneRingVertices( fPtr, vIndex, oneRingV, oneRingF );
            std::vector< Vertex * > creaseV;
            subdiv::surf::extractByEdgeTag( oneRingV, subdiv::surf::EDGE_CREASE, creaseV );
            const unsigned numCreases = creaseV.size( );
            subdiv::surf::extractByEdgeTag( oneRingV, subdiv::surf::EDGE_SMOOTHCREASE, creaseV );
            const unsigned numSmoothCreases = creaseV.size( ) - numCreases;
            const unsigned numAllCreases = creaseV.size();
            std::vector< Vertex * > nonManifV;
            subdiv::surf::extractByEdgeCategory( oneRingV, subdiv::surf::ECAT_NONMANIFOLD, nonManifV );
            if ( numAllCreases == 1 ) {
                vPtr->setTag( subdiv::surf::VERTEX_DART );
            }
            else if ( numAllCreases == 2 ) {
                if ( numSmoothCreases == 1 and numCreases == 1 )
                    vPtr->setTag( subdiv::surf::VERTEX_DART );
                if ( nonManifV.size( ) == 1 )
                    vPtr->setTag( subdiv::surf::VERTEX_DART );
            }
            else if ( numAllCreases > 2 ) {
                if ( numSmoothCreases > 0 )
                    vPtr->setTag( subdiv::surf::VERTEX_SMOOTHCORNER );
                else 
                    vPtr->setTag( subdiv::surf::VERTEX_CORNER );
            }
        }
    }
#endif

    return;
}

//------------------------------------------------------------------------------
/// Build the tree structure
template< typename V, typename E, typename F, typename FT > 
void subdiv::surf::Mesh< V, E, F, FT >::implantFacetTrees( )
{
    const FacetIterator fBegin = this->fBegin();
    const FacetIterator fEnd = this->fEnd();

    // instantiate Facet Trees, i.e. allocate an empty tree for each
    // top level facet
    for ( FacetIterator fIter = fBegin; fIter != fEnd; ++fIter ) {
        // collect facet vertices
        VecVPtrNV vec;
        std::copy( (*fIter)->verticesBegin(), (*fIter)->verticesEnd(), vec.begin() );

        // create adjoint facet tree object
        FTree * ft = new FTree( vec, 0, NULL );
        FTL_VERIFY( ft != NULL );
        ft->setIndex( (*fIter)->index( ) );  //TEMP
        std::pair< FacetTreeIterator, bool > result = 
            facetTrees_.insert( std::make_pair( *fIter, ft ) );
        FTL_VERIFY( result.second );
    }
    
    // copy neighbor information from facets to facet trees
    // (we have to do this _after_ all facet trees were allocated)
    for ( FacetIterator fIter = fBegin; fIter != fEnd; ++fIter ) {
        FTree * ft = facetTrees_[ *fIter ];

        // mimick non-manifold neighbours
        for ( unsigned i = 0; i < numNeighbors; i++ )
            ft->neighbor( i ) = facetTrees_[ (*fIter)->neighbor( i ) ];

        // mimick manifold neighbours
        for ( unsigned i = 0; i < numNeighbors; i++ )
            ft->manifoldNeighbor( i ) = facetTrees_[ (*fIter)->manifoldNeighbor( i ) ];

        // copy edge tags
        for ( unsigned e = 0; e < numEdges; ++e ) {
            const bool localEdgeTag = true;
            ft->setEdgeTag( e, (*fIter)->getEdgeTag( e, localEdgeTag ) );
        }
    }

    return;
}

//------------------------------------------------------------------------------
/// Subdivide mesh
template< typename V, typename E, typename F, typename FT >  
template< typename SUBDIV >
void subdiv::surf::Mesh< V, E, F, FT >::subdivide( SUBDIV * subdiv ) 
{
    // subdivide all elements
    for ( FacetIterator it = this->fBegin(); it != this->fEnd(); ++it ) {
        facetTrees_[ *it ]->subdivide( subdiv );
    }
    return;
}

//------------------------------------------------------------------------------
/// Write mesh to output file
template< typename V, typename E, typename F, typename FT > 
void subdiv::surf::Mesh< V, E, F, FT >::writeMeshWithIndex( std::ostream & smf,
                                                            std::ostream * tg,
                                                            const NamedAddCoord addCoord,
                                                            std::ostream * vdat )
{
    typedef typename std::vector< Vertex * >        VecVPtr;
    typedef typename VecVPtr::iterator              VecVPtrIter;
    typedef typename std::vector< FTree * >         VecFTPtr;
    typedef typename VecFTPtr::iterator             VecFTPtrIter;

    // facet range
    const FacetIterator fBegin = facets_.begin();
    const FacetIterator fEnd   = facets_.end();

    // index vertices and count vertices and facets
    unsigned numVerticesMesh, numLeafs;
    this->updateIndices_( numVerticesMesh, numLeafs );

    // write header
    corlib::SmfHead smfHead;
    smfHead.setElementShape( myShape );
    smfHead.setElementNumPoints( numVertices );
    smfHead.write( smf );
    smf << numVerticesMesh << " " << numLeafs << std::endl;

    // write control vertex coordinates
    const VertexIterator vBegin = vertices_.begin();
    const VertexIterator vEnd = vertices_.end();
    for ( VertexIterator vIt = vBegin; vIt != vEnd; ++vIt ) {
        smf << *vIt << std::endl;
    }
    // write new vertex coordinates
    for ( FacetIterator f = fBegin; f != fEnd; ++f ) {
        VecVPtr newVerts;
        facetTrees_[ *f ]->addNewVertices( std::back_inserter( newVerts ) );
        for ( VecVPtrIter vIt = newVerts.begin( ); vIt != newVerts.end( ); ++vIt ) {
            smf << *vIt << std::endl;
        }
    }

    // write vertex data on additional vertex data
    if ( not addCoord.empty() and vdat ) {
        // write header
        *vdat <<  numVerticesMesh << " " << addCoord.size() << std::endl;
        unsigned numTotalDataPerVertex = 0;
        for ( NamedAddCoord::const_iterator ac = addCoord.begin(); ac != addCoord.end(); ++ac ) {
            *vdat << ac->first << " " << ac->second << std::endl;
            numTotalDataPerVertex += ac->second;
        }
        FTL_VERIFY( numTotalDataPerVertex <= Vertex::getNumAddCoords() );

        // write coord
        unsigned offset = Vertex::dim;
        for ( NamedAddCoord::const_iterator ac = addCoord.begin(); ac != addCoord.end(); ++ac ) {
            const unsigned numDataPerVertex = ac->second;
            // control vertices
            for ( VertexIterator vIt = vBegin; vIt != vEnd; ++vIt ) {
                const typename Vertex::Point & coord = (*vIt)->levBegin()->second;
                for ( unsigned a = 0; a < numDataPerVertex; ++a ) {
                    *vdat << coord[ offset + a ] << " ";
                }
                *vdat << std::endl;
            }
            // new vertices
            for ( FacetIterator f = fBegin; f != fEnd; ++f ) {
                VecVPtr newVerts;
                facetTrees_[ *f ]->addNewVertices( std::back_inserter( newVerts ) );
                for ( VecVPtrIter vIt = newVerts.begin( ); vIt != newVerts.end( ); ++vIt ) {
                    const typename Vertex::Point & coord = (*vIt)->levBegin()->second;
                    for ( unsigned a = 0; a < numDataPerVertex; ++a ) {
                        *vdat << coord[ offset + a ] << " ";
                    }
                    *vdat << std::endl;
                }
            }
            // update
            offset += numDataPerVertex;
        }
    }

    // write leafs
    for ( FacetIterator f = fBegin; f != fEnd; ++f ) {
        VecFTPtr leafs;
        facetTrees_[ *f ]->addLeafs( std::back_inserter( leafs ) );
        for ( VecFTPtrIter lIter = leafs.begin(); lIter != leafs.end(); ++lIter ) {
            for ( unsigned i = 0; i < numVertices; i++ ) {
                smf << (*lIter)->vertex( i )->index( ) << " ";
            }
            smf << std::endl;
        }
    }

    // tag file
    if ( tg ) {
        // count number of tagged vertices and collect tagged edges
        unsigned numTaggedVertices = 0;
        for ( VertexIterator vIt = vBegin; vIt != vEnd; ++vIt ) {
            if ( (*vIt)->tag( ) != subdiv::surf::VERTEX_NOTAG ) {
                numTaggedVertices += 1;
            }
        }
        for ( FacetIterator f = fBegin; f != fEnd; ++f ) {
            VecVPtr newVerts;
            facetTrees_[ *f ]->addNewVertices( std::back_inserter( newVerts ) );
            for ( VecVPtrIter vIt = newVerts.begin( ); vIt != newVerts.end( ); ++vIt ) {
                if ( (*vIt)->tag( ) != subdiv::surf::VERTEX_NOTAG ) {
                    numTaggedVertices += 1;
                }
            }
        }

        // collect/create only tagged edges of leafs
        typedef std::multimap< Edge *, std::vector< int >, EdgeLess > EdgeMap;
        typedef typename Edge::EdgeTag                                EdgeTag;
        typedef typename EdgeMap::iterator                            EdgeMapIter;

        EdgeMap taggedEdges;
        for ( FacetIterator f = fBegin; f != fEnd; ++f ) {
            VecFTPtr leafs;
            facetTrees_[ *f ]->addLeafs( std::back_inserter( leafs ) );
            for ( VecFTPtrIter lIter = leafs.begin(); lIter != leafs.end(); ++lIter ) {
                for ( unsigned e = 0; e < numEdges; ++e ) {
                    EdgeTag eTag = (*lIter)->getEdgeTag( e );
                    if ( eTag != subdiv::surf::EDGE_NOTAG ) {
                        Edge * edge = new Edge( (*lIter)->vertex( e ),
                                                (*lIter)->vertex( next_[ e ] ) );
                        edge->setTag( eTag );
                        std::vector< int > extraInts;
                        if ( eTag == subdiv::surf::EDGE_SMOOTHCREASE ) {
                            const unsigned facetId0 = (*lIter)->index( );
                            extraInts.push_back( facetId0 );
                            FTree * fn = (*lIter)->manifoldNeighbor( e );
                            if ( fn ) {
                                const unsigned facetId1 = fn->index( );
                                extraInts.push_back( facetId1 );
                            }
                        }
                        taggedEdges.insert( std::make_pair( edge, extraInts ) );
                    }
                }
            }
        }

        // write file header
        *tg << "TG" << std::endl
            << numTaggedVertices << " " << taggedEdges.size( ) << std::endl;

        // write tagged vertex IDs and their tags
        for ( VertexIterator vIt = vBegin; vIt != vEnd; ++vIt ) {
            if ( (*vIt)->tag( ) != subdiv::surf::VERTEX_NOTAG ) {
                *tg << (*vIt)->index( ) << " "
                    << static_cast< int >( (*vIt)->tag() ) << std::endl;
            }
        }
        for ( FacetIterator f = fBegin; f != fEnd; ++f ) {
            VecVPtr newVerts;
            facetTrees_[ *f ]->addNewVertices( std::back_inserter( newVerts ) );
            for ( VecVPtrIter vIt = newVerts.begin( ); vIt != newVerts.end( ); ++vIt ) {
                if ( (*vIt)->tag() != subdiv::surf::VERTEX_NOTAG ) {
                    *tg << (*vIt)->index( ) << " "
                        << static_cast< int >( (*vIt)->tag() ) << std::endl;
                }
            }
        }
        
        // write tagged edges
        const EdgeMapIter teBegin = taggedEdges.begin( );
        const EdgeMapIter teEnd = taggedEdges.end( );
        for ( EdgeMapIter teIt = teBegin; teIt != teEnd; ++teIt ) {
            Edge * te = teIt->first;
            // vertex IDs and edge tag
            *tg << te->first( )->index( ) << " "
                << te->second( )->index( ) << " "
                << static_cast< int >( te->getTag( ) ) << " ";
            // extra integers
            std::copy( teIt->second.begin( ), teIt->second.end( ),
                       std::ostream_iterator< int >( *tg, " " ) );
            *tg << std::endl;
        }

        // clean-up
        for ( EdgeMapIter teIt = teBegin; teIt != teEnd; ++teIt )
            delete teIt->first;
        taggedEdges.clear( );
    }

    return;
}

//------------------------------------------------------------------------------
/// Update neighbor relations of the facets connected to the edge
template< typename V, typename E, typename F, typename FT > 
void subdiv::surf::Mesh< V, E, F, FT >::addOrMatch_( Vertex * v1, Vertex * v2, Facet * f )
{
    // create edge based on start and end node
    Edge * e = new Edge( v1, v2 ); 
    FTL_VERIFY( e != NULL );
    // associate facet from which the nodes originate
    e->setFacet( f, NULL );
    // set default edge tag
    e->setTag( subdiv::surf::EDGE_NOTAG );
    e->setCategory( subdiv::surf::ECAT_BOUNDARY );

    // (try to) insert into edge set
    const std::pair< EdgeSetConstIterator, bool > ei = edges_.insert( e );
    
    // edge already appeared in the EdgeSet (ei's bool is false)
    if ( not ei.second ) {
        // delete the duplicate edge
        delete e;
    
        // Given two edges that are the same
        Edge * eo = *(ei.first);
        FTL_VERIFY( eo != NULL );

        // Get first attached facet to edge
        Facet * fo = eo->getFacet( );

        // The edge was set at read-in and has not any facet attached
        if ( not fo ) {
            eo->setFacet( f, NULL );  // EDGE_NOTAG
            eo->setCategory( subdiv::surf::ECAT_BOUNDARY );
        }
        // Given two edges that are the same, set up the face adjacency and
        else if ( f != fo ) {
            // add facet to existing edge's facets
            eo->insertFacet( f, NULL );

            // check if "old face" already knows its opposite face
            if ( fo->opposite( eo->first(), eo->second() ) ) {
                // non-manifold geometry,
                // i.e. more than two faces join at an edge
                eo->setCategory( subdiv::surf::ECAT_NONMANIFOLD );

                // facets are stored in ring of opposite neighbors
                //    fo -->  f1 --> ... ---> fn
                //    ^                       |
                //    |                       |
                //    `-----------------------'
                // find facet #fn whose opposite is starting point #fo
                Facet * fn = fo->opposite( eo->first(), eo->second() );
                //std::cout<<fn->index()<<std::endl;
                while ( fn->opposite( eo->first(), eo->second() ) != fo )
                    fn = fn->opposite( eo->first(), eo->second() );
                //std::cout<<fn->index()<<std::endl;

                // insert #f between #fn and #fo
                //    fo -->  f1 --> ... ---> fn ---> f
                //    ^                               |
                //    |                               |
                //    `-------------------------------'
                FTL_VERIFY( f != fn );
                fn->opposite( eo->first(), eo->second() ) = f;
                f->opposite( eo->first(), eo->second() ) = fo;
            }
            else {
                // set neighbor information
                //    fo ---> f
                //    ^       |
                //    |       |
                //    `-------'
                eo->setCategory( subdiv::surf::ECAT_MANIFOLD );
                FTL_VERIFY( f != fo );

                fo->opposite( eo->first(), eo->second() ) = f;
                f->opposite( eo->first(), eo->second() ) = fo;
            }
        }
    }  // end of "edge already known"

    return;
}

//------------------------------------------------------------------------------
/// Update facet tree and vertex indices
template< typename V, typename E, typename F, typename FT > 
void subdiv::surf::Mesh< V, E, F, FT >::updateIndices_( unsigned & numVerticesMesh,
                                                        unsigned & numLeafs )
{
    typedef typename std::vector< Vertex * >        VecVPtr;
    typedef typename VecVPtr::iterator              VecVPtrIter;
    typedef typename std::vector< FTree * >         VecFTPtr;
    typedef typename VecFTPtr::iterator             VecFTPtrIter;

    // facet range
    const FacetIterator fBegin = facets_.begin();
    const FacetIterator fEnd   = facets_.end();

    // index vertices and count vertices and facets
    numVerticesMesh = vertices_.size( );
    numLeafs = 0;
    for ( FacetIterator f = fBegin; f != fEnd; ++f ) {
        if ( facetTrees_.find( *f ) != facetTrees_.end() ) {
            VecVPtr newVerts;
            facetTrees_[ *f ]->addNewVertices( std::back_inserter( newVerts ) );
            for ( VecVPtrIter v = newVerts.begin(); v != newVerts.end(); ++v ) {
                (*v)->setIndex( numVerticesMesh++ );
            }
            VecFTPtr leafs;
            facetTrees_[ *f ]->addLeafs( std::back_inserter( leafs ) ); 
            for ( VecFTPtrIter lIter = leafs.begin(); lIter != leafs.end(); ++lIter ) {
                (*lIter)->setIndex( numLeafs++ );
            }
        }
    }

    // done
    return;
}

//------------------------------------------------------------------------------
template< typename V, typename E, typename F, typename FT > 
void subdiv::surf::Mesh< V, E, F, FT >::orientateManifolds_( )
{
    // collect alignment status of each manifold
    std::vector< std::pair< bool, unsigned > > alManifolds;
    unsigned numFacets = 0;
    
    // collection of checked facets w.r.t. orientation
    std::map< Facet *, unsigned > alFacets;

    // loop facets
    for ( FacetIterator fIt = facets_.begin( ); fIt != facets_.end( ); ++fIt ) {
        // pointer to facet
        Facet * f = *fIt;
        
        // if facet has not yet been treated
        if ( alFacets.find( f ) == alFacets.end() ) {
            // entered NEW manifold
            bool manOrientation = false;
            alFacets.insert( std::make_pair( f, alManifolds.size() ) );

            // initialise collection of facets of manifold to be parsed
            // this also initialises the reference orientation of the new manifold
            std::set< Facet * > mFacets;
            mFacets.insert( f );

            // collect (and possibly reverse orientation) of facets
            // on current manifold
            while ( not mFacets.empty() ) {
                // collect new arrivals in facet collection which are neigbours
                // to current active set #mFacets
                std::vector< Facet * > nFacets;
                typedef typename std::set< Facet * >::iterator SetFPtrIter;
                for ( SetFPtrIter mIt = mFacets.begin(); mIt != mFacets.end(); ++mIt ) {
                    // facet on manifold
                    Facet * fm = *mIt;
                    // loop neigbours of latter
                    for ( unsigned n = 0; n < F::numNeighbors; ++n ) {
                        Facet * fn = fm->neighbor( n );
                        // only those connected over manifold edge
                        if ( fn and fm->isManifold( n ) ) {
                            // vertices on edge
                            const int vn0 = fn->localVertexIndex( fm->vertex( n ) );
                            const int vn1 = fn->localVertexIndex( fm->vertex( next_[ n ] ) );
                            // set manifold neighborhood
                            fm->manifoldNeighbor( n ) = fn;
                            if ( vn1 == next_[ vn0 ] )  // reversed oriented
                                fn->manifoldNeighbor( vn0 ) = fm;
                            else  // same orientation
                                fn->manifoldNeighbor( vn1 ) = fm;
                            // align
                            if ( alFacets.find( fn ) == alFacets.end( ) ) {
                                if ( vn1 == next_[ vn0 ] ) {  // reversed oriented
                                    fn->reverseOrientation( );
                                    manOrientation = true;  // set alignment status
                                }
                                nFacets.push_back( fn );
                                alFacets.insert( std::make_pair( fn, alManifolds.size( ) ) );
                            }

                        }
                    }
                }
                // update : copy newly arrived manifold facets to active set
                mFacets.clear();
                mFacets.insert( nFacets.begin(), nFacets.end() );
            }
            // store manifold and its alignment status
            alManifolds.push_back( std::make_pair( manOrientation,
                                                   alFacets.size( ) - numFacets ) );
            numFacets = alFacets.size( );
        }

    }
    FTL_VERIFY( facets_.size() == alFacets.size() );

    // inform user
    if ( verbose_ ) {
        const unsigned numAligned =
            std::count_if( alManifolds.begin( ), alManifolds.end( ),
                           std::bind( &std::pair< bool, unsigned >::pair::first, std::placeholders::_1 ) );

        std::cout << "Mesh contains " << alManifolds.size( )
                  << " manifold(s) of which " << numAligned
                  << " needed alignment" << std::endl;
        std::cout << "    Facets per manifold are ";
        std::transform( alManifolds.begin( ), alManifolds.end( ),
                        std::ostream_iterator< unsigned >( std::cout, " " ),
                        std::bind( &std::pair< bool, unsigned >::pair::second, std::placeholders::_1 ) );
        std::cout << std::endl;
    }
    
    return;
}

//------------------------------------------------------------------------------
/// Temp debugging job
template < typename V, typename E, typename F, typename FT > 
void subdiv::surf::Mesh< V, E, F, FT >::tempJob( ) 
{
    //a) print neighbors
    // FacetIterator  it, itb=fBegin(), ite=fEnd();    
    // for (it=itb; it!=ite; ++it) {
    //     std::cout<<(*it)->index()<<" | ";
        
    //     for (unsigned i = 0 ; i < numNeighbors; i++)
    //         if ((*it)->neighbor(i))
    //             std::cout<<(*it)->neighbor(i)->index()<<" ";
    //     std::cout<<std::endl;
    // }

    //b) call temp job in facet tree
    this -> iterateOverFacetTrees(std::mem_fun(&FTree::tempJob));

    return;
}
