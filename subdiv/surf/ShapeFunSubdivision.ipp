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
// Tensor indices     \alpha\beta=  11,  12,  21,  22
// C indices                        00,  01,  10,  11
// Access   2*\alpha+\beta           0,   1    2    3
// 3-Voigt C-indices                 0    1    1    2
template< typename FACET, typename SUBDIV >
const unsigned subdiv::surf::ShapeFunSubdivision< FACET, SUBDIV >::
voigtForward[localDim][localDim] = { { 0, 1 }, { 1, 2 } };


//------------------------------------------------------------------------------
template< typename FACET, typename SUBDIV >
subdiv::surf::ShapeFunSubdivision< FACET, SUBDIV >::ShapeFunSubdivision(
    Facet * f, const bool independentPatch
    ) : facet_( f ), indiePatch_( NULL )
{

    // collect patch
    subdiv::surf::collectNRingFacets( numRings, facet_, patchFacets_, patchVertices_ );
    patchFacets_.insert( patchFacets_.begin( ), facet_ ); // at _front_
    patchVertices_.insert( facet_->verticesBegin( ), facet_->verticesEnd( ) );

    // if subdivision shape functions use extra rings (necessary for certain rules)
    if ( numExtraRings > 0 ) {
        const unsigned numTotalRings = numRings + numExtraRings;
        // collect 2-ring of #facet_
        subdiv::surf::collectNRingFacets( numTotalRings, facet_, extraFacets_, extraVertices_ );
        // remove facets contained in #patchFacets_
        const VecFPtrIter extraFacetsEnd =
            std::remove_if( extraFacets_.begin( ), extraFacets_.end( ),
                            std::bind( corlib::contains< VecFPtrIter, Facet * >,
                                         patchFacets_.begin( ), patchFacets_.end( ), std::placeholders::_1 ) );
        extraFacets_.erase( extraFacetsEnd, extraFacets_.end( ) );
        // remove vertices contained in #patchVertices_
        for ( SetVPtrIter v = patchVertices_.begin(); v != patchVertices_.end(); ++v )
            extraVertices_.erase( *v ); 
    }

    // make independent copy of local patch
    if ( independentPatch ) {
        this->liberatePatch_( );
    }

    return;
}



//------------------------------------------------------------------------------
template< typename FACET, typename SUBDIV >
void subdiv::surf::ShapeFunSubdivision< FACET, SUBDIV >::evaluate(
    const VecLDim & xi,
    Vec & phi
    )
{
    // determine subdivison matrix
    Mat subMat( numFunctions, patchVertices_.size( ) );
    VecLDim subXi;
    MatLDimLDim subJac;
    this->computeSubdivisionMatrix_( xi, subMat, subXi, subJac );

    // determine shape functions over patch
    VecNF subPhi;
    spline_.evaluate( subXi, subPhi );
    phi = subPhi * subMat;

    return;
}

//------------------------------------------------------------------------------
template< typename FACET, typename SUBDIV >
void subdiv::surf::ShapeFunSubdivision< FACET, SUBDIV >::evaluateGradient(
    const VecLDim & xi,
    Mat & dPhiDXi
    )
{
    // determine subdivison matrix
    Mat subMat;
    VecLDim subXi;
    MatLDimLDim subJac;
    this->computeSubdivisionMatrix_( xi, subMat, subXi, subJac );
    
    // determine 1st derivative shape functions over patch
    MatLDimNF subDPhiDXi;
    spline_.evaluateGradient( subXi, subDPhiDXi );
    dPhiDXi = ( subJac.traspose( ) * subDPhiDXi ) * subMat;
    
    return;
}

//------------------------------------------------------------------------------
template< typename FACET, typename SUBDIV >
void subdiv::surf::ShapeFunSubdivision< FACET, SUBDIV >::evaluateHessian(
    const VecLDim & xi,
    Mat & ddPhiDDXi
    )
{
    // determine subdivison matrix
    Mat subMat;
    VecLDim subXi;
    MatLDimLDim subJac;
    this->computeSubdivisionMatrix_( xi, subMat, subXi, subJac );
    
    // determine 2nd derivative of shape functions over patch
    MatVecNFLDimLDim subDDPhiDDXi;
    spline_.evaluateHessian( subXi, subDDPhiDDXi );
    MatSDimNF temp2;
    for ( unsigned f = 0; f < numFunctions; ++f ) {
        MatLDimLDim ddPhiDDXi_f;
        for ( unsigned i = 0; i < localDim; ++i )
            for ( unsigned j = 0; j < localDim; ++j )
                ddPhiDDXi_f( i, j ) = subDDPhiDDXi( i, j )( f );
        ddPhiDDXi_f = subJac.transpose( ) * ( ddPhiDDXi_f * subJac );
        // store
        temp2( voigtForward[0][0], f ) = ddPhiDDXi_f( 0, 0 );
        temp2( voigtForward[0][1], f ) = ddPhiDDXi_f( 0, 1 );
        temp2( voigtForward[1][1], f ) = ddPhiDDXi_f( 1, 1 );
    }
    ddPhiDDXi = temp2 * subMat;

    return;
}

//------------------------------------------------------------------------------
template< typename FACET, typename SUBDIV >
void subdiv::surf::ShapeFunSubdivision< FACET, SUBDIV >::evaluateGradHess(
    const VecLDim & xi,
    Vec & phi,
    Mat & dPhiDXi,
    Mat & ddPhiDDXi
    )
{
    // determine subdivison matrix
    Mat subMat;
    VecLDim subXi;
    MatLDimLDim subJac;
    this->computeSubdivisionMatrix_( xi, subMat, subXi, subJac );

    // determine shape functions over patch
    VecNF subPhi;
    spline_.evaluate( subXi, subPhi );
    phi = subPhi.transpose() * subMat;

    // determine 1st derivative of shape functions over patch
    MatLDimNF subDPhiDXi;
    spline_.evaluateGradient( subXi, subDPhiDXi );
    dPhiDXi = ( subJac.transpose( ) * subDPhiDXi ) * subMat;

    // determine 2nd derivative of shape functions over patch
    MatVecNFLDimLDim subDDPhiDDXi;
    spline_.evaluateHessian( subXi, subDDPhiDDXi );
    MatSDimNF temp2;
    for ( unsigned f = 0; f < numFunctions; ++f ) {
        MatLDimLDim ddPhiDDXi_f;
        for ( unsigned i = 0; i < localDim; ++i )
            for ( unsigned j = 0; j < localDim; ++j )
                ddPhiDDXi_f( i, j ) = subDDPhiDDXi( i, j )( f );
        ddPhiDDXi_f = subJac.transpose( ) * ( ddPhiDDXi_f * subJac );
        // store
        temp2( voigtForward[0][0], f ) = ddPhiDDXi_f( 0, 0 );
        temp2( voigtForward[0][1], f ) = ddPhiDDXi_f( 0, 1 );
        temp2( voigtForward[1][1], f ) = ddPhiDDXi_f( 1, 1 );
    }
    ddPhiDDXi = temp2 * subMat;

    return;
}

//------------------------------------------------------------------------------
template< typename FACET, typename SUBDIV >
void subdiv::surf::ShapeFunSubdivision< FACET, SUBDIV >::setPositions(
    const Eigen::MatrixXd & coords
    )
{
    FTL_VERIFY( coords.rows() == Vertex::dim );
    FTL_VERIFY( coords.cols() == patchVertices_.size( ) );
    for ( SetVPtrIter v = patchVertices_.begin( ); v != patchVertices_.end( ); ++v ) {
        const unsigned i = std::distance( patchVertices_.begin( ), v );
        const int lev = (*v)->maxActiveLevel( );
        Point_ coord = (*v)->point( lev );
        coord.head( Vertex::dim ) = coords.col( i );
        (*v)->addNewLevel( coord, lev );  // overwrites level
    }
    return;
}

//------------------------------------------------------------------------------
template< typename FACET, typename SUBDIV >
void subdiv::surf::ShapeFunSubdivision< FACET, SUBDIV >::evaluate(
    const unsigned vIndex,
    Vec & coeff
    )
{
    FTL_VERIFY( vIndex < numVertices );

    // create patch
    Mesh_ patch;
    this->formPatchMesh_( patch );

    // get corresponding facet tree
    FacetTree_ * ft = patch.facetTree( *(patch.fBegin( )) );
    // and its tangents
    typename Subdiv::MapVPtrD effCoeff;
    subdiv_.computeLimit( ft, vIndex, 0, effCoeff );

    // extract gradient w.r.t. vertex positions
    coeff.resize( this->numPatchVertices( ) );
    auto coeffIt = coeff.data(); // starting Eigen 3.4 it should be coeff.begin()
    typedef typename Mesh_::VertexIterator VertexIter;
    for ( VertexIter v = patch.vBegin( ); v != patch.vEnd( ); ++v, ++coeffIt ) {
        if ( effCoeff.find( *v ) != effCoeff.end( ) ) {
            *(coeffIt) = effCoeff.find( *v )->second;
        }
        else  {
            *(coeffIt) = 0.;
        }
    }

    // done
    return;
}


//------------------------------------------------------------------------------
template< typename FACET, typename SUBDIV >
void subdiv::surf::ShapeFunSubdivision< FACET, SUBDIV >::evaluateGradient(
    const unsigned vIndex,
    Mat & coeff
    )
{
    FTL_VERIFY( vIndex < numVertices );

    // create patch
    Mesh_ patch;
    this->formPatchMesh_( patch );

    // get corresponding facet tree
    FacetTree_ * ft = patch.facetTree( *(patch.fBegin( )) );
    // and its tangents
    typename Subdiv::MapVPtrD effCoeff0, effCoeff1;
    subdiv_.computeTangents( ft, vIndex, 0, effCoeff0, effCoeff1 );

    // extract gradient w.r.t. vertex positions
    coeff.resize( localDim, this->numPatchVertices() );
    typedef typename Mesh_::VertexIterator VertexIter;
    for ( VertexIter v = patch.vBegin( ); v != patch.vEnd( ); ++v ) {
        const unsigned i = std::distance( patch.vBegin( ), v );

        if ( effCoeff0.find( *v ) != effCoeff0.end( ) )
            coeff( 0, i ) = effCoeff0.find( *v )->second;
        else
            coeff( 0, i ) = 0.;

        if ( effCoeff1.find( *v ) != effCoeff1.end( ) )
            coeff( 1, i ) = effCoeff1.find( *v )->second;
        else
            coeff( 1, i ) = 0.;
    }

    // done
    return;
}


//------------------------------------------------------------------------------
template< typename FACET, typename SUBDIV >
void subdiv::surf::ShapeFunSubdivision< FACET, SUBDIV >::write( std::ostream & os )
{
    os << "Facet " << facet_->index( ) << std::endl;

    os << "    vertex IDs : ";
    std::transform( facet_->verticesBegin(), facet_->verticesEnd(),
                    std::ostream_iterator< unsigned >( os, " " ),
                    std::bind( &Vertex::index, std::placeholders::_1 ) );
    os << std::endl;

    os << "    patch vertex IDs : ";
    std::transform( patchVertices_.begin(), patchVertices_.end(),
                    std::ostream_iterator< unsigned >( os, " " ),
                    std::bind( &Vertex::index, std::placeholders::_1 ) );
    os << std::endl;

    os << "    patch facet IDs : ";
    std::transform( patchFacets_.begin(), patchFacets_.end(),
                    std::ostream_iterator< unsigned >( os, " " ),
                    std::bind( &Facet::index, std::placeholders::_1 ) );
    os << std::endl;

    os << "    neighbour facet IDs : ";
    corlib::transform_if( facet_->neighborsBegin( ), facet_->neighborsEnd( ),
                          std::ostream_iterator< unsigned >( os, " " ),
                          std::bind( &Facet::index, std::placeholders::_1 ),
                          std::bind( std::not_equal_to< Facet * >( ), std::placeholders::_1,
                                       static_cast< Facet * >( NULL ) ) );
    os << std::endl;

    os << "    all (non-)manifold neighbour IDs : ";
    for ( unsigned n = 0; n < numNeighbors; n++ ) {
        std::vector< Facet * > all;
        facet_->neighborsAtEdge( n, all );
        std::transform( all.begin(), all.end(),
                        std::ostream_iterator< unsigned >( os, " " ),
                        std::bind( &Facet::index, std::placeholders::_1 ) );
    }
    os << std::endl;

    os << "    manifold neighbour facet IDs : ";
    corlib::transform_if( facet_->manNeighborsBegin( ), facet_->manNeighborsEnd( ),
                          std::ostream_iterator< unsigned >( os, " " ),
                          std::bind( &Facet::index, std::placeholders::_1 ),
                          std::bind( std::not_equal_to< Facet * >( ), std::placeholders::_1,
                                       static_cast< Facet * >( NULL ) ) );
    os << std::endl;

    os << "    edge tags : ";
    for ( unsigned e = 0; e < numEdges; ++e ) {
        os << facet_->getEdgeTag( e ) << " ";
    }
    os << std::endl;

    os << "    extra vertex IDs : ";
    std::transform( extraVertices_.begin(), extraVertices_.end(),
                    std::ostream_iterator< unsigned >( os, " " ),
                    std::bind( &Vertex::index, std::placeholders::_1 ) );
    os << std::endl;

    os << "    extra facet IDs : ";
    std::transform( extraFacets_.begin(), extraFacets_.end(),
                    std::ostream_iterator< unsigned >( os, " " ),
                    std::bind( &Facet::index, std::placeholders::_1 ) );
    os << std::endl;

    return;
}

//------------------------------------------------------------------------------
template< typename FACET, typename SUBDIV >
void subdiv::surf::ShapeFunSubdivision< FACET, SUBDIV >::liberatePatch_( )
{
    // create patch copy
    indiePatch_ = new Mesh_( );
    this->formPatchMesh_( *indiePatch_ );

    // assign global node IDs to independent patch
    typename Mesh_::VertexIterator w = indiePatch_->vBegin();
    for ( SetVPtrIter v=patchVertices_.begin(); v!=patchVertices_.end(); ++v, ++w )
        (*w)->setIndex( (*v)->index() );
    for ( SetVPtrIter v=extraVertices_.begin(); v!=extraVertices_.end(); ++v, ++w )
        (*w)->setIndex( (*v)->index() );

    // assign pointers to patch copy
    facet_ = *(indiePatch_->fBegin());
    const unsigned numPatchFacets = patchFacets_.size();
    std::copy( indiePatch_->fBegin(), indiePatch_->fBegin()+numPatchFacets,
               patchFacets_.begin() );
    const unsigned numPatchVertices = patchVertices_.size();
    patchVertices_.clear();
    std::copy( indiePatch_->vBegin(), indiePatch_->vBegin()+numPatchVertices,
               std::inserter( patchVertices_, patchVertices_.begin() ) );
    FTL_VERIFY( patchVertices_.size() == numPatchVertices );

    if ( numExtraRings > 0 ) {
        std::copy( indiePatch_->fBegin()+numPatchFacets, indiePatch_->fEnd(),
                   extraFacets_.begin() );
        const unsigned numExtraVertices = extraVertices_.size();
        extraVertices_.clear();
        std::copy( indiePatch_->vBegin()+numPatchVertices, indiePatch_->vEnd(),
                   std::inserter( extraVertices_, extraVertices_.begin() ) );
        FTL_VERIFY( extraVertices_.size() == numExtraVertices );
    }
}

//------------------------------------------------------------------------------
template< typename FACET, typename SUBDIV >
void subdiv::surf::ShapeFunSubdivision< FACET, SUBDIV >::writePatch_(
    std::ostream & smf,
    std::ostream & tg
    )
{
    // write head
    corlib::SmfHead smfHead;
    smfHead.setElementShape( myShape );
    smfHead.setElementNumPoints( numVertices );
    smfHead.write( smf );
    smf << patchVertices_.size( ) + extraVertices_.size() << " "
        << patchFacets_.size( ) + extraFacets_.size( ) << std::endl;

    // write vertex co-ordinates
    for ( SetVPtrIter v = patchVertices_.begin( ); v != patchVertices_.end( ); ++v ) {
        smf << *v << std::endl;
    }
    for ( SetVPtrIter v = extraVertices_.begin( ); v != extraVertices_.end( ); ++v ) {
        smf << *v << std::endl;
    }

    // write facet connectivity
    for ( VecFPtrIter f = patchFacets_.begin(); f != patchFacets_.end(); ++f ) {
        for ( unsigned i = 0; i < numVertices; ++i ) {
            Vertex * v = (*f)->vertex( i );
            FTL_VERIFY( patchVertices_.find( v ) != patchVertices_.end() );
            smf << std::distance( patchVertices_.begin(), patchVertices_.find( v ) ) << " ";
        }
        smf << std::endl;
    }
    for ( VecFPtrIter f = extraFacets_.begin(); f != extraFacets_.end(); ++f ) {
        for ( unsigned i = 0; i < numVertices; ++i ) {
            Vertex * v = (*f)->vertex( i );
            unsigned vId;
            if ( patchVertices_.find( v ) != patchVertices_.end() )
                vId = std::distance( patchVertices_.begin(), patchVertices_.find( v ) );
            else
                vId = patchVertices_.size() +
                    std::distance( extraVertices_.begin(), extraVertices_.find( v ) );
            smf << vId << " ";
        }
        smf << std::endl;
    }

    // collect only tagged vertices
    typedef std::vector< std::pair< Vertex *, unsigned > > VertexIndexVector;
    VertexIndexVector taggedVertices;
    for ( SetVPtrIter v = patchVertices_.begin( ); v != patchVertices_.end( ); ++v ) {
        if ( (*v)->tag() != subdiv::surf::VERTEX_NOTAG ) {
            const unsigned vId = std::distance( patchVertices_.begin( ), v );
            taggedVertices.push_back( std::make_pair( *v, vId ) );
        }
    }
    for ( SetVPtrIter v = extraVertices_.begin( ); v != extraVertices_.end( ); ++v ) {
        if ( (*v)->tag() != subdiv::surf::VERTEX_NOTAG ) {
            const unsigned vId = patchVertices_.size() +
                std::distance( extraVertices_.begin( ), v );
            taggedVertices.push_back( std::make_pair( *v, vId ) );
        }
    }

    // collect only tagged edges
    typedef subdiv::surf::Edge< Vertex, Facet >                   Edge;
    typedef typename Edge::EdgeTag                                EdgeTag;
    typedef typename Edge::Less                                   EdgeLess;
    typedef std::multiset< Edge *, EdgeLess >                     EdgeMSet;
    typedef typename EdgeMSet::iterator                           EdgeMSetIter;

    EdgeMSet taggedEdges;
    for ( VecFPtrIter f = patchFacets_.begin(); f != patchFacets_.end(); ++f ) {
        for ( unsigned e = 0; e < numEdges; ++e ) {
            if ( (*f)->getEdgeTag( e ) != subdiv::surf::EDGE_NOTAG ) {
                Vertex * v0 = (*f)->vertex( e );
                Vertex * v1 = (*f)->vertex( subdiv::surf::ShapeProp< myShape >::next[ e ] );
                Edge * edge = new Edge( v0, v1 );
                const EdgeTag eTag = (*f)->getEdgeTag( e );
                // revert "smooth crease" to crease, if manifold neighbour is not in patch
                if ( eTag == subdiv::surf::EDGE_SMOOTHCREASE ) {
                    Facet * fn = (*f)->manifoldNeighbor( e );
                    if ( fn ) {
                        if ( corlib::contains( patchFacets_.begin(), patchFacets_.end(), fn ) or
                             corlib::contains( extraFacets_.begin(), extraFacets_.end(), fn ) ) {  // "smooth crease"
                            edge->setTag( eTag );
                            edge->setFacet( *f, fn );
                        }
                        else {                                 // "smooth crease" but reaching out of patch
                            edge->setTag( subdiv::surf::EDGE_CREASE );
                            edge->setFacet( *f, NULL );
                        }
                    }
                    else {                                    // true crease attached to "smooth crease"
                        edge->setTag( subdiv::surf::EDGE_CREASE );
                        edge->setFacet( *f, NULL );
                    }
                }
                else {                                        // ordinary true crease
                    edge->setTag( eTag );
                    edge->setFacet( *f, NULL );
                }
                // store edge
                taggedEdges.insert( edge  );
            }
        }
    }

    for ( VecFPtrIter f = extraFacets_.begin(); f != extraFacets_.end(); ++f ) {
        for ( unsigned e = 0; e < numEdges; ++e ) {
            if ( (*f)->getEdgeTag( e ) != subdiv::surf::EDGE_NOTAG ) {
                Vertex * v0 = (*f)->vertex( e );
                Vertex * v1 = (*f)->vertex( subdiv::surf::ShapeProp< myShape >::next[ e ] );
                Edge * edge = new Edge( v0, v1 );
                const EdgeTag eTag = (*f)->getEdgeTag( e );
                // revert "smooth crease" to crease, if manifold neighbour is not in patch
                if ( eTag == subdiv::surf::EDGE_SMOOTHCREASE ) {
                    Facet * fn = (*f)->manifoldNeighbor( e );
                    if ( fn ) {
                        if ( corlib::contains( patchFacets_.begin(), patchFacets_.end(), fn ) or
                             corlib::contains( extraFacets_.begin(), extraFacets_.end(), fn ) ) {  // "smooth crease"
                            edge->setTag( eTag );
                            edge->setFacet( *f, fn );
                        }
                        else {                                 // "smooth crease" but reaching out of patch
                            edge->setTag( subdiv::surf::EDGE_CREASE );
                            edge->setFacet( *f, NULL );
                        }
                    }
                    else {                                    // true crease attached to "smooth crease"
                        edge->setTag( subdiv::surf::EDGE_CREASE );
                        edge->setFacet( *f, NULL );
                    }
                }
                else {                                        // ordinary true crease
                    edge->setTag( eTag );
                    edge->setFacet( *f, NULL );
                }
                // store edge
                taggedEdges.insert( edge );
            }
        }
    }

    // write tag file
    tg << "TG" << std::endl
       << taggedVertices.size() << " " << taggedEdges.size() << std::endl;
    // ... here come the vertex tags
    typedef typename VertexIndexVector::iterator VertexIndexVectorIter;
    typedef typename Vertex::VertexTag           VertexTag;
    for ( VertexIndexVectorIter v = taggedVertices.begin( ); v != taggedVertices.end( ); ++v ) {
        VertexTag vTag = (*v).first->tag( );
        if ( vTag == subdiv::surf::VERTEX_SMOOTHCORNER or
             vTag == subdiv::surf::VERTEX_SMOOTHCREASE ) {
            Vertex * vPtr = (*v).first;
            // find (tagged) edges which contain vertex
            // edges are _not_ unique
            unsigned numSmoothCreases = 0;
            unsigned numCreases = 0;
            for ( EdgeMSetIter te = taggedEdges.begin(); te != taggedEdges.end(); ++te ) {
                // get range of coincident edges in multi-set
                std::pair< EdgeMSetIter, EdgeMSetIter > coinEdges =
                    taggedEdges.equal_range( *te );
                if ( (*te)->first() == vPtr or (*te)->second() == vPtr ) {
                    bool hasSmoothCrease = false;
                    bool hasCrease = false;
                    for ( EdgeMSetIter cte = coinEdges.first; cte != coinEdges.second; ++cte ) {
                        if ( (*cte)->getTag( ) == subdiv::surf::EDGE_SMOOTHCREASE )
                            hasSmoothCrease = true;
                        else if ( (*cte)->getTag( ) == subdiv::surf::EDGE_CREASE )
                            hasCrease = true;
                    }
                    if      ( hasSmoothCrease )  numSmoothCreases += 1;
                    else if ( hasCrease )        numCreases += 1;
                }
                te = --(coinEdges.second);
            }
            assert( ( numCreases + numSmoothCreases ) > 0 );
            
            // modify vertex tag acc. to edge environment
            if ( vTag == subdiv::surf::VERTEX_SMOOTHCORNER ) {
                if ( numSmoothCreases == 0 )
                    vTag = subdiv::surf::VERTEX_CORNER;
            }
            else if ( vTag == subdiv::surf::VERTEX_SMOOTHCREASE ) {
                if ( numSmoothCreases < 2 )
                    vTag = subdiv::surf::VERTEX_CREASE;
            }
        }

        // write
        tg << (*v).second << " "
           << static_cast< int >( vTag ) << std::endl;
    }
    // ... followed by the edge tags
    for ( EdgeMSetIter e = taggedEdges.begin(); e != taggedEdges.end(); ++e ) {
        Vertex * v0 = (*e)->first( );
        Vertex * v1 = (*e)->second( );
        unsigned iv0, iv1;
        if ( numExtraRings > 0 ) {
            iv0 = patchVertices_.find( v0 ) != patchVertices_.end( ) ?
                std::distance( patchVertices_.begin(), patchVertices_.find( v0 ) ) :
                patchVertices_.size( ) +
                std::distance( extraVertices_.begin(), extraVertices_.find( v0 ) );
            iv1 = patchVertices_.find( v1 ) != patchVertices_.end( ) ?
                std::distance( patchVertices_.begin(), patchVertices_.find( v1 ) ) :
                patchVertices_.size( ) +
                std::distance( extraVertices_.begin(), extraVertices_.find( v1 ) );
        }
        else {
            iv0 =
                std::distance( patchVertices_.begin(), patchVertices_.find( v0 ) );
            iv1 =
                std::distance( patchVertices_.begin(), patchVertices_.find( v1 ) );
        }

        // write
        tg << iv0 << " "
           << iv1 << " "
           << static_cast< int >( (*e)->getTag( ) );
        if ( (*e)->getTag( ) == subdiv::surf::EDGE_SMOOTHCREASE ) {
            Facet * f0 = ((*e)->facetsBegin( ))->first;
            Facet * f1 = ((*e)->facetsBegin( ))->second;
            assert( f1 );
            unsigned if0, if1;
            if ( numExtraRings > 0 ) {
                if0 = 
                    corlib::contains( patchFacets_.begin(), patchFacets_.end(), f0 ) ?
                    corlib::position( patchFacets_.begin(), patchFacets_.end(), f0 ) :
                    patchFacets_.size( ) +
                    corlib::position( extraFacets_.begin(), extraFacets_.end(), f0 );
                if1 = 
                    corlib::contains( patchFacets_.begin(), patchFacets_.end(), f1 ) ?
                    corlib::position( patchFacets_.begin(), patchFacets_.end(), f1 ) :
                    patchFacets_.size( ) +
                    corlib::position( extraFacets_.begin(), extraFacets_.end(), f1 );
            }
            else {
                if0 =
                    corlib::position( patchFacets_.begin(), patchFacets_.end(), f0 );
                if1 =
                    corlib::position( patchFacets_.begin(), patchFacets_.end(), f1 );
            }
            tg << " " << if0 << " " << if1;
        }
        tg << std::endl;
    }

    // clean-up taggedEdges
    std::for_each( taggedEdges.begin(), taggedEdges.end(), corlib::deleteFunctor( ) );
    taggedEdges.clear( );

    // done
    return;
}

//------------------------------------------------------------------------------
template< typename FACET, typename SUBDIV >
unsigned subdiv::surf::ShapeFunSubdivision< FACET, SUBDIV >::neededSubdivisions_(
    const VecLDim & xi
    )
{
    // determine minimal distance to facet boundary
    const double minDistance = subdiv::surf::shapeMinDistanceToBoundary< myShape >( xi );

    // number of needed subdivisions n such that xi lies within regular patch:
    //    (1/2)^n  <  minDistance(xi)
    return static_cast< unsigned >( std::floor( log2( 1. / minDistance ) ) + 1 );
}

//------------------------------------------------------------------------------
template< typename FACET, typename SUBDIV >
void subdiv::surf::ShapeFunSubdivision< FACET, SUBDIV >::formPatchMesh_(
    Mesh_ & patch
    )
{
    // generate patch streams
    std::stringstream smf, tg;
    this->writePatch_( smf, tg );  // this contains #patchFacets_ and #extraFacets_
    smf.seekg( 0 );
    tg.seekg( 0 );

    // generate patch mesh
    patch.readSmf( smf );
    patch.readTags( tg );
    patch.buildFacetTopology( );
    //patch.iterateOverFacets( std::bind( &Facet::write, std::placeholders::_1, std::ref( std::cout ) ) );
    patch.implantFacetTrees( );
    //patch.iterateOverFacetTrees( std::bind( &FacetTree_::write, std::placeholders::_1, std::ref( std::cout ) ) );

    return;
}

//------------------------------------------------------------------------------
template< typename FACET, typename SUBDIV >
void subdiv::surf::ShapeFunSubdivision< FACET, SUBDIV >::computeSubdivisionMatrix_(
    const VecLDim & xi,
    Mat & subMat,
    VecLDim & subXi,
    MatLDimLDim & subJac
    )
{
    // enlarge Vertex's point size to hold generalised co-ordinates for picking matrix
    const unsigned oldNumAddCoords = Vertex::getNumAddCoords( );
    Vertex::setNumAddCoords( patchVertices_.size() + extraVertices_.size() );

    // generate patch mesh
    Mesh_ patch;
    this->formPatchMesh_( patch );

#ifdef SHAPEFUN_DEBUG_WRITEPATCHES
    const std::pair< unsigned, unsigned > writePatchesRange( 24, 30 );
    const std::string patchesFolder( "patches" );
    if ( facet_->index() >= writePatchesRange.first and
         facet_->index() < writePatchesRange.second ) {
        std::stringstream makePatchesFolder;
        makePatchesFolder << "mkdir -p " << patchesFolder;
        FTL_VERIFY( ! system( makePatchesFolder.str().c_str() ) );

        std::stringstream smfName;
        smfName << patchesFolder << "/" << "patch-"
                << std::setfill('0') << std::setw(4) << facet_->index()
                << ".0.smf";
        std::ofstream patchSmf( smfName.str().c_str() );
        std::stringstream tgName;
        tgName << patchesFolder << "/"<< "patch-"
               << std::setfill('0') << std::setw(4) << facet_->index()
               << ".0.tg";
        std::ofstream patchTg( tgName.str().c_str() );

        patch.writeMeshWithIndex( patchSmf, &patchTg );

        patchSmf.close();
        patchTg.close();
    }
#endif

    // sanity checks
    const unsigned patchNumFacets = std::distance( patch.fBegin( ), patch.fEnd( ) );
    FTL_VERIFY( patchNumFacets == patchFacets_.size( ) + extraFacets_.size() );
    const unsigned patchNumVertices = std::distance( patch.vBegin( ), patch.vEnd( ) );
    FTL_VERIFY( patchNumVertices == patchVertices_.size( ) + extraVertices_.size() );
    FTL_VERIFY_DESCRIPTIVE( patchVertices_.size( ) <= Vertex::getNumAddCoords(),
                            "With patchNumVertices=%d, Vertex::pointSize_=%d+%d\n",
                            patchNumVertices, Vertex::dim, Vertex::getNumAddCoords() );

    // determine number of needed subdivisions
    const unsigned numSub = this->neededSubdivisions_( xi );

    // loop over vertices and set Kronecker delta in additional variable slots
    for ( typename Mesh_::VertexIterator v = patch.vBegin(); v != patch.vEnd(); ++v ) {
        Point_ extCoord = (*v)->point( 0 );
        for ( unsigned i = Vertex::dim; i < extCoord.size(); ++i ) {
            if ( i - Vertex::dim == (*v)->index() )
                extCoord( i ) = 1.;
            else
                extCoord( i ) = 0.;
        }
        (*v)->addNewLevel( extCoord, 0 );
    }
    
    // perform subdivision
    for ( unsigned k = 0; k < numSub; ++k ) {
        // subdivide all elements
        for ( typename Mesh_::FacetIterator f = patch.fBegin(); f != patch.fEnd(); ++f ) {
            patch.subdivideFacetTree( *f, &subdiv_ );
        }
#ifdef SHAPEFUN_DEBUG_WRITEPATCHES
        if ( facet_->index() >= writePatchesRange.first and
             facet_->index() < writePatchesRange.second ) {
            std::stringstream smfName;
            smfName << patchesFolder << "/" << "patch-"
                    << std::setfill('0') << std::setw(4) << facet_->index()
                    << "." << k+1 << ".smf";
            std::ofstream patchSmf( smfName.str().c_str() );
            std::stringstream tgName;
            tgName << patchesFolder << "/"<< "patch-"
                   << std::setfill('0') << std::setw(4) << facet_->index()
                   << "." << k+1 << ".tg";
            std::ofstream patchTg( tgName.str().c_str() );
            
            patch.writeMeshWithIndex( patchSmf, &patchTg );
            
            patchSmf.close();
            patchTg.close();
        }
#endif
    }

    // find child facet in which evaluation point is located
    // REMEMBER: first facet is facet of interest (the centre facet)
    FacetTree_ * subFacet = patch.facetTree( *patch.fBegin() );
    subXi = xi;
    subJac = Eigen::MatrixXd::Identity( localDim, localDim );
    for ( unsigned k = 0; k < numSub; ++k ) {
        unsigned child;
        subdiv::surf::ShapeChildCoordinate< myShape, Facet >()( facet_, subXi, subJac, child );
        subFacet = subFacet->child( child );
    }
    FTL_VERIFY( subFacet );

    // collect regular vertex neighbourhood
    typedef typename std::array<Vertex *, numFunctions>     VecVPtrNF;
    typedef typename VecVPtrNF::iterator                    VecVPtrNFIter;
    VecVPtrNF regSubVertices;
    subdiv::surf::collectManifoldRegularVertices< FacetTree_ >( subFacet, regSubVertices );

    // extract subdivision matrix
    subMat.resize( numFunctions, patchVertices_.size( ) );
    for ( unsigned i = 0; i < numFunctions; ++i ) {
        const int lev = regSubVertices[ i ]->maxActiveLevel( );
        Point_ coord = regSubVertices[ i ]->point( lev );
        for ( unsigned d = 0; d < patchVertices_.size( ); ++d )
            subMat( i, d ) = coord( Vertex::dim + d );
    }
    // each row must sum to 1
    FTL_VERIFY( corlib::fuzzyEqual( corlib::sumOfEntries( subMat ), numFunctions ) );

#ifdef SHAPEFUN_DEBUG_WRITEPATCHES
    // write shape functions
    if ( facet_->index() >= writePatchesRange.first and
         facet_->index() < writePatchesRange.second ) {

        // collect leafs
        typedef typename std::vector< FacetTree_ * >    VecFTPtr;
        typedef typename VecFTPtr::iterator             VecFTPtrIter;

        VecFTPtr leafs;
        for ( typename Mesh_::FacetIterator f = patch.fBegin(); f != patch.fEnd(); ++f ) {
            FacetTree_ * ft = patch.facetTree( *f );
            ft->addLeafs( std::back_inserter( leafs ) ); 
        }

        // collect vertices
        typedef typename std::set< Vertex * >           SetVPtr;
        typedef typename std::set< Vertex * >::iterator SetVPtrIterator;

        SetVPtr vertices;
        for ( VecFTPtrIter l = leafs.begin(); l != leafs.end(); ++l ) {
            vertices.insert( (*l)->verticesBegin( ),
                             (*l)->verticesEnd( ) );
        }

        // loop vertices
        for ( unsigned d = 0; d < patchVertices_.size( ); ++d ) {
            const unsigned i = Vertex::dim + d;

            // open file
            std::stringstream funName;
            funName << patchesFolder << "/" << "patch-"
                    << std::setfill('0') << std::setw(4) << facet_->index()
                    << "." << numSub << ".fun" << d;
            std::ofstream patchFun( funName.str().c_str() );
            patchFun << "TG" << std::endl
                     << vertices.size() << " 0" << std::endl;
            
            // loop nodes
            for ( SetVPtrIterator w = vertices.begin(); w != vertices.end(); ++w ) {
                const int lev = (*w)->maxActiveLevel( );
                patchFun << std::distance( vertices.begin(), w ) << " "
                         << (*w)->point( lev )[ i ] << std::endl;
            }
            
            // close file
            patchFun.close();
            
        }
    }
#endif
    
    // reset number of additional variables on vertices
    Vertex::setNumAddCoords( oldNumAddCoords );

    return;
}
