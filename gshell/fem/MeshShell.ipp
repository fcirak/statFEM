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
template< typename GELEMENT, typename SHAPEFUN >
template< typename SMESH >
void gshell::fem::MeshShell< GELEMENT, SHAPEFUN >::initialise( SMESH & sMesh )
{
    // types
    typedef typename SMESH::Vertex              SVertex;
    typedef typename SMESH::Facet               SFacet;
    typedef typename SMESH::VertexIterator      SVertexConstIterator;
    typedef typename SMESH::FacetIterator       SFacetConstIterator;

    typedef std::map< unsigned, NodeType * >              NodesByIdMap;
    typedef typename NodesByIdMap::const_iterator         NodesByIdMapConstIter;

    typedef typename ShapeFun::SetVPtr                       SetVPtr;
    typedef typename ShapeFun::SetVPtrCIter                  SetVPtrConstIter;

    // copy nodes (create duplicates)
    NodesByIdMap nodesById;
    const SVertexConstIterator snb = sMesh.vBegin();
    const SVertexConstIterator sne = sMesh.vEnd();
    for ( SVertexConstIterator sn=snb; sn!=sne; ++sn ) {
        NodeType * gNode = new NodeType();
        FTL_VERIFY( gNode );
        gNode->setId( (*sn)->index() );
        gNode->setCoordinates( (*sn)->giveCoordinates() );
        nodes_.push_back( gNode );
        nodesById.insert( std::make_pair( gNode->giveId(), gNode ) );
    }

    // connect elements
    const SFacetConstIterator seb = sMesh.fBegin();
    const SFacetConstIterator see = sMesh.fEnd();
    for ( SFacetConstIterator se=seb; se!=see; ++se ) {
        // create subdivision shape function
    	ShapeFun sf( *se );
        // set supporting nodes
        SetVPtrConstIter pVBegin = sf.patchVerticesBegin();
        SetVPtrConstIter pVEnd   = sf.patchVerticesEnd();
        std::vector< NodeType * > supportNodes;
        for ( SetVPtrConstIter v=pVBegin; v!=pVEnd; ++v ) {
            const unsigned nodeId = (*v)->index();
            const NodesByIdMapConstIter gNodeIt = nodesById.find( nodeId );
            FTL_VERIFY( gNodeIt != nodesById.end() );
            NodeType * gNode = gNodeIt->second;
            FTL_VERIFY( gNode );
            supportNodes.push_back( gNode );
        }

        // create g-shell element
        ElementType * gElement = new ElementType();
        FTL_VERIFY( gElement );
        gElement->setSupportNodes( supportNodes );
        // set vertices (maybe not needed)
        for ( unsigned n=0; n<SFacet::numVertices; ++n ) {
            const unsigned nodeId = (*se)->vertex( n )->index();
            const NodesByIdMapConstIter gNodeIt = nodesById.find( nodeId );
            FTL_VERIFY( gNodeIt != nodesById.end() );
            NodeType * gNode = gNodeIt->second;
            FTL_VERIFY( gNode );
            gElement->setNodePtr( n, gNode );
        }

        // store g-shell element
        elements_.push_back( gElement );
    }

    return;
}

//------------------------------------------------------------------------------
template< typename GELEMENT, typename SHAPEFUN >
template< typename SMESH, typename QUADRATURE >
void gshell::fem::MeshShell< GELEMENT, SHAPEFUN >::populateShapeFunctionCache( SMESH & sMesh,
                                                                     const QUADRATURE & quadrature )
{
    // types
    typedef typename SMESH::Facet               SFacet;
    typedef typename SMESH::FacetIterator       SFacetConstIterator;

    // loop over surface facets and shell elements at the same time
    // and populate shape function evaluation points by
    // looping thru quadrature points
    const SFacetConstIterator seb = sMesh.fBegin();
    const SFacetConstIterator see = sMesh.fEnd();
    typename MeshBasic_::ElementVec_::iterator ge = elements_.begin();
    for ( SFacetConstIterator se = seb; se != see; ++se,++ge ) {
        // create subdivision shape function
        ShapeFun sf( *se );

        // create g-shell element
        ElementType * gElement = *ge;
        gElement->populateShapeFunctionCache( &sf, quadrature );
    }

    return;
}


//------------------------------------------------------------------------------
//! Free all dynamically allocated memory 
template< typename GELEMENT, typename SHAPEFUN >
gshell::fem::MeshShell< GELEMENT, SHAPEFUN >::~MeshShell()
{
    this -> clearLinks();
    return;
}

//------------------------------------------------------------------------------
template< typename GELEMENT, typename SHAPEFUN >
void gshell::fem::MeshShell< GELEMENT, SHAPEFUN >::readLinks( std::istream & is )
{
    // read header
    unsigned numLinksExternal, numLinksInternal;
    is >> numLinksExternal >> numLinksInternal;

    // read external links
    for ( unsigned i = 0; i < numLinksExternal; ++i ) {

        // get node and element
        unsigned nodeId, elementId;
        is >> nodeId >> elementId;

        NodeType * node = this -> getNodePointer( nodeId );
        FTL_VERIFY( node );

        ElementType * elem = this -> getElementPointer( elementId );
        FTL_VERIFY( elem );

        // create external link
        LinkExt * linkExt = new LinkExt( node, elem );
        // and store
        linksExt_.push_back( linkExt );

    }

    // read internal links
    for ( unsigned i = 0; i < numLinksInternal; ++i ) {

        // get node and element
        unsigned nodeId, elementId0, elementId1;
        is >> nodeId >> elementId0 >> elementId1;

        NodeType * node = this -> getNodePointer( nodeId );
        FTL_VERIFY( node );

        ElementType * elem0 = this -> getElementPointer( elementId0 );
        FTL_VERIFY( elem0 );
        ElementType * elem1 = this -> getElementPointer( elementId1 );
        FTL_VERIFY( elem1 );

        // create external link
        LinkInt * linkInt = new LinkInt( node, elem0, elem1 );
        // and store
        linksInt_.push_back( linkInt );

    }

    // done
    return;
}

//------------------------------------------------------------------------------
template< typename GELEMENT, typename SHAPEFUN >
void gshell::fem::MeshShell< GELEMENT, SHAPEFUN >::clearLinks()
{
    std::for_each( linksExt_.begin(), linksExt_.end(), corlib::deleteFunctor() );
    linksExt_.clear();

    std::for_each( linksInt_.begin(), linksInt_.end(), corlib::deleteFunctor() );
    linksInt_.clear();

    return;
}

//------------------------------------------------------------------------------
template< typename GELEMENT, typename SHAPEFUN >
void gshell::fem::MeshShell< GELEMENT, SHAPEFUN >::addLinkExt( const LinkExt & linkExt )
{
    linksExt_.push_back( new LinkExt( linkExt ) );
    return;
}

//------------------------------------------------------------------------------
template< typename GELEMENT, typename SHAPEFUN >
void gshell::fem::MeshShell< GELEMENT, SHAPEFUN >::addLinkInt( const LinkInt & linkInt )
{
    linksInt_.push_back( new LinkInt( linkInt ) );
    return;
}

//------------------------------------------------------------------------------
template< typename GELEMENT, typename SHAPEFUN >
template< typename OP >
OP gshell::fem::MeshShell< GELEMENT, SHAPEFUN >::iterateOverLinksExt( OP op )
{
    typename VecLinkExt_::iterator begin = linksExt_.begin();
    typename VecLinkExt_::iterator end   = linksExt_.end();
    return std::for_each( begin, end, op );
}

//------------------------------------------------------------------------------
template< typename GELEMENT, typename SHAPEFUN >
template< typename OP >
OP gshell::fem::MeshShell< GELEMENT, SHAPEFUN >::iterateOverLinksInt( OP op )
{
    typename VecLinkInt_::iterator begin = linksInt_.begin();
    typename VecLinkInt_::iterator end   = linksInt_.end();
    return std::for_each( begin, end, op );
}
