// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//
//                 Computational Structural Mechanics Lab
//                         University of Cambridge
//
// This software is copyrighted and all rights are retained by the University of
// Cambridge. All use, disclosure, and/or reproduction of any part not expressly
// authorised by the University of Cambridge is prohibited. (C) 2010
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

//! @file ConcentratedLoad.ipp

//------------------------------------------------------------------------------
//! Constructor given the function object representing the body force
template< typename ELEMENT >
template< typename MESH >
beam::fem::ConcentratedLoad< ELEMENT >::ConcentratedLoad( std::istream& inp,
                                                          const MESH& mesh )
  : commentChar_( '#' ),
    factor_( 1.0 )
{
    // skip comment lines
    skipCommentLines( inp );

    // get number of constraints
    unsigned numConLoads = 0;
    inp >> numConLoads;

    // go through constraints
    for ( unsigned c=0; c<numConLoads; ++c ) {

        // skip comment lines
        skipCommentLines( inp );

        // read in data of node
        unsigned nodeNum;
        inp >> nodeNum;
        Vec3_ forces;
        for ( unsigned i=0; i<3; ++i )
            inp >> forces[ i ];
        Vec3_ torques;
        for ( unsigned i=0; i<3; ++i )
            inp >> torques[ i ];

        // node range
        typedef typename MESH::NodeConstIterator NodeConstIterator;
        const NodeConstIterator nb = mesh.nodesBegin();
        const NodeConstIterator ne = mesh.nodesEnd();
            
        // find node pointer
        NodeType* np = NULL;
        for ( NodeConstIterator n=nb; n!=ne; ++n ) {
            if ( (*n)->giveId() == nodeNum ) {
                np = *n;
                break;
            }
        }
        FTL_VERIFY( np!=NULL &&
                    "Could not find node with given node ID" );

        // element range
        typedef typename MESH::ElementConstIterator ElementConstIterator;
        const ElementConstIterator eb = mesh.elementsBegin();
        const ElementConstIterator ee = mesh.elementsEnd();
            
        // find element of which node is a vertex
        ElementType* ep = NULL;
        VecLDim_ locCoor;
        for ( ElementConstIterator e=eb; e!=ee; ++e ) {
            typename ElementType::NodeIterator nIter = (*e)->nodesBegin();
            if ( *nIter == np ) {
                ep = (*e);
                locCoor = Eigen::VectorXd::Zero( localDim_ );
                break;
            }
            if ( *(nIter+1) == np ) {
                ep = (*e);
                locCoor = Eigen::VectorXd::Constant( localDim_, 1.0 );
                break;
            }
        }
        FTL_VERIFY( ep!=NULL &&
                    "Could not find element for node" );

        // store concentrated load
        loads_.insert( std::make_pair( ep, std::make_tuple( locCoor, forces, torques ) ) );

    }
    FTL_VERIFY( std::distance( loads_.begin(), loads_.end() )==static_cast< int >( numConLoads ) &&
                "Number of concentrated loads does not match" );

    return;
}

//------------------------------------------------------------------------------
//! Skip comment lines
template< typename ELEMENT >
void beam::fem::ConcentratedLoad< ELEMENT >::skipCommentLines( std::istream& inp ) const
{
    bool comment = true;
    while ( comment ) {
        char next;
        inp >> next;
        if ( next == commentChar_ ) {
            std::string commentLine;
            std::getline( inp, commentLine );
        }
        else {
            comment = false;
            inp.putback( next );
        }
    }
    return;
}

//------------------------------------------------------------------------------
//! Basic operation to add concentrated loads onto nodal forces via element call
template< typename ELEMENT >
void beam::fem::ConcentratedLoad< ELEMENT >::operator()( ELEMENT * ep )
{
    // find element in map of concentrated loads
    const LoadMapIter_ l = loads_.find( ep );

    // add concentrated loads
    if ( l != loads_.end() ) {
        ep->concentratedLoad( std::get< 0 >( l->second ),
                              std::get< 1 >( l->second ),
                              std::get< 2 >( l->second ),
                              factor_ );
    }

    return;
}

//------------------------------------------------------------------------------
template< typename ELEMENT >
std::ostream& beam::fem::ConcentratedLoad< ELEMENT >::print( std::ostream& os ) const
{
    // range
    const LoadMapConstIter_ lb = loads_.begin();
    const LoadMapConstIter_ le = loads_.end();
        
    // loop loads and print content
    for ( LoadMapConstIter_ el = lb; el != le; ++el ) {
        const unsigned loc =
            static_cast< unsigned >( std::get< 0 >( el -> second )[ 0 ] );
        
        typename ELEMENT::NodeConstIterator iter = (el -> first) -> nodesBegin();
        std::advance( iter, loc );

        os << (*iter) -> giveId() << " : "
           << std::get< 0 >( el->second ) << " "
           << std::get< 1 >( el->second ) << " "
           << std::get< 2 >( el->second );
        os << std::endl;
    }
    return os;
}
