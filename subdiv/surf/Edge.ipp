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
template< typename V, typename F >
subdiv::surf::Edge< V, F >::Edge( VertexType * v1, VertexType * v2 ) :
    v1_( v1 ),
    v2_( v2 ),
    cat_( subdiv::surf::ECAT_UNKNOWN ),
    tag_( subdiv::surf::EDGE_NOTAG )
{
    facets_.push_back( std::make_pair( static_cast< FacetType * >( NULL ),
                                       static_cast< FacetType * >( NULL ) ) );
    return;
}

//------------------------------------------------------------------------------
template< typename V, typename F >
void subdiv::surf::Edge< V, F >::insertFacet( FacetType * f, FacetType * fm,
                                            bool overwrite )
{
    FacetPairIter fp =
        std::find_if( facets_.begin( ), facets_.end( ),
                      std::bind( std::equal_to< FacetType * >(),
                                 std::bind( &FacetPair::pair::first, 
                                 std::placeholders::_1 ),
                                 f ) );
    if ( fp == facets_.end( ) )
        facets_.push_back( FacetPair( f, fm ) );
    else if ( overwrite )
        fp->second = fm;
    return;
}

//------------------------------------------------------------------------------
template< typename V, typename F >
std::ostream & subdiv::surf::Edge< V, F >::write( std::ostream & os ) const
{
    os << "Edge" << std::endl;
        
    os << "    vertex IDs : " << v1_->index( ) << " " << v2_->index( ) << std::endl;

    os << "    category : " << cat_ << std::endl;

    os << "    tag : " << tag_ << std::endl;
        
    os << "    attached facets and thier manifold facets : " << std::endl;
    for ( unsigned f = 0; f < facets_.size(); ++f ) {
        os << "        " << facets_[ f ].first->index( );
        if ( facets_[ f ].second )
            os << " " << facets_[ f ].second->index( )  << std::endl;
        else
            os << " NULL" << std::endl;
    }

    return os;
}
