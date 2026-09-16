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
template< typename FACET >
subdiv::surf::Walker< FACET >::Walker( Facet * f, Vertex * v ) :
    facet_( f ), vertex_( v )
{
    return;
}

//------------------------------------------------------------------------------
template< typename FACET >
subdiv::surf::Walker< FACET >::Walker( Facet * f, const int vi ) :
    facet_( f ), vertex_( NULL )
{
    FTL_VERIFY( facet_ );
    vertex_ = facet_->vertex( vi );
    return;
}

//------------------------------------------------------------------------------
template< typename FACET >
subdiv::surf::Walker< FACET >::~Walker( )
{
    vertex_ = NULL;
    facet_ = NULL;
    return;
}

//------------------------------------------------------------------------------
template< typename FACET >
void subdiv::surf::Walker< FACET >::set( Facet * f, const int vi )
{
    FTL_VERIFY( f );
    facet_ = f;
    vertex_ = facet_->vertex( vi );
    return;
}

//------------------------------------------------------------------------------
template< typename FACET >
subdiv::surf::Walker< FACET > subdiv::surf::Walker< FACET >::next( )
{
    const int localIndex = facet_->localVertexIndex( vertex_ );
    const int nextLocalIndex = ShapeProp< myShape >::next[ localIndex ];
    vertex_ = facet_->vertex( nextLocalIndex );
    return Walker< FACET >( facet_, vertex_ );
}

//------------------------------------------------------------------------------
template< typename FACET >
subdiv::surf::Walker< FACET > subdiv::surf::Walker< FACET >::previous( )
{
    const int localIndex = facet_->localVertexIndex( vertex_ );
    const int prevLocalIndex = ShapeProp< myShape >::previous[ localIndex ];
    vertex_ = facet_->vertex( prevLocalIndex );
    return Walker< FACET >( facet_, vertex_ );
}

//------------------------------------------------------------------------------
template< typename FACET >
subdiv::surf::Walker< FACET > subdiv::surf::Walker< FACET >::opposite( )
{
    const int localIndex = facet_->localVertexIndex( vertex_ );
    const int nextLocalIndex = ShapeProp< myShape >::next[ localIndex ];
    vertex_ = facet_->vertex( nextLocalIndex );
    facet_ = facet_->neighbor( localIndex );
    return Walker< FACET >( facet_, vertex_ );
}

//------------------------------------------------------------------------------
/// Move the walker in a combination of directions
template< typename FACET >
template< unsigned N >
void subdiv::surf::Walker< FACET >::move( const std::array< subdiv::surf::move, N > & dir )
{
    for ( unsigned i = 0; i < N; i++ ) {
        const subdiv::surf::move d = dir[ i ];

        const int localIndex = facet_->localVertexIndex( vertex_ );        
        const int nextLocalIndex = ShapeProp< myShape >::next[ localIndex ];
        const int prevLocalIndex = ShapeProp< myShape >::previous[ localIndex ];

        switch ( d ) {
        case subdiv::surf::NEXT:
            vertex_ = facet_->vertex( nextLocalIndex );
            break;
        case subdiv::surf::PREVIOUS:
            vertex_ = facet_->vertex( prevLocalIndex );
            break;
        case subdiv::surf::OPPOSITE:
            vertex_ = facet_->vertex( nextLocalIndex );
            facet_ = facet_->neighbor( localIndex );
            break;
        }        
    }
    return;
}
