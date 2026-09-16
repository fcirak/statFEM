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

//! @author Burkhard Bornemann, Kosala Bandara
//! @date   2010

#ifndef subdiv_apps_taggingfun_h 
#define subdiv_apps_taggingfun_h

#include <corlib/Constraints.hpp>

#include <Eigen/Core>
#include <corlib/eigenX.hpp>

//==============================================================================
// declarations
namespace app{

    template< typename VERTEX >
    void taggingVertex( VERTEX * v );

    template< typename FACET >
    void taggingEdge( FACET * f );

}


//==============================================================================
// definitions

//------------------------------------------------------------------------------
template< typename VERTEX >
void app::taggingVertex( VERTEX * v )
{
    namespace eigenX = corlib::eigenX;

    eigenX::VectorSd< 3 > x = v->giveCoordinates( );

    // bridge
    if ( corlib::fuzzyEqual( x[ 1 ], -5., 1.e-5 ) or
         corlib::fuzzyEqual( x[ 1 ],  5., 1.e-5 ) ) {
        const double rSq = x[ 0 ]*x[ 0 ] + x[ 2 ]*x[ 2 ];
        if ( corlib::fuzzyEqual( rSq, 100., 1.e-4 ) or
             corlib::fuzzyEqual( rSq,   0., 1.e-5 ) ) {
            v->setTag( subdiv::surf::VERTEX_CORNER );
        }
    }

    return;
}

//------------------------------------------------------------------------------
template< typename FACET >
void app::taggingEdge( FACET * f )
{
    namespace eigenX = corlib::eigenX;

    const enum corlib::shape myShape = FACET::myShape;
    const unsigned numEdges          = FACET::numEdges;
    const unsigned numVertices       = FACET::numVertices;
    const std::array< int, numVertices > & next = 
        subdiv::surf::ShapeProp< myShape >::next;

    typedef typename FACET::Vertex        Vertex;
    typedef eigenX::VectorSd< 3 >         Vec3;

    // loop facet edges
    for ( unsigned e = 0; e < numEdges; ++e ) {

        // information on vertices of edge
        Vertex * v0 = f->vertex( e );
        const Vec3 x0 = v0->giveCoordinates( );
        Vertex * v1 = f->vertex( next[ e ] );
        const Vec3 x1 = v1->giveCoordinates( );


        // "lateral kinks"
        if ( ( corlib::fuzzyEqual( x0[ 0 ], -10., 1.e-5 ) or
               corlib::fuzzyEqual( x0[ 0 ],  10., 1.e-5 ) ) and
             ( corlib::fuzzyEqual( x1[ 0 ], -10., 1.e-5 ) or
               corlib::fuzzyEqual( x1[ 0 ],  10., 1.e-5 ) ) ) {
            f->setEdgeTag( e, subdiv::surf::EDGE_CREASE );
            if ( v0->tag( ) == subdiv::surf::VERTEX_NOTAG )
                v0->setTag( subdiv::surf::VERTEX_CREASE );
            if ( v1->tag( ) == subdiv::surf::VERTEX_NOTAG )
                v1->setTag( subdiv::surf::VERTEX_CREASE );   
        }

    }

    return;
}
    
#endif


