// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//
//                 Computational Structural Mechanics Lab
//                         University of Cambridge
//
// This software is copyrighted and all rights are retained by the CSMLab. All
// use, disclosure, and/or reproduction of any part not expressly authorized by
// F Cirak is prohibited. (C) 2011.
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//

//! @author Kosala Bandara, Burkhard Bornemann
//! @date   2011


#ifndef subdiv_curve_walker_h
#define subdiv_curve_walker_h

#include <cstdlib>
#include <iterator>
#include <functional>
#include <algorithm>

#include <corlib/verify.hpp>

#include <subdiv/curve/ShapeDim.hpp>

//------------------------------------------------------------------------------
namespace subdiv{
    namespace curve{

        //! \brief Enumerator for walker movement
        enum move {
            NEXT,               //!< next
            PREVIOUS,           //!< previous
            OPPOSITE            //!< opposite
        };

        template< typename EDGE >
        class Walker;
    }
}

//==============================================================================
/// Walker navigates over the mesh.
///
/// It stores a edge and one of its vertices and updates this pair
/// due to #next(), #previous() and #opposite() steps.
///
/// These operations rely on a <i>consistently oriented</i> faces.
///
/// <i>Walk, don't run.</i>
///
/// \tparam  EDGE         Edge (tree) type
template< typename EDGE >
class subdiv::curve::Walker
{
public:
    typedef EDGE                                           Edge;
    typedef typename Edge::Vertex                          Vertex;

    static const corlib::shape myShape      = Edge::myShape;

public:
    /// Constructor
    ///
    /// \param[in]    e         Pointer to edge
    /// \param[in]    v         Pointer to vertex of edge
    Walker( Edge * e, Vertex * v ): edge_( e ), vertex_( v ) { }
    
    /// Constructor with local vertex index
    ///
    /// \param[in]    f         Pointer to edge
    /// \param[in]    vi        Edge-wise vertex index
    Walker( Edge * e, const int vi ): edge_( e ), vertex_( NULL )
    {
        FTL_VERIFY( edge_ );
        vertex_ = edge_->vertex( vi );
    }


    //@name Access status
    //@{

    /// Return incident face
    Edge * getEdge( ) { return edge_; }
    
    /// Return incident vertex
    Vertex * getVertex( ) { return vertex_; }
    //@}


    //@name Walking steps
    //@{

    /// Walk one edge ahead. This updates the walker position.
    /// \return    A walker starting at new position.
    Walker< EDGE > next( )
    {
        const int localIndex = edge_->localVertexIndex( vertex_ );
        edge_ = edge_->neighbor( localIndex );
        const int localIndexN = edge_->localVertexIndex( vertex_ );
        const int nextLocalIndex = ShapeProp< myShape >::next[ localIndex ];
        vertex_ = edge_->vertex( nextLocalIndex );
        return Walker< EDGE >( edge_, vertex_ );

    }

    /// Walk one edge back. This updates the walker position.
    /// \return    A walker starting at new position.
    Walker< EDGE > previous( )
    {
        const int localIndex = edge_->localVertexIndex( vertex_ );
        const int prevLocalIndex = ShapeProp< myShape >::previous[ localIndex ];
        vertex_ = edge_->vertex( prevLocalIndex );
        edge_ = edge_->neighbor( prevLocalIndex );
        return Walker< EDGE >( edge_, vertex_ );
    }

    /// Jump to other side of edge. This updates the walker position.
    /// \return    A walker starting at new position.
    Walker< EDGE > opposite( )
    {
        const int localIndex = edge_->localVertexIndex( vertex_ );        
        const int nextLocalIndex = ShapeProp< myShape >::next[ localIndex ];
        vertex_ = edge_->vertex( nextLocalIndex );
        return Walker< EDGE >( edge_, vertex_ );        
    }

    /// Move the walker in a combination of directions
    void move( const unsigned n, const subdiv::curve::move * dir )
    {
        for (unsigned i = 0; i < n; i++ ){
            if (!edge_)
                break;

            subdiv::curve::move d = dir[i];
            
            const int localIndex = edge_->localVertexIndex( vertex_ );        
            const int nextLocalIndex = ShapeProp< myShape >::next[ localIndex ];
            const int prevLocalIndex = ShapeProp< myShape >::previous[ localIndex ];

            if (d== subdiv::curve::NEXT){
                edge_ = edge_->neighbor( localIndex );
                const int localIndexN = edge_->localVertexIndex( vertex_ );
                const int nextLocalIndex = ShapeProp< myShape >::next[ localIndexN ];
                vertex_ = edge_->vertex( nextLocalIndex );

            } else if (d ==subdiv::curve::PREVIOUS){
                vertex_ = edge_->vertex( prevLocalIndex );
                edge_ = edge_->neighbor( prevLocalIndex );
            } else if (d ==subdiv::curve ::OPPOSITE)
                vertex_ = edge_->vertex( nextLocalIndex );
        }
        return;        
    }
    //@}

private:
    /// incident edge
    Edge *                  edge_;
    /// incident vertex
    Vertex *                vertex_;
};

#endif
