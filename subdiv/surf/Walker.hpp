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

#ifndef subdiv_surf_walker_h
#define subdiv_surf_walker_h

#include <cstdlib>
#include <iterator>
#include <functional>
#include <algorithm>
#include <array>

#include <corlib/verify.hpp>
#include <subdiv/surf/ShapeDim.hpp>


//------------------------------------------------------------------------------
namespace subdiv{
    namespace surf{

        //! \brief Enumerator for walker movement
        enum move {
            NEXT,               //!< next
            PREVIOUS,           //!< previous
            OPPOSITE            //!< opposite
        };

        template< typename FACET >
        class Walker;
    }
}

//==============================================================================
/// Walker navigates over the mesh.
///
/// It stores a facet and one of its vertices and updates this pair
/// due to #next(), #previous() and #opposite() steps.
///
/// These operations rely on <i>consistently oriented</i> faces.
///
/// <i>Walk, don't run.</i>
///
/// \tparam  FACET         Facet (tree) type
template< typename FACET >
class subdiv::surf::Walker
{
public:
    typedef FACET                                           Facet;
    typedef typename Facet::Vertex                          Vertex;

    static const corlib::shape myShape      = Facet::myShape;

public:
    /// Constructor
    ///
    /// \param[in]    f         Pointer to facet
    /// \param[in]    v         Pointer to vertex of facet
    Walker( Facet * f, Vertex * v );

    /// Constructor with local vertex index
    ///
    /// \param[in]    f         Pointer to facet
    /// \param[in]    vi        Facet-wise vertex index
    Walker( Facet * f, const int vi );

    /// Destructor
    ~Walker( );

    /// Set new position with facet-wise vertex index
    ///
    /// \param[in]    f         Pointer to facet
    /// \param[in]    vi        Facet-wise vertex index
    void set( Facet * f, const int vi );

    //@name Access status
    //@{

    /// Return incident face
    Facet * getFacet( ) { return facet_; }
    
    /// Return incident vertex
    Vertex * getVertex( ) { return vertex_; }

    //@}

    //@name Walking steps
    //@{

    /// Walk one edge ahead. This updates the walker position.
    ///
    /// Sketch with [f] = incident face, [v] = incident vertex
    ///\verbatim
    ///                                       * 
    ///     |       |                 |       ^
    ///     |  [f]  |       next      |  [f]  |
    ///     |       |       ===>      |       |
    ///   -[v]----->*--             --*------[v]-
    ///     |       |                 |       |
    ///
    ///\endverbatim
    ///
    /// \return    A walker starting at new position.
    Walker< FACET > next( );

    /// Walk one edge back. This updates the walker position.
    ///
    /// Sketch with [f] = incident face, [v] = incident vertex
    ///\verbatim
    ///                              [v]    
    ///     |       |                 |       |
    ///     |  [f]  |     previous    |  [f]  |
    ///     |       |       ===>      V       |
    ///   -[v]----->*--             --*-------*--
    ///     |       |                 |       |
    ///
    ///\endverbatim
    ///
    /// \return    A walker starting at new position.
    Walker< FACET > previous( );

    /// Jump to other side of edge. This updates the walker position.
    ///
    /// Sketch with [f] = incident face, [v] = incident vertex
    ///\verbatim
    ///                                     
    ///     |  [f]  |     opposite    |       |
    ///     |       |       ===>      |       |
    ///   -[v]----->*--             --*<-----[v]-
    ///     |       |                 |       |
    ///     |       |                 |  [f]  |
    ///
    ///\endverbatim
    ///
    /// \return    A walker starting at new position.
    Walker< FACET > opposite( );


    /// Move the walker in a combination of directions
    ///
    /// \tparam       N           Number of movement directions
    /// \param[in]    dir         Sequence of movement directions    
    template< unsigned N >
    void move( const std::array< subdiv::surf::move, N > & dir );
    //@}

protected:
    /// incident facet
    Facet *                   facet_;
    /// incident vertex
    Vertex *                  vertex_;
};

#include "Walker.ipp"

#endif
