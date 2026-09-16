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

#ifndef subdiv_surf_edge_h
#define subdiv_surf_edge_h

#include <iostream>
#include <vector>
#include <utility>
#include <functional>
#include <string>

#include <subdiv/surf/MeshTags.hpp>
//------------------------------------------------------------------------------
namespace subdiv{
    namespace surf{

        template< typename V, typename F >
        class Edge;

    }
}


//------------------------------------------------------------------------------
/// Edge object. It stores two pointers to its start #v1_ and end vertex #v2_. It
/// also stores pointers to the facets #facets_ which share the edge.
/// Each facet in #facets_ can have a manifold neighbour. The manifold neighbour
/// is linked to the facet to establish a "smooth crease" from the facet to
/// its neighbour.
///
///\verbatim
///
///       v1--------------v2
///             [this]
///             
///\endverbatim
///
/// \tparam   V     Vertex type
/// \tparam   F     Facet type
template< typename V, typename F >
class subdiv::surf::Edge
{
public:
    typedef V                                            VertexType;
    typedef F                                            FacetType;
    typedef Edge< V, F >                                 EdgeType;

    typedef enum subdiv::surf::edgeTag                   EdgeTag;
    typedef enum subdiv::surf::edgeCategory              EdgeCategory;

    /// An edge-related facet (first) and its manifold neighbour (second).
    /// If the facet does not have a manifold neighbour, the second
    /// entry is NULL.
    typedef std::pair< FacetType *, FacetType * >        FacetPair;

private:
    typedef std::vector< FacetPair >                     FacetPairVec_;

public:
    typedef typename FacetPairVec_::iterator             FacetPairIter;

public:
    /// Constructor
    ///
    /// \param[in]    v1     Start vertex
    /// \param[in]    v2     End vertex
    Edge( VertexType * v1, VertexType * v2 );

    /// Destructor
    ~Edge( ) { }

    /// Add an facet and its tag, if not already. Does _not_ overwrite
    /// neighbouring facet, if nor explicitly specified
    ///
    /// \param[in]   eTag        Edge tag
    /// \param[in]   f           Facet pointer to add
    /// \param[in]   fm          Its manifold neighbour, or NULL if none
    /// \param[in]   overwrite   Overwrite manifold neighbour facet, if true
    void insertFacet( FacetType * f, FacetType * fm,
                      bool overwrite = false );

    /// Number of edge forms
    unsigned numFacets( ) const { return facets_.size( ); }
    FacetPairIter facetsBegin( ) { return facets_.begin( ); }
    FacetPairIter facetsEnd( ) { return facets_.end( ); }

    /// Return pointer to start vertex
    VertexType * first( ) const { return v1_; }
    /// Return pointer to end vertex
    VertexType * second( ) const { return v2_; }

    /// Get pointer to neighbouring facet
    /// (With i=0 its the edge's master facet)
    FacetType * getFacet( const unsigned i = 0 ) const { return facets_[ i ].first; }
    /// Set pointer to a neighbouring facet
    void setFacet( FacetType * f, FacetType * fm )
    {
        facets_[ 0 ].first = f; 
        facets_[ 0 ].second = fm;
        return;
    }

    /// Less-than comparison functor for two edges
    ///
    /// \param[in]    e1   First edge
    /// \param[in]    e2   Second edge
    /// \return            True if first edge is "less" than second edge
    struct Less :
        public std::binary_function< const EdgeType *, const EdgeType *, bool >
    {
        /// Constructor
        Less( ) { }
 
       /// Main function
        bool operator()( const EdgeType * e1, const EdgeType * e2 ) const 
        {
        	 const VertexType * m1 = std::min( e1->first(), e1->second() );
        	        const VertexType * m2 = std::min( e2->first(), e2->second() );
        	        return ( ( m1 < m2 ) or
        	                 ( ( m1 == m2 ) and
        	                   ( std::max( e1->first(), e1->second() ) <
        	                     std::max( e2->first(), e2->second() ) ) ) );

            //return EdgeType::less( e1, e2 );
        }
    };

    /// Give tag
    EdgeTag getTag( ) const { return tag_; }
    /// Set tag
    void setTag( const EdgeTag eTag ) { tag_ = eTag; return; }

    /// Give category
    EdgeCategory getCategory( ) const { return cat_; }
    /// Set category
    void setCategory( const EdgeCategory eCat ) { cat_ = eCat; return; }

    /// Write something about my vertices, tag dreams, etc.
    std::ostream & write( std::ostream & os ) const;

private:
    VertexType *    v1_;         ///< First vertex
    VertexType *    v2_;         ///< Second vertex
    EdgeCategory    cat_;        ///< Edge category
    EdgeTag         tag_;        ///< Edge tag
    FacetPairVec_   facets_;     ///< Facets attached to edge
};

#include "Edge.ipp"

#endif 
