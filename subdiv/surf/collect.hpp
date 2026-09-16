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

#ifndef subdiv_surf_collect_h
#define subdiv_surf_collect_h

#include <vector>
#include <set>
#include <map>
#include <functional>

#include <subdiv/surf/MeshTags.hpp>
#include <subdiv/surf/ShapeDim.hpp>
#include <subdiv/surf/Walker.hpp>
#include <subdiv/surf/Subdivision.hpp>

//------------------------------------------------------------------------------
namespace subdiv{
    namespace surf{
        
        //======================================================================
        /// Collect one-ring neighbourhood of facets and vertices of a facet
        /// on a manifold. The neighbourhood of the #facet must be a simple
        /// manifold but can be irregular.
        ///
        /// Example: The facets [0] to [7] are the one-neighbourhood
        ///          facets of [this] facet. And the vertices (*)
        ///          are the one-neighbourhood vertices of [this] facet.
        ///          The shown example one-ring is closed.
        ///\verbatim
        ///
        ///   |       |       |       |
        /// -(*)-----(*)-----(*)-----(*)-
        ///   |       |       |       |
        ///   |  [6]  |  [5]  |  [4]  |
        ///   |       |       |       |
        /// -(*)------+-------+------(*)-
        ///   |       |       |       |
        ///   |  [7]  | [this]|  [3]  |
        ///   |       |       |       |
        /// -(*)------+-------+------(*)-
        ///   |       |       |       |
        ///   |  [0]  |  [1]  |  [2]  |
        ///   |       |       |       |
        /// -(*)-----(*)-----(*)-----(*)-
        ///   |       |       |       |
        ///
        ///\endverbatim
        ///
        /// \tparam      FACET       Facet type
        ///
        /// \param[in]   facet       The centre facet
        /// \param[out]  oneRingF    One-ring of facets (counter-clockwise)
        /// \param[out]  oneRingV    One-ring of vertices (unordered)
        /// \param[out]  closedRing  True if one-ring is closed, otherwise false
        template< typename FACET >
        void collectManifoldOneRingFacets( FACET * facet,
                                           std::vector< FACET * > & oneRingF,
                                           std::set< typename FACET::Vertex * > & oneRingV,
                                           bool & closedRing );
        
        //======================================================================
        /// Collect one-ring neighbourhood of facets and vertices of a facet.
        /// The neighbour of the #facet can be non-manifold.
        ///
        /// \tparam      FACET       Facet type
        ///
        /// \param[in]   facet       The centre facet
        /// \param[out]  oneRingF    One-ring of facets (unordered)
        /// \param[out]  oneRingV    One-ring of vertices (unordered)
        template< typename FACET >
        void collectOneRingFacets( FACET * facet,
                                   std::vector< FACET * > & oneRingF,
                                   std::set< typename FACET::Vertex * > & oneRingV );

        //======================================================================
        /// Collect N-ring neighbourhood of facets and vertices of a facet
        ///
        /// \tparam      FACET       Facet type
        ///
        /// \param[in]   order       Order N>0
        /// \param[in]   facet       The centre facet
        /// \param[out]  oneRingF    N-ring of facets (unordered)
        /// \param[out]  oneRingV    N-ring of vertices (unordered)
        template< typename FACET >
        void collectNRingFacets( const unsigned order,
                                 FACET * facet,
                                 std::vector< FACET * > & nRingF,
                                 std::set< typename FACET::Vertex * > & nRingV );

        //======================================================================
        /// Collect one-ring of vertices around a vertex on a (multi-petal) manifold,
        /// i.e. the multi-petal manifold stretches across a non-manifold edge.
        /// The vertex lies on a <i>single</b> non-manifold edge across which 
        /// multi-petal manifold may be established.
        ///
        /// <h4>Note:</h4>
        /// - The facets can be differently oriented.
        /// - The vertex can be irregular.
        /// - The neighbourhoor can be open.
        /// - It works on ordinary single-petal manifolds too.
        ///
        /// Example: The 6 vertices indicated with <tt>(*)</tt> are the 1-neighbourhood
        ///          of regular vertex <tt>[i]</tt>
        ///\verbatim
        ///
        ///   `|  `|
        ///  -(*)-(*)-
        ///    |`  |`
        ///    | ` | `
        ///   `|  `|  `|
        ///  -(*)-[i]-(*)-
        ///    |`  |`  |`
        ///      ` | ` |
        ///       `|  `|
        ///      -(*)-(*)-
        ///        |`  |`
        ///
        ///\endverbatim
        ///
        /// Example: The 4 vertices indicated with <tt>(*)</tt> are the 1-neighbourhood
        ///          of a regular vertex <tt>[i]</tt>
        ///\verbatim
        ///
        ///    |     |     |
        ///   -+----(*)----+-
        ///    |     |     |
        ///  -(*)---[i]---(*)-
        ///    |     |     |
        ///   -+----(*)----+-
        ///    |     |     |
        ///
        ///\endverbatim
        ///
        /// Example: The 5 vertices indicated with <tt>(0)</tt> to <tt>(4)</tt> are the
        ///          open 1-neighbourhood of the vertex <tt>(i)</tt>. Vertex <tt>(0)</tt>
        ///          and <tt>(4)</tt> lie on the boundary of the manifold. They will
        ///          appears as first and last in #oneRingV. 
        ///\verbatim
        ///
        ///      ` /   ` /
        ///     -(3)---(2)-
        ///      / `   / `
        ///   ` /   ` /[f]` /
        ///  =(4)===(i)---(1)-
        ///          ``   / `
        ///           `` /
        ///            (0)-
        ///             ``
        ///      
        ///\endverbatim
        ///
        /// Example: The 6 vertices indicated with <tt>(0)</tt> to <tt>(5)</tt> are the
        ///          closed multi-petal 1-neighbourhood of the vertex <tt>(i)</tt>, if
        ///          the bold edges <tt>==</tt> and <tt>``</tt> indicate non-manifold edges.
        ///          The multi-petal connectivity across non-manifold edges is set
        ///          uni-laterally on petal containing facet <tt>[f]</tt>. 
        ///\verbatim
        ///
        ///      ` /   ` /
        ///     -(3)---(2)-
        ///      / `   / `
        ///   ` /   ` /[f]` /
        ///  =(4)===(i)---(1)-
        ///     `   /``   / `
        ///      ` /  `` /
        ///     -(5)---(0)-
        ///      / `   / `
        ///
        ///\endverbatim
        ///
        /// \param[in]   ft                Facet to which vertex is attached
        /// \param[in]   ivtx              Local index of vertex w.r.t. facet
        /// \param[in]   getNeighborFun    Facet's member function to get hold of its manifold neighbors
        /// \param[out]  oneRingV          One-ring of vertices. The one-ring is ordered, i.e. in case
        ///                                of open neighbourhood the first and last vertices lie on the
        ///                                boundary of the manifold.
        /// \param[out]  oneRingF          One-ring of facets, i.e. facets which share vertex and are
        ///                                on same manifold. The one-ring is ordered, i.e. in case 
        ///                                of open neighborhood the first and last facets touch the
        ///                                boundary of the manifold.
        /// \return                        True if closed neighbour
        template< typename FACET >
        bool collectManifoldOneRingVertices( FACET * ft,
                                             const int ivtx,
                                             FACET * & (FACET::*getNeighborFun)( const int ),
                                             std::vector< typename FACET::Vertex * > & oneRingV,
                                             std::vector< FACET * > & oneRingF,
                                             const bool stopOnCrease = false );

        //======================================================================
        /// Collect one-ring of vertices and facets around a vertex
        ///
        /// This method tolerates non-manifold topology.
        ///
        /// \tparam      FACET       Facet type
        ///
        /// \param[in]   ft          Facet to which vertex is attached
        /// \param[in]   ivtx        Local index of vertex w.r.t. facet
        /// \param[out]  oneRingV    One-ring of vertices
        /// \param[out]  oneRingF    One-ring of factes
        template< typename FACET >
        void collectOneRingVertices( FACET * ft,
                                     const int ivtx,
                                     std::map< typename FACET::Vertex *,
                                               std::pair< enum subdiv::surf::edgeTag,
                                                          enum subdiv::surf::edgeCategory > > & oneRingV,
                                     std::vector< FACET * > & oneRingF );

        //======================================================================
        /// Extract vertices which fall in given edgeTag group
        ///
        /// \tparam      VERTEX      Vertex type
        ///
        /// \param[in]   oneRingV    Vertices with kind-of-edge information
        /// \param[in]   eTag        Select by this given edeg tag
        /// \param[out]  vertices    The vertices in sub-group
        template< typename VERTEX >
        void extractByEdgeTag( const std::map< VERTEX *,
                                               std::pair< enum subdiv::surf::edgeTag,
                                                          enum subdiv::surf::edgeCategory > > & oneRingV,
                               const enum subdiv::surf::edgeTag eTag,
                               std::vector< VERTEX * > & vertices );

        //======================================================================
        /// Extract vertices which fall in given edgeCategory group
        ///
        /// \tparam      VERTEX      Vertex type
        ///
        /// \param[in]   oneRingV    Vertices with kind-of-edge information
        /// \param[in]   eCat        Select by this given edge category
        /// \param[out]  vertices    The vertices in sub-group
        template< typename VERTEX >
        void extractByEdgeCategory( const std::map< VERTEX *,
                                                    std::pair< enum subdiv::surf::edgeTag,
                                                               enum subdiv::surf::edgeCategory > > & oneRingV,
                                    const enum subdiv::surf::edgeCategory eCat,
                                    std::vector< VERTEX * > & vertices );

        //======================================================================
        /// Collect regular vertices around refined central facet on a manifold
        ///
        /// Please note: Vector of regular vertices has to follow same pattern
        ///              like the order of the related splines.
        ///
        /// Please note: Facets in manifold must have same orientation.
        ///
        /// \param[in]   centreFacet       Facet at centre of regular patch
        /// \param[out]  regularVertices   The surrounding regular vertices
        template< typename FACET, subdiv::surf::method METHOD >
        void collectManifoldRegularVertices(
            FACET * centreFacet,
            std::array< typename FACET::Vertex *,
                    Subdivision< METHOD, FACET >::numFunctions > & regularVertices
            );

        /// Collect regular vertices around refined central facet
        /// (specialisation for Catmull-Clark)
        ///
        /// Catmull-Clark:
        ///\verbatim
        ///    |       |       |       |  
        /// -(12)----(13)----(14)----(15)-
        ///    |       |       |       |  
        ///    |       |       |       |  
        ///    |       |       |       |  
        /// --(8)-----(9)----(10)----(11)-
        ///    |       |*******|       |  
        ///    |       |*******|       |  
        ///    |       |*******|       |  
        /// --(4)-----(5)-----(6)-----(7)-
        ///    |       |       |       |  
        ///    |       |       |       |  
        ///    |       |       |       |  
        /// --(0)-----(1)-----(2)-----(3)-
        ///    |       |       |       |
        ///\endverbatim
        template< typename FACET >
        void collectManifoldRegularVertices(
            FACET * centreFacet,
            std::array< typename FACET::Vertex *,
                                   Subdivision< subdiv::surf::CATMULL_CLARK,
                                                FACET >::numFunctions > & regularVertices
            );

        /// Collect regular vertices around refined central facet
        /// (specialisation for Loop)
        ///
        /// Loop:
        ///\verbatim
        ///
        ///     (8)-(11)      
        ///      |`  |`       
        ///      | ` | `      
        ///      |  `|  `     
        ///     (4)-(7)-(10)  
        ///      |`  |`  |`   
        ///      | ` |*` | `  
        ///      |  `|**`|  ` 
        ///     (1)-(3)-(6)-(9)
        ///       `  |`  |`  |
        ///        ` | ` | ` |
        ///         `|  `|  `|
        ///         (0)-(2)-(5)
        ///
        ///\endverbatim
        template< typename FACET >
        void collectManifoldRegularVertices(
            FACET * centreFacet,
            std::array< typename FACET::Vertex *,
                                   Subdivision< subdiv::surf::LOOP,
                                                FACET >::numFunctions > & regularVertices
            );

    }
}

//------------------------------------------------------------------------------
#include "collect.ipp"

#endif
        
