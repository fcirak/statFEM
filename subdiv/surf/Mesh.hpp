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

#ifndef subdiv_surf_mesh_h
#define subdiv_surf_mesh_h

#include <vector>
#include <map>
#include <set>
#include <iostream>
#include <iomanip>
#include <utility>
#include <algorithm>

#include <corlib/Shape.hpp>
#include <corlib/misc.hpp>

#include <subdiv/surf/ShapeDim.hpp>
#include <subdiv/surf/FacetTree.hpp>

//------------------------------------------------------------------------------
namespace subdiv{
    namespace surf{

        template< typename V, typename E, typename F, typename FT >
        class Mesh;

    }
}

//------------------------------------------------------------------------------
/// A surface mesh made of polygonal facets and their vertices.
///
/// The surface mesh, which can be non-manifold, is constructed by
/// reading the control vertices and facets. The facets are copied
/// to facet tress (FTree) which can refine itself. The refinement level
/// can differ by one level from one facet tree to its neighbour.
/// The refinement of the facet tress generates new vertices, which are
/// kept in the facet trees. The facets are either triangles or
/// quadrilaterals.
///
/// \tparam  V     Vertex type
/// \tparam  E     Edge type
/// \tparam  F     Facet type
/// \tparam  FT    Facet tree type
template< typename V, typename E, typename F, typename FT >
class subdiv::surf::Mesh
{
public:
    static const corlib::shape   myShape       = FT::myShape;
    static const unsigned        numVertices   = FT::numVertices;
    static const unsigned        numNeighbors  = FT::numNeighbors;
    static const unsigned        numEdges      = FT::numEdges;

    // typedefs
    typedef V                                               Vertex;
    typedef E                                               Edge;
    typedef F                                               Facet;
    typedef FT                                              FTree;

    typedef typename Edge::Less                             EdgeLess;
    typedef std::set< Edge *, EdgeLess >                    EdgeSet;
    
    typedef typename FTree::VecVPtrNV                       VecVPtrNV;
    typedef std::map< Facet *, FTree * >                    FTreeMap;
    
    // iterators
    typedef typename std::vector< Vertex * >::iterator      VertexIterator;
    typedef typename EdgeSet::iterator                      EdgeIterator;
    typedef typename std::vector< Facet * >::iterator       FacetIterator;
    typedef typename FTreeMap::iterator                     FacetTreeIterator;
    typedef typename EdgeSet::const_iterator                EdgeSetConstIterator;

    typedef std::vector< std::pair< std::string, unsigned > > NamedAddCoord;
    
public:
    /// Default constructor reading the mesh
    ///
    /// \param[in]   smf   SMF mesh file stream
    Mesh( std::istream & smf );

    /// Empty constructor --- you have to read mesh later
    Mesh( ) : verbose_( false ) { }

    /// Destructor
    ~Mesh( );

    /// @name Generic application of functors to containers
    //@{
    /// Iterate over vertices
    template< typename OP >
    OP iterateOverVertices( OP op );
    /// Iterate over edges
    template< typename OP >
    OP iterateOverEdges( OP op );
    /// Iterate over facets
    template< typename OP >
    OP iterateOverFacets( OP op );
    /// Iterate over facetTrees
    template< typename OP >
    OP iterateOverFacetTrees( OP op );

    VertexIterator vBegin( ) { return vertices_.begin( ); }
    VertexIterator vEnd( ) { return vertices_.end( ); }
    EdgeIterator eBegin( ) { return edges_.begin( ); }
    EdgeIterator eEnd( ) { return edges_.end( ); }
    FacetIterator fBegin( ) { return facets_.begin( ); }
    FacetIterator fEnd( ) { return facets_.end( ); }
    FacetTreeIterator ftBegin( ) { return facetTrees_.begin( ); }
    FacetTreeIterator ftEnd( ) { return facetTrees_.end( ); }
    //@}

    /// @name Set-up methods
    //@{

    /// Read SMF file
    ///
    /// SMF file format
    ///\verbatim
    ///     ! header line
    ///     ! ...
    ///     <numVertices>  <numFacets>
    ///     <vertexCoord0> <vertexCoord1> <vertexCoord2>
    ///     <vertexCoord0> <vertexCoord1> <vertexCoord2>
    ///     ...            ...            ...
    ///     <vertexId> <vertexId> <vertexId> [ <vertexId> ]
    ///     <vertexId> <vertexId> <vertexId> [ <vertexId> ]
    ///     ...        ...        ...        ...
    ///\endverbatim
    ///
    /// \param[in]   smf    SMF mesh file stream
    void readSmf( std::istream & smf );

    /// Read tag file
    ///
    /// Tag file format:
    ///\verbatim
    ///     TG
    ///     <numVertexTags>  <numEdgeTags>
    ///     <vertexId>  <vertexTag>
    ///     <vertexId>  <vertexTag>
    ///     ...        ...
    ///     <vertexId> <vertexId>  <edgeTag>  [ <facetId> <facetId> ]
    ///     <vertexId> <vertexId>  <edgeTag>  [ <facetId> <facetId> ]
    ///     ...        ...
    ///\endverbatim
    ///
    /// Hint:
    /// - This method can _only_ be called _before_ #buildFacetTopology().
    /// - The vertexTag and edgeTag are integers which are converted to
    ///   enum #subdiv::surf::vertexTag or enum #subdiv::surf::edgeTag.
    /// - The optional \tt <facetId> appears only for #subdiv::surf::EDGE_SMOOTHCREASE
    ///   edges and indicates which facets (i.e. their manifolds) are considered
    ///   like a ordinary (i.e. "smooth") spline interpolation across
    ///   the edge. The edge is a border of both manifolds/facets.
    ///
    /// \param[in]  tg      Tag file stream
    void readTags( std::istream & tg );
   
    /// Build the neighbourhood topology, align orientation of (true)  manifold
    void buildFacetTopology( );

    /// Build facet trees
    ///
    /// This method generates a facet tree for each facet object.
    /// The "adjoint" facet tree is stored in the map #facetTrees_.
    void implantFacetTrees( );

    /// Get facet tree
    ///
    /// \param[in]  facet   A pointer to a facet in the mesh
    /// \return             Pointer to adjoint facet tree
    FTree * facetTree( Facet * facet )
    {
        return facetTrees_[ facet ];
    }

    //@}

    /// Subdivide mesh
    ///
    /// \param[in]  subdiv   Subdivision scheme
    template< typename SUBDIV >
    void subdivide( SUBDIV * subdiv );

    /// Subdivide a facet tree
    ///
    /// \param[in]  facet   Facet whose facet tree is refined
    /// \param[in]  subdiv  Subdivision scheme
    template< typename SUBDIV >
    void subdivideFacetTree( Facet * facet, SUBDIV * subdiv )
    {
        facetTrees_[ facet ]->subdivide( subdiv );
        return;
    }

    /// Set verbosity
    void setVerbose( ) { verbose_ = true; return; }

    /// Unset verbosity
    void unsetVerbose( ) { verbose_ = true; return; }

    /// Write mesh to output files setting indices in vertices and facets/leafs
    /// 
    /// This method is considerably faster than #writeMesh(), especially
    /// for large meshes.
    ///
    /// Format of VDAT stream
    ///\verbatim
    ///     <numVertices> <numData>
    ///     <nameFirstData>  <numFirstDataPerVertex>
    ///     <nameSecondData> <numSecondDataPerVertex>
    ///     ...              ...
    ///     <vertexFirstDat0> <vertexFirstDat1> <vertexFirstDat2> ...
    ///     <vertexFirstDat0> <vertexFirstDat1> <vertexFirstDat2> ...
    ///     ...               ...               ...               ...
    ///     <vertexSecondDat0> <vertexSecondDat1> ...
    ///     <vertexSecondDat0> <vertexSecondDat1> ...
    ///     ...                ...                ...
    ///\endverbatim
    ///
    /// \param[out]     smf       SMF file mesh stream
    /// \param[out]     tg        TG file tag stream (skipped if NULL)
    /// \param[in]      addCoord  Names and number of components of additional
    ///							  vertex data to write to vdat
    /// \param[out]     vdat      VDAT file (skipped if NULL)
    void writeMeshWithIndex( std::ostream & smf, std::ostream * tg = NULL,
                             const NamedAddCoord addCoord = NamedAddCoord(),
                             std::ostream * vdat = NULL );

    /// Temp debugging job
    void tempJob( );

private:
    /// Add vertex
    ///
    /// \param[in]   v   Vertex to add
    void addVertex_( Vertex * v ) { vertices_.push_back( v ); }

    /// Add facet
    ///
    /// \param[in]   f   Facet to add
    void addFacet_( Facet * f ) { facets_.push_back( f ); }

    /// Update neighbor relations of the facets connected to the edge
    /// determined by its first and second vertices. Store edge in #edges_
    /// if treated first time.
    ///
    /// \param[in]      v1     First vertex on edge
    /// \param[in]      v2     Second vertex on edge
    /// \param[in,out]  f      Facet (is linked to opposite face, if edge already in #edges_)
    void addOrMatch_( Vertex * v1, Vertex * v2, Facet * f );

    /// Update facet tree and vertex indices
    ///
    /// \param[out]  numVerticesMesh  Number of all vertices including new ones due to refinement
    /// \param[out]  numLeafs         Number of all leafs, i.e. refined facets
    void updateIndices_( unsigned & numVerticesMesh, unsigned & numLeafs );

    /// Make sure every manifold is consistently oriented
    ///
    /// This method acts on the #facets_ and relies on an existing
    /// neighborhood relations build by #addOrMatch_()
    void orientateManifolds_( );

private:
    /// permutators to find next node on cell
    static const std::array< int, corlib::ShapeTraits< FT::myShape >::numVertices >  next_;
    /// permutators to find previous node on cell
    static const std::array< int, corlib::ShapeTraits< FT::myShape >::numVertices >  previous_;

    std::vector< Vertex * >    vertices_;        ///< vertices
    std::vector< Facet * >     facets_;          ///< (control) facets
    FTreeMap                   facetTrees_;      ///< facet tree map
    EdgeSet                    edges_;           ///< edges

    bool                       verbose_;         ///< verbosity status
};

//------------------------------------------------------------------------------
#include "Mesh.ipp"

#endif

