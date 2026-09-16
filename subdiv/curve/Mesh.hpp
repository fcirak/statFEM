//
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

#ifndef subdiv_curve_mesh_h
#define subdiv_curve_mesh_h

#include <vector>
#include <map>
#include <set>
#include <iostream>
#include <iomanip>
#include <utility>
#include <algorithm>
#include <cassert>
#include <cstdlib>
#include <algorithm>
#include <functional>
#include <iterator>
#include <array>

#include <corlib/verify.hpp>
#include <corlib/Shape.hpp>
#include <corlib/SmfHead.hpp>

#include <subdiv/curve/helpers.hpp>
#include <subdiv/curve/BuildChildren.hpp>
#include <subdiv/curve/ShapeDim.hpp>
#include <subdiv/curve/ShapeDim.hpp>
//------------------------------------------------------------------------------
namespace subdiv{
    namespace curve{
        template <typename V, typename E, typename ET>
	class Mesh;	    

        // namespace detail_{
        //     std::string GetFileName(std::string smfFilename, std::string type, const unsigned n);
        // }
    }
}
// //--------------------------------------------------------------------------
// std::string subdiv::curve::detail_::GetFileName(std::string smfFilename, 
//                                           std::string type, const unsigned n)
// {
//     std::string outputFileName = smfFilename.substr( 0, smfFilename.find( ".smf" ) ) + 
//         type + "_" + boost::lexical_cast< std::string >( n ) + ".smf";
    
//     return outputFileName;
// }

//-------------------------------------------------------------------------- 
template <typename V, typename E, typename ET>
class subdiv::curve::Mesh  {
public:

    static const corlib::shape   myShape       = ET::myShape;
    static const unsigned        numVertices   = ET::numVertices; 
    static const unsigned        numNeighbors  = ET::numNeighbors; 
    
    static const unsigned        dim           = V::dim; 
    static const unsigned        localDim      = ET::localDim; 

    // typedefs
    typedef V                                               Vertex;
    typedef E                                               Edge;
    typedef ET                                              ETree;

    typedef typename ET::VecVPtrNV                          VecVPtrNV;
    typedef std::vector<Vertex* >                           VecVPtr; 
    typedef std::vector<ETree* >                            VecETreePtr; 

    typedef std::map<Edge*, ETree* >                        ETreeMap;
    
    // iterators
    typedef typename std::vector<Vertex* >::iterator        VertexIterator;
    typedef typename std::set<Vertex* >::iterator           VertexSetIterator;
    typedef typename std::vector<Edge* >::iterator          EdgeIterator;
    typedef typename std::vector<ETree* >::iterator         EdgeTreeVecIterator;
    typedef typename ETreeMap::iterator                     EdgeTreeIterator;    
    
    typedef std::vector< std::pair< std::string, unsigned > > NamedAddCoord;

public:
    Mesh(std::istream & smf);

    ~Mesh();

    //! @name Generic application of functors to containers
    //@{
    //! Iterate over vertices
    template<typename OP>  OP iterateOverVertices( OP op ); 
    //! Iterate over edge
    template<typename OP>  OP iterateOverEdges( OP op ); 
    //! Iterate over edgeTrees
    template<typename OP>  OP iterateOverEdgeTrees( OP op ); 

    VertexIterator vBegin( ) { return vertices_.begin( ); }
    VertexIterator vEnd( ) { return vertices_.end( ); }
    EdgeIterator eBegin( ) { return edges_.begin( ); }
    EdgeIterator eEnd( ) { return edges_.end( ); }
    EdgeTreeIterator etBegin( ) { return edgeTrees_.begin( ); }
    EdgeTreeIterator etEnd( ) { return edgeTrees_.end( ); }
    //@}
    
    /// Get edge tree
    ETree * edgeTree( Edge * edge )
    {
        return edgeTrees_[ edge ];
    }

    //! Read edge Tags
    void readTags( std::istream & tg );

    //! Build the tree structure
    void buildTopology();

    //! Subdivide mesh
    template< typename SUBDIV >
    void subdivide( SUBDIV * subdiv ); 

    /// Write mesh to output files setting indices in vertices and edges/leafs
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
    /// \param[in]      addCoord  Names of additional vertex co-ordinates to write to vdat
    /// \param[out]     vdat      VDAT file (skipped if NULL)
    void writeMesh(std::ostream& out,
                   std::ostream * tg = NULL,
                   const NamedAddCoord addCoord = NamedAddCoord(),
                   std::ostream * vdat = NULL );

private:
    /// Collect all vertices
    void collectETreesVertices_(VecVPtr & vVec, VecETreePtr & eVec);

    //! Update facet tree and vertex indices
    void updateIndices_();

    //! Add vertex
    void addVertex_(Vertex* v) {vertices_.push_back(v);}
    
    //! Add edge
    void addEdge_(Edge* e) {edges_.push_back(e);}

protected:
    int                        refinementLevel_;

private:
    static const std::array<int,ShapeDim<ET::myShape>::numVertices>       next_;
    static const std::array<int,ShapeDim<ET::myShape>::numVertices>       previous_;

    std::vector<Vertex* >      vertices_;        // vertices
    std::vector<Edge* >        edges_;           // edges	
    ETreeMap                   edgeTrees_;       // edge tree map

};
#include <subdiv/curve/Mesh.ipp>

#endif

