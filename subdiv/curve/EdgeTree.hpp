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

#ifndef subdiv_curve_edge_tree_h
#define subdiv_curve_edge_tree_h

#include <cstdlib>
#include <iterator>
#include <utility>
#include <set>
#include <algorithm>
#include <array>

#include <corlib/misc.hpp>
#include <corlib/verify.hpp>
#include <corlib/linalg.hpp>

#include <subdiv/curve/ShapeDim.hpp>
#include <subdiv/curve/MeshTags.hpp>
#include <subdiv/curve/CompareFunctor.hpp>
#include <subdiv/curve/Subdivision.hpp>
#include <subdiv/curve/BuildChildren.hpp>
#include <subdiv/curve/helpers.hpp>

namespace subdiv{
    namespace curve{
        template < corlib::shape SHAPE, typename V>
        class EdgeTree;
    }
}
        
template < corlib::shape SHAPE, typename V>
class subdiv::curve::EdgeTree {
public:

    static const corlib::shape myShape      = SHAPE;
    static const unsigned      numVertices  = subdiv::curve::ShapeDim<myShape>::numVertices; 
    static const unsigned      numNeighbors = subdiv::curve::ShapeDim<myShape>::numNeighbors;
    static const unsigned      numChildren  = subdiv::curve::ShapeDim<myShape>::numChildren; 

    static const unsigned      dim          = V::dim; 
    static const unsigned      localDim     = subdiv::curve::ShapeDim<myShape>::dim; 
    static const unsigned      numDeriv2    = localDim + subdiv::curve::Factorial<localDim>::result;

    typedef V                                               Vertex;

    typedef subdiv::curve::EdgeTree<SHAPE, V>     ETree;
    typedef subdiv::curve::refinementStatus       RefinementStatus;

    typedef std::array<Vertex*, numVertices>      VecVPtrNV;
    typedef std::array<ETree*, numChildren>       VecETPtrNC;
    typedef std::array<ETree*, numNeighbors>      VecETPtrNN;

    typedef std::vector<Vertex*>                  VecVPtr;
    typedef std::vector<ETree*>                   VecETreePtr;
    typedef typename VecETPtrNC::iterator         VecETPtrNCItr;
    typedef typename VecVPtrNV::iterator          VecVPtrNVItr;

public:	
    EdgeTree(VecVPtrNV vec, int depth, ETree* parent); 

    ~EdgeTree();
    
    /** @name Interface */
    //@{
    /// Get local vertex index by vertex pointer
    ///
    /// \param[in]   v   Pointer to vertex
    /// \return          Local index of vertex
    unsigned localVertexIndex( const Vertex * v ) const;

    /// Get vertex
    Vertex*& vertex(const int i){return vertices_[i];}

    /// Get neighbor
    ETree*& neighbor(const int i){return neighbors_[i];}

    /// Get parent
    ETree*& parent(){return parent_;}

    /// Get child
    ETree*& child(int i){return children_[i];}

    /// Get neighbor opposite to vertex v
    EdgeTree<SHAPE, V>*& opposite( const Vertex* v);

    /// Return iterator to beginning of vertices
    VecVPtrNVItr vBegin( ) { return vertices_.begin( ); }

    /// Return iterator to end of vertices
    VecVPtrNVItr vEnd( ) { return vertices_.end( ); }

    /// Get depth of refinement
    int depth()const {return depth_;}
    
    /// Return Index
    unsigned index() const {return index_;} 
    
    /// Set index
    void setIndex(const unsigned i) {index_ = i;} 
        
    /// Check if is a leaf (has no children)
    bool isLeaf();
    //@}
    
    /** @name Operations */
    //@{
    /// Subdivide this edge tree
    template< typename SUBDIV >
    void subdivide(const int level, SUBDIV * subdiv);

    /// Collect leafs to a vector
    template <typename  OutputIterator>
    void addLeafs(OutputIterator iter);

    /// Collect all leafs to a vector
    template <typename  OutputIterator>
    void addLeafsAll(OutputIterator iter);

    /// Collect all new vertices
    template< typename OUTITER >
    void addNewVertices( OUTITER vtxs );

    /// Collect all new vertices up to a level
    template< typename OUTITER >
    void addNewVertices(const int level, OUTITER vtxs );

    /// Find extordinary vertex if any
    Vertex* extVertex();
    //@}
    
private:
    /// Refine the edge vertex
    template< typename SUBDIV >
    void newCoordEdgeVertices_(Vertex*& v,  SUBDIV * subdiv );
    
    /// Refine the existing vertices
    template< typename SUBDIV >
    void newCoordExistingVertices_( SUBDIV * subdiv );

    /// Clear new vertices
    void clearNewVertices_( ) { newVertices_.clear( ); }
    
    // /// Reverse Subedivide  all vertices
    // void unrefineVertices_( const int level );
    
protected:
    /// Subdivide edge tree leaf
    template< typename SUBDIV >
    void subdivideLeaf_( SUBDIV * subdiv );
    
    /// Return iterator to beginning of children
    VecETPtrNCItr cBegin( ) { return children_.begin( ); }

    /// Return iterator to end of children
    VecETPtrNCItr cEnd( ) { return children_.end( ); }

    /// Clear new vertices
    void clearNewVertices( ) { newVertices_.clear( ); }

    /// Clear children
    void clearChildren( ) { children_.clear( ); }

private:
    static const std::array<int, subdiv::curve::ShapeDim<SHAPE>::numVertices > next_;
    static const std::array<int, subdiv::curve::ShapeDim<SHAPE>::numVertices > previous_;

    VecVPtrNV                  vertices_ ;           // vertices
    ETree*                     parent_;              // parent pointer 
    VecETPtrNC                 children_;            // child pointers     
    VecETPtrNN                 neighbors_;           // neighbors
    unsigned                   index_;               // index
    int                        depth_;               // depth of refinement

protected:
    VecVPtr                    newVertices_;         // newly allocated vertices used by children

}; 
    
#include <subdiv/curve/EdgeTree.ipp>

#endif
