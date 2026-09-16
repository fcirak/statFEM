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

#ifndef subdiv_curve_edge_h
#define subdiv_curve_edge_h

#include <cstdlib>
#include <array>

#include <corlib/verify.hpp>
#include <subdiv/curve/ShapeDim.hpp>

namespace subdiv{
    namespace curve{
        template <corlib::shape SHAPE, typename V>
        class Edge;
    }
}
        

template <corlib::shape SHAPE, typename V>
class subdiv::curve::Edge {    
public:
    static const corlib::shape myShape           = SHAPE;
    static const unsigned      numVertices       = subdiv::curve::ShapeDim<myShape>::numVertices; 
    static const unsigned      numNeighbors      = subdiv::curve::ShapeDim<myShape>::numNeighbors;

private:
    typedef V                                    Vertex_;
    typedef subdiv::curve::Edge<myShape, V>      Edge_;

    typedef std::array<Vertex_*, numVertices>    VecVPtrNV_;
    typedef std::array<Edge*, numNeighbors>      VecEPtrNN_;

public:
    Edge(){}

    Edge(VecVPtrNV_ vec): vertices_(vec){
        for (unsigned i = 0; i < numNeighbors; i++)
            neighbors_[i]=NULL;
    } 

    ~Edge(){}
    
    /** @name Interface */
    //@{
    Vertex_* vertex(const unsigned i ){return vertices_[i];}
    inline Edge<SHAPE, V>*& opposite( const Vertex_*);
    inline Edge<SHAPE, V>*& neighbor(int i){return neighbors_[i];}
    //@}
    
    /** @name Output */
    //@{
    void printConnectivityData(std::ostream& os) {
        for (unsigned i = 0; i < numVertices; i++)
            os << vertices_[i]->index() << " ";
    }
    
    //! Return Index
    unsigned index() const {return index_;} 
    void setIndex(const unsigned i) {index_ = i;} 
    //@}

private:
    static const std::array<int,ShapeDim<SHAPE>::numVertices>       next_;
    static const std::array<int,ShapeDim<SHAPE>::numVertices>       previous_;
    unsigned          index_;            // vertex index       

    VecVPtrNV_        vertices_;         // vertices
    VecEPtrNN_        neighbors_;        // neighboring edges
};

#include <subdiv/curve/Edge.ipp>
    
#endif
