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

#ifndef subdiv_curve_vertex_h
#define subdiv_curve_vertex_h

#include <map>
#include <set>
#include <utility>
#include <algorithm>
#include <iostream>

#include <Eigen/Core>
#include <corlib/eigenX.hpp>

#include <subdiv/curve/MeshTags.hpp>
#include <subdiv/curve/CompareFunctor.hpp>
#include <subdiv/curve/VertexProp.hpp>
#include <subdiv/curve/helpers.hpp>

namespace subdiv{
    namespace curve{
        template<unsigned N>
	    class Vertex; 	    
        
	    namespace eigenX = corlib::eigenX;
    }
}

template<unsigned N>
class subdiv::curve::Vertex {

public:
    static const unsigned     dim              = 3; 
    //to speed up computation of subdivision matrix
    static const unsigned     numVariables     = 3 + N;

    typedef subdiv::curve::vertexTag           VertexTag;
    typedef eigenX::VectorSd<numVariables>     VecP;

    typedef subdiv::curve::VertexProp<Vertex>                VertexProp;
    typedef std::map<int, VertexProp, std::greater<int> >    VertexPropMap; 
    typedef typename VertexPropMap::iterator                 LevelIterator;

public:	
    Vertex(int level=0);
    
    ~Vertex();    
    
    /// iterators
    LevelIterator levBegin() {return vertexPropMap_.begin();}
    LevelIterator levEnd() {return vertexPropMap_.end();}
    
    /// Index
    unsigned index() const {return index_;} 
    void setIndex(const unsigned i) {index_ = i;} 
    
    /// Vertex tag
    void setTag(VertexTag tg, int level = 0){ vertexPropMap_[level].tag = tg; }
    VertexTag tag(const int level){ return vertexPropMap_[level].tag; } 

    /// Coordinate points
    VecP point(){ return vertexPropMap_.begin()->second.point; } 
    VecP point(const int level){return vertexPropMap_[level].point;} 
    void setPoint(VecP v, const int level=0){ vertexPropMap_[level].point = v;}

    // active level
    int maxActiveLevel(){return vertexPropMap_.begin()->first;}
    int minActiveLevel()
    {
        LevelIterator lEnd = vertexPropMap_.end();
        lEnd--;
        return lEnd->first;
    }

    // /// Last unrefinement level
    // int unRefLevel( ) { return unrefLevel_; }
    // void setUnRefLevel( const int level ) { unrefLevel_ = level; }

    //! Add new refinement level
    void addNewLevel(const VecP vec, int level);
        
    /// Remove level
    void removeLevel(const int level)
    {
        LevelIterator itr = vertexPropMap_.find( level );
        assert( itr != vertexPropMap_.end( ) );
        vertexPropMap_.erase(itr);
        return;
    }

    // /// Reset vertex at particular level
    // void reset(const int level = 0)
    // {
    //     //std::cout<<index_<<" "<<this->maxActiveLevel()<<" |";//TEMP
    //     for (int l = level + 1; l < this->maxActiveLevel() + 1; l++)
    //         this->removeLevel(l);
    //     //std::cout<<this->maxActiveLevel()<<std::endl;//TEMP
    //     unrefLevel_ = std::numeric_limits<int>::max();
    //     return;
    // }

    // io
    friend void operator>>(std::istream& is, 
                           subdiv::curve::Vertex<N>* v) {
        VecP coor;
        is >> coor[0] >> coor[1]  >> coor[2];
        v->setPoint(coor);
    }

    friend std::ostream& operator<<(std::ostream& os, 
                                    subdiv::curve::Vertex<N>* v) {
        VecP coor = v->point();
        os << coor[0] << " "<< coor[1] << " "<< coor[2];
        return os;
    }

private:
    unsigned          index_;        // vertex index       
    //int               unrefLevel_;
    VertexPropMap     vertexPropMap_;
};
    
#include <subdiv/curve/Vertex.ipp>

#endif 
