using namespace subdiv::curve;

//! Constructor
template<unsigned N>
Vertex<N>::Vertex(int level)//:
//    unrefLevel_(std::numeric_limits< int >::max())
{
    VertexProp vp(level);
    vertexPropMap_[level] = vp;
    
    VecP newpoint;
    newpoint.setZero();
    vertexPropMap_[level].point = newpoint;
    vertexPropMap_[level].tag = subdiv::curve::NOTAG_VERTEX;
    index_ = 0;
}

template<unsigned N>
Vertex<N>::~Vertex()
{

}


//! Add new refinement level
template<unsigned N>
void Vertex<N>::addNewLevel(const VecP vec, int level)
{
    LevelIterator it = vertexPropMap_.find(level);	    
    if (it != vertexPropMap_.end())
        (*it).second.point = vec;
    else {
        VertexProp vp(level);
        vp.point = vec;
        vertexPropMap_.insert(typename VertexPropMap::value_type(level, vp));  
    }
    return;
}	


