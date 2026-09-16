using namespace subdiv::curve;

//! Local index of next vertex
template <corlib::shape SHAPE, typename V>
const std::array<int,ShapeDim<SHAPE>::numVertices>
Edge<SHAPE, V>::next_ = ShapeProp<SHAPE>::next;

//! Local index of previous vertex
template <corlib::shape SHAPE, typename V>
const std::array<int,ShapeDim<SHAPE>::numVertices>
Edge<SHAPE, V>::previous_ = ShapeProp<SHAPE>::previous;

//! Give opposite edge
template <corlib::shape SHAPE, typename V>
inline Edge<SHAPE, V>*&
Edge<SHAPE, V>::opposite(const Vertex_* v)
{    
    FTL_VERIFY(v != NULL);
    
    unsigned indexOpposite = numNeighbors;
    for (unsigned i = 0; i < numVertices; i++){
        if (v == vertices_[i]){
            indexOpposite = i;
            break;
        }
    }
    FTL_VERIFY(indexOpposite!= numNeighbors);
    return neighbors_[indexOpposite];
}

