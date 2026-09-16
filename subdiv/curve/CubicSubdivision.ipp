#include <cassert>

#include <corlib/verify.hpp>
#include <subdiv/curve/MeshTags.hpp>
#include <subdiv/curve/collect.hpp>

using namespace subdiv::curve;

//! Local index of next vertex
template <typename ET>
const std::array<int, ET::numVertices>
Subdivision<CUBIC, ET>::next_ = ShapeProp<corlib::LINE>::next;

//! Local index of previous vertex
template <typename ET>
const std::array<int, ET::numVertices>
Subdivision<CUBIC, ET>::previous_ = ShapeProp<corlib::LINE>::previous;

template<typename ET>
void Subdivision<CUBIC, ET>::refineEdge(ETree* et, int level, VecP& posNew)
{
    Vertex* v0 = et->vertex(0);
    Vertex* v1 = et->vertex(1);
    VecP pos0 = v0->point(level);
    VecP pos1 = v1->point(level);

    posNew += pos0;
    posNew += pos1;
    posNew /= 2;    

    return;
}

template<typename ET>
void Subdivision<CUBIC, ET>::refineVertex(ETree* et, int ivtx, int level, VecP& posNew)
{
    Vertex* v=et->vertex(ivtx);
    std::vector<Vertex*> vtxs;
    
    if (v->tag(level)== subdiv::curve::NOSUB_VERTEX)
        posNew =  et->vertex(ivtx)->point(level);
    else {
        subdiv::curve::CollectOneRing(et, ivtx, vtxs);
        assert(vtxs.size()==2);
        
        posNew += vtxs[0]->point(level);
        posNew += vtxs[1]->point(level);	    
        posNew *= 0.125;
        
        VecP p0 = v->point(level);
        p0 *= 0.75;	     
        posNew += p0;
    }	   
    return;
}

/// Return tangent vectors to limit surce
template<typename ET>
void Subdivision<CUBIC, ET>::computeTangents( ETree * et, const int ivtx, 
                                             const int level, MapVPtrD & coeff)
{


}

