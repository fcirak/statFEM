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

#ifndef subdiv_curve_collect_h
#define subdiv_curve_collect_h

#include <vector>
#include <set>

#include <subdiv/curve/BuildChildren.hpp>
#include <subdiv/curve/ShapeDim.hpp>

namespace subdiv{
    namespace curve{
        /// Collect one ring of a vertex (ordered)
        template< typename ETREE >
        void CollectOneRing(ETREE *et, int ivtx, 
                            std::vector<typename ETREE::Vertex*>& oneRingV);

        /// Collect one ring of an edge (unordered)
        template< typename ETREE >
        void CollectOneRingEdge(ETREE *et, std::vector<ETREE*>& oneRingE,
                                std::set<typename ETREE::Vertex*>& oneRingV);

    }
}

// collect
template< typename ETREE >
void subdiv::curve::CollectOneRing(ETREE *et, int ivtx, 
                                std::vector<typename ETREE::Vertex*>& oneRingV)
{
    oneRingV.push_back(et->vertex(ShapeProp<corlib::LINE>::next[ivtx]));

    ETREE* en = et->neighbor(ivtx);
    assert(en!=NULL);
    int j = localVtxNumber(en, et->vertex(ivtx));
    oneRingV.push_back(en->vertex(ShapeProp<corlib::LINE>::next[j]));

    return;
}


template< typename ETREE >
void subdiv::curve::CollectOneRingEdge(ETREE *et, std::vector<ETREE*>& oneRingE, 
                                std::set<typename ETREE::Vertex*>& oneRingV)
{
    
    typedef ETREE                           ETree;
    typedef typename ETREE::VecVPtrNVItr    VecVPtrNVItr;
    typedef std::vector<ETREE*>             VecETPtr;
    typedef typename VecETPtr::iterator     VecETPtrItr;

    const unsigned numNeighbors = ETree::numNeighbors;

    for ( unsigned n = 0; n < numNeighbors; n++ ) {
        ETree * en = et->neighbor( n );
        if ( en )
            oneRingE.push_back( en );
    }

    // collect vertices
    for ( VecETPtrItr e = oneRingE.begin(); e != oneRingE.end(); ++e ) {
        oneRingV.insert( (*e)->vBegin( ), (*e)->vEnd( ) );
    }
    for ( VecVPtrNVItr v = et->vBegin( ); v != et->vEnd( ); ++v ) {
        oneRingV.erase( *v );
    }

    return;
}


#endif
