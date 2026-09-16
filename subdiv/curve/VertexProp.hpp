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

#ifndef subdiv_curve_vertex_prop_h
#define subdiv_curve_vertex_prop_h

#include <Eigen/Core>
#include <corlib/eigenX.hpp>

namespace subdiv{
    namespace curve{
        template <typename V>
	    class VertexProp;	    

        namespace eigenX = corlib::eigenX;
    }
}

template <typename V>
class subdiv::curve::VertexProp  {
public:
    typedef V                            Vertex;

    typedef typename Vertex::VecP        VecP;
    typedef subdiv::curve::vertexTag     VertexTag;

public:
    VertexProp(int lv=0)
        :level(lv), tag(subdiv::curve::NOTAG_VERTEX),
         point(eigenX::VectorSd< Vertex::numVariables >::Zero()){ }

public:
    int                level;
    VertexTag          tag;
    VecP               point;
};

#endif
