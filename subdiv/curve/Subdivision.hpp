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

#ifndef subdiv_curve_subdivision_h
#define subdiv_curve_subdivision_h

#include <iostream>
#include <string>
#include <utility>
#include <vector>
#include <map>
#include <array>

#include <corlib/Shape.hpp>

#include <subdiv/curve/ShapeDim.hpp>
#include <subdiv/curve/MeshTags.hpp>
#include <subdiv/curve/BuildChildren.hpp>


namespace subdiv{
    namespace curve{	        

        //! \brief Enumerator for description of subdivision method
        enum method {
            CUBIC
        };
        
        //! Subdivision rules 
        template<method METHOD, typename ET> 
        class Subdivision;

        //! Specialisation for cubic
        template<typename ET> 
        class Subdivision<CUBIC, ET>;
    }
}
 
       
//----------------------------------------------------------------------
template<typename ET>
class subdiv::curve::Subdivision<subdiv::curve::CUBIC, ET>
{
public:
    static const subdiv::curve::method myMethod     = subdiv::curve::CUBIC;
    static const corlib::shape myShape           = corlib::LINE;          
    static const unsigned     degree             = 3;
    static const unsigned     ctrlPtsRing        = 1;
    static const unsigned     supportRing        = 2;
    static const unsigned     numVertices        = 2; 
    static const unsigned     numVerticesSN      = 4; 
    
    typedef ET                                      ETree;
    typedef typename ET::Vertex                     Vertex; 
    
    typedef typename Vertex::VecP                   VecP;
    
    typedef std::map< Vertex *, double >            MapVPtrD;
    typedef typename MapVPtrD::iterator             MapVPtrDIter;
    
    typedef subdiv::curve::vertexTag                VertexTag;
    
public:            
    //empty constructor
    Subdivision(){}
    
    //refine edge in given mesh
    void refineEdge(ETree* et, int level, VecP& posNew);
    
    //refine vertex in given mesh
    void refineVertex(ETree* et, int ivtx, int level, VecP& posNew);
    
    /// Return tangent vectors to limit surface
    void computeTangents( ETree * ft, const int ivtx, 
                          const int level, MapVPtrD & coeff);
    
private:
    static const std::array<int, ET::numVertices>       next_;
    static const std::array<int, ET::numVertices>       previous_;
    
};




//------------------------------------------------------------------------------

#include <subdiv/curve/CubicSubdivision.ipp>

#endif
