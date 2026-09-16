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

#ifndef subdiv_curve_build_children_h
#define subdiv_curve_build_children_h

#include <corlib/Shape.hpp>
#include <corlib/verify.hpp>

#include <subdiv/curve/ShapeDim.hpp>
#include <subdiv/curve/MeshTags.hpp>

// helper function
namespace subdiv {
    namespace curve {

        //======================================================================
        /// Make children of edge
        ///
        /// \tparam   ET      Edge (tree) type
        template<typename ET>
        class BuildChildren;

        //! Get local index of a vertex within the edge tree
        template<typename V, typename ET>
        int localVtxNumber(ET* et, V* vtx);

	//--------------------------------------------------------------------------
        //! Get local index of a vertex within edge
        template<typename V, typename ET>
        int localVtxNumber(ET *et, V* vtx){
            int j = (et->vertex(0) == vtx) ? 0 : ((et->vertex(1) == vtx) ? 1: 2);    
            FTL_VERIFY(j!=2);
            return j;
            
        }
    }
}

//--------------------------------------------------------------------------
/// Make children of edge
template<typename ET>
class subdiv::curve::BuildChildren
{
public:
    typedef ET                                           ETree;
    typedef typename ETree::Vertex                       Vertex;
    
public:
    /// Make children
    void operator()( ETree * et, Vertex * ev);
    
private:
            /// Permutation: Gives local index of next node depending on current node index
    static const std::array< int, ET::numVertices >       next_;
    /// Permutation: Gives local index of previous node depending on current node index
    static const std::array< int, ET::numVertices >       previous_;
};
#include <subdiv/curve/BuildChildren.ipp>

#endif
