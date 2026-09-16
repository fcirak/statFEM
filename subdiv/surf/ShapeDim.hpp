//
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//
//                 Computational Structural Mechanics Lab
//                         University of Cambridge
//
// This software is copyrighted and all rights are retained by the CSMLab. All
// use, disclosure, and/or reproduction of any part not expressly authorized by
// F Cirak is prohibited. (C) 2010.
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//

#ifndef subdiv_surf_shapedim_h
#define subdiv_surf_shapedim_h

#include <iostream>
#include <string>
#include <array>

#include <corlib/Shape.hpp>

namespace subdiv{
    namespace surf{

        //======================================================================
        //! gives access to the inherent properties of each shape
        template< corlib::shape SHAPE >
        struct ShapeProp
        {
            static const std::array< int, corlib::ShapeTraits< SHAPE >::numVertices > next;
            static const std::array< int, corlib::ShapeTraits< SHAPE >::numVertices > previous;
        };

        template<>
        const std::array< int, 3 >
        ShapeProp< corlib::TRIANGLE >::next = { {1, 2, 0} };

        template<> 
        const std::array< int, 3 > 
        ShapeProp< corlib::TRIANGLE >::previous = { {2, 0, 1} };

        template<> 
        const std::array< int, 4 > 
        ShapeProp< corlib::QUADRILATERAL >::next = { {1, 2, 3, 0} };

        template<> 
        const std::array< int, 4 > 
        ShapeProp< corlib::QUADRILATERAL >::previous = { {3, 0, 1, 2} };

    } // end namespace surf
} // end namespace subdiv

#endif
