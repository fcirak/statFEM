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

#ifndef subdiv_curve_shapedim_h
#define subdiv_curve_shapedim_h

#include <iostream>
#include <string>
#include <array>

#include <corlib/Shape.hpp>

namespace subdiv{
    namespace curve{
	//! gives access to the inherent dimension of each shape
	template<corlib::shape SHAPE> 
        struct ShapeDim;

        //! gives access to the inherent properties of each shape
	template<corlib::shape SHAPE> 
        struct ShapeProp
        {
            static const std::array<int, 
                ShapeDim<SHAPE>::numVertices>      next;
            static const std::array<int, 
                ShapeDim<SHAPE>::numVertices>      previous;
        };

	template<> 
        struct ShapeDim<corlib::LINE> 
	{
	    static const unsigned dim           = 1; 
	    static const unsigned numVertices   = 2;
	    static const unsigned numChildren   = 2;
	    static const unsigned numNeighbors  = 2;
	};

        template<> 
        const std::array<int,2> 
        subdiv::curve::ShapeProp<corlib::LINE>::next = { {1, 0} };

        template<> 
        const std::array<int,2> 
        subdiv::curve::ShapeProp<corlib::LINE>::previous = { {1, 0} };

    } // end namespace curve
} // end namespace subdiv

#endif
