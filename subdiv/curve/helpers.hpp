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

#ifndef subdiv_curve_helpers_h
#define subdiv_curve_helpers_h

#include <cmath> 

#include <corlib/linalg.hpp>
#include <corlib/Shape.hpp>

#include <subdiv/curve/ShapeDim.hpp>

#include <Eigen/Core>
#include <corlib/eigenX.hpp>

namespace subdiv{
    namespace curve{
    namespace eigenX = corlib::eigenX;
    
    template<unsigned P>
	class Factorial; 	    

	template <unsigned N, unsigned R>
	class NcR;

    template<corlib::shape SHAPE>
    double criticalDimension(const Eigen::VectorXd& xi);
    //double criticalDimension(const eigenX::VectorSd<ShapeDim<SHAPE>::dim>& xi);
    
    }
}


// Factorial
namespace subdiv{
    namespace curve{
        template <unsigned P>
        class Factorial {
        public:
            static const unsigned result=P*Factorial<P-1>::result;
        };
    }
}

namespace subdiv{
    namespace curve{
	template <>
	class Factorial<0> {
	public:
	    static const unsigned result = 1;
	};
    }
}


// helper templates for computing nCr
namespace subdiv{
    namespace curve{
        template <unsigned N, unsigned R>
        class NcR {
        public:
            static const unsigned n1 = Factorial<N>::result;
            static const unsigned n2 = Factorial<N-R>::result;
            static const unsigned n3 = Factorial<R>::result;
            
            static const unsigned result= n1/(n2*n3);            
        };
    }
}


//!find critical dimension in line element
namespace subdiv{
    namespace curve{
	template <>
        double criticalDimension<corlib::LINE>(const Eigen::VectorXd & xi)
        {    
            double xiC;
            //FTL_VERIFY( xi.rows( ) == ShapeDim<SHAPE>::dim );

            double allDims[] = {xi(0), 1 - xi(0)};
            xiC = *std::min_element(allDims, allDims + 2);

            return xiC;
            
        }
    }
}
#endif
