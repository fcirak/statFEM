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

#ifndef subdiv_curve_compare_functor_h
#define subdiv_curve_compare_functor_h

#include <utility>

namespace subdiv{
    namespace curve{
	template<typename T>
	class CompareFunctor;	
        
        // comparison functor for pairs
        template <typename N>
        class CompareFunctor<std::pair<int,N> > {
        public:
            bool operator ()(const std::pair<int,N>& lhs,const std::pair<int,N>& rhs)const 
                {
                    int lhsI = lhs.first;
                    int rhsI = rhs.first;

                    N lhsN = lhs.second;
                    N rhsN = rhs.second;

                    return ( lhsI < rhsI or
                             (lhsI == rhsI and
                              lhsN < rhsN ) );
                 
                }
         
        };

	// template<typename T>
	// class CompareFunctorId;	

        // // id comparison functor for pointers
        // template<typename T>
        // class CompareFunctorId 
        // {
        // public:
        //     bool operator ()(const T& lhs,const T& rhs)const 
        //         {
        //             std::cout<<lhs->index()<<" "<<rhs->index()<<std::endl;
        //             return lhs->index() < rhs->index();
        //         }
        // };

    }
}        
#endif

