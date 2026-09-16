// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//
//                 Computational Structural Mechanics Lab
//                         University of Cambridge
//
// This software is copyrighted and all rights are retained by the University of
// Cambridge. All use, disclosure, and/or reproduction of any part not expressly
// authorised by the University of Cambridge is prohibited. (C) 2010
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

//! @file misc.hpp

#ifndef beam_fem_misc_h
#define beam_fem_misc_h

//! System includes
#include <functional>
#include <numeric>
#include <cmath>
#include <type_traits>
//! Corlib includes
#include <corlib/verify.hpp>

//==============================================================================
// declarations
namespace beam {
    namespace fem {

        //----------------------------------------------------------------------
        /// Helper function to get inner product of a vector with itself
        template<typename VEC>
        double norm2Squared( const VEC & v )
        {
            static_assert( (std::is_same<typename VEC::value_type, double>::value),
                                   "Value type has to be double");
            return std::inner_product( v.begin(), v.end(), v.begin(), 0. );
        }

        //----------------------------------------------------------------------
        /// Helper function to get Euclidian norm of vector and is bindable
        template<typename VEC>
        double norm2( const VEC & v )
        {
            static_assert( (std::is_same<typename VEC::value_type, double>::value),
                                   "Value type has to be double" );
            return std::sqrt( beam::fem::norm2Squared<VEC>( v ) );
        }

        //----------------------------------------------------------------------
        /// Helper function to get max norm of vector and is bindable
        template<typename VEC>
        double normMax( const VEC & v )
        {
            static_assert( (std::is_same<typename VEC::value_type, double>::value),
                                   "Value type has to be double" );
            double(*dabs)(double) = static_cast<double(*)(double)>( &std::abs );
            return 
                dabs( *std::max_element( v.begin(), v.end(),
                                         std::bind( std::less<double>(),
                                                      std::bind( dabs, std::placeholders::_1 ),
                                                      std::bind( dabs, std::placeholders::_2 ) ) ) );
        }

        //----------------------------------------------------------------------
        /// Multiplication (eg scalar multiplication) of an item of type S
        /// (eg scalar) with an instance of type V (eg vector) resulting in an
        /// object of type W (eg vector again)
        ///
        /// \tparam W  Type of result
        /// \tparam S  Type of first term
        /// \tparam V  Type of second term
        template<typename W, typename S, typename V>
        struct Multiply : std::function<W(const S&, const V&)>
        {
            /// Return \f$w = s * v\f$
            W operator()( const S & s, const V & v ) { return (s * v); }
        };

    }
}

#endif
