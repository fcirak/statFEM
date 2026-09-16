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

//! @author Burkhard Bornemann, Kosala Bandara
//! @date   2010

#ifndef gshell_fem_misc_h
#define gshell_fem_misc_h

#include <functional>

//==============================================================================
// declarations
namespace gshell {
    namespace fem {

        //----------------------------------------------------------------------
        /// Helper function to allocate object (with e.g. std::generate)
        template< typename T >
        T * newPointer() { return new T; }

        //----------------------------------------------------------------------
        /// Helper function to get inner product of two vectors
        template<typename VEC>
        double norm2Squared( const VEC & v )
        {
            static_assert(( std::is_same<typename VEC::value_type,
                                                 double>::value ));
            return std::inner_product( v.begin(), v.end(), v.begin(), 0. );
        }

        //----------------------------------------------------------------------
        /// Helper function to get Euclidian norm of vector and is bindable
        template<typename VEC>
        double norm2( const VEC & v )
        {
            static_assert(( std::is_same<typename VEC::value_type,
                                                 double>::value ));
            return std::sqrt( gshell::fem::norm2Squared<VEC>( v ) );
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
