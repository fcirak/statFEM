// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//
//                 Computational Structural Mechanics Lab
//                         University of Cambridge
//
// This software is copyrighted and all rights are retained by the University of
// Cambridge. All use, disclosure, and/or reproduction of any part not expressly
// authorised by the University of Cambridge is prohibited. (C) 2010
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

//! @file fuzzyEqual.hpp

#ifndef corlib_fuzzyEqual_h
#define corlib_fuzzyEqual_h
//------------------------------------------------------------------------------
#include <cmath>

namespace corlib{

    //! Comparison of two doubles \f$ a \approx b \f$
    //! up to a tolerance
    //!
    //! \param[in]  a      Number a
    //! \param[in]  b      Number b
    //! \param[in]  tol    Tolerance 
    //! \return            Result
    bool fuzzyEqual( const double & a, const double & b, 
                     const double & tol = 1.e-8 )
    {
        return ( std::fabs( a - b ) < tol );
    }

    //! Comparison of two doubles \f$ a \stackrel{<}{\approx} b \f$
    //! up to a tolerance
    //!
    //! \param[in]  a      Number a
    //! \param[in]  b      Number b
    //! \param[in]  tol    Tolerance 
    //! \return            Result
    bool fuzzyLessEqual( const double & a, const double & b, 
                         const double & tol = 1.e-8 )
    {
        return ( ( std::fabs( a - b ) < tol ) or ( a < b ) );
    }

    //! Comparison of two doubles \f$ a \stackrel{>}{\approx} b \f$
    //! up to a tolerance
    //!
    //! \param[in]  a      Number a
    //! \param[in]  b      Number b
    //! \param[in]  tol    Tolerance 
    //! \return            Result
    bool fuzzyGreaterEqual( const double & a, const double & b, 
                            const double & tol = 1.e-8 )
    {
        return ( ( std::fabs( a - b ) < tol ) or ( a > b ) );
    }

}

#endif
