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

#ifndef gshell_app_load_h 
#define gshell_app_load_h

#include <Eigen/Core>
#include <corlib/eigenX.hpp>

//==============================================================================
// declarations
namespace app{
 
    namespace eigenX = corlib::eigenX;
   
    eigenX::VectorSd<3> pressure(
            const eigenX::VectorSd<3> x,
            const double thickness
    );

}

//==============================================================================
// definitions

//------------------------------------------------------------------------------
/// Constant body force in global z-direction
corlib::eigenX::VectorSd<3> app::pressure(
    const eigenX::VectorSd<3> x,
    const double thickness
    )
{
    const double t = thickness;   // m
    const double pRes = 1.e-3;  // N/m^2

    eigenX::VectorSd<3> p;  // N/m^3
    p.setZero();
//    p[ 2 ] = - pRes / t;
    return p;
}

#endif
