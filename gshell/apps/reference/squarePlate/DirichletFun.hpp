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

#ifndef gshell_apps_dirichletfun_h 
#define gshell_apps_dirichletfun_h

#include <Eigen/Core>

#include <corlib/eigenX.hpp>
#include <corlib/Constraints.hpp>

//==============================================================================
// declarations
namespace app{

    namespace eigenX = corlib::eigenX;
    
    corlib::NodalConstraint dirichletSimply(
        const eigenX::VectorSd<3> X
        );

    corlib::NodalConstraint enforceZeroShears(
        const eigenX::VectorSd<3> X
        );

}

//==============================================================================
// definitions

//------------------------------------------------------------------------------
corlib::NodalConstraint app::dirichletSimply(
    const eigenX::VectorSd<3> X
    )
{
    corlib::NodalConstraint constraint;

    // simply supported
    if ( corlib::fuzzyEqual( X[ 0 ],  0., 1.e-5 ) or
         corlib::fuzzyEqual( X[ 1 ],  0., 1.e-5 ) or
         corlib::fuzzyEqual( X[ 0 ], 10., 1.e-5 ) or
         corlib::fuzzyEqual( X[ 1 ], 10., 1.e-5 ) ) {
        constraint.setComponent( 0, 0. );
        constraint.setComponent( 1, 0. );
        constraint.setComponent( 2, 0. );
        //constraint.setComponent( 3, 0. );
        //constraint.setComponent( 4, 0. );
    } 

    return constraint;
}

//------------------------------------------------------------------------------
corlib::NodalConstraint app::enforceZeroShears(
    const eigenX::VectorSd<3> X
    )
{
    corlib::NodalConstraint constraint;

    constraint.setComponent( 3, 0. );
    constraint.setComponent( 4, 0. );
    
    return constraint;
}

#endif


