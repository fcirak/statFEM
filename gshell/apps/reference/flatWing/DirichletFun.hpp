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

#ifndef gshell_app_dirichletfun_h 
#define gshell_app_dirichletfun_h

#include <Eigen/Core>
#include <corlib/eigenX.hpp>
#include <corlib/Constraints.hpp>
#include "RotationFlapping.hpp"

//==============================================================================
// declarations
namespace app{

    namespace eigenX = corlib::eigenX;

    corlib::NodalConstraint supportNothing(
        const eigenX::VectorSd<3> x
        );

    corlib::NodalConstraint supportFlapping(
        const eigenX::VectorSd<3> x, 
        const double period,
        const double time,
        const double timeStepSize
        );

}

//==============================================================================
// definitions

//------------------------------------------------------------------------------
corlib::NodalConstraint app::supportNothing(
    const eigenX::VectorSd<3> x
    )
{
    corlib::NodalConstraint constraint;
    return constraint;
}

//------------------------------------------------------------------------------
corlib::NodalConstraint app::supportFlapping(
    const eigenX::VectorSd<3> x, 
    const double period,
    const double time,
    const double timeStepSize
    )
{
    namespace ublas = boost::numeric::ublas;

    // initialise return object
    corlib::NodalConstraint constraint;

    //               ^ y
    //    ___________|
    //   /           |
    //  <          --+------> x
    //   `  wing   / |
    // tip `______/   root
    
    const bool onBoundary = ( x[ 0 ] > -1.0 );

    if ( onBoundary ) {

        typedef eigenX::MatrixSd<3,3> Mat3x3;
        const Mat3x3 rot0 = flappingRotMatrix( period, time, NULL );  // t_n
        const Mat3x3 rot1 = flappingRotMatrix( period, time+timeStepSize, NULL );  // t_{n+1}

        typedef eigenX::VectorSd<3> Vec3;
        const Vec3 increment = ( rot1 * x ) - ( rot0 * x );

        for ( unsigned d = 0; d < 3; ++d )
            constraint.setComponent( d, increment[ d ] );
    }

    return constraint;
}




#endif


