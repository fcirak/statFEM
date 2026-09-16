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

#ifndef gshell_app_rotationflapping_h 
#define gshell_app_rotationflapping_h

#include <corlib/linalg.hpp>

//==============================================================================
// declarations
namespace app{

    namespace eigenX = corlib::eigenX;

    eigenX::MatrixSd< 3, 3 > flappingRotMatrix(
        const double strokePeriod,
        const double time,
        eigenX::MatrixSd< 3, 3 > * rotMatDeriv
        );

    eigenX::VectorSd< 3 > flappingDisplacements(
        const eigenX::VectorSd< 3 > & x,
        const double & period,
        const double & time
        );

    eigenX::VectorSd< 3 > flappingVelocities(
        const eigenX::VectorSd< 3 > & x,
        const double & period,
        const double & time
        );

};

//==============================================================================
// definitions

//------------------------------------------------------------------------------
// Rotation matrix due to prescribed flapping motion
corlib::eigenX::MatrixSd< 3, 3 > app::flappingRotMatrix(
    const double strokePeriod,
    const double time,
    eigenX::MatrixSd< 3, 3 > * rotMatDeriv
    )
{
    namespace eigenX = corlib::eigenX;
    
    typedef eigenX::VectorSd< 3 >       Vec3;
    typedef eigenX::MatrixSd< 3, 3 >    Mat3x3;

    // constants
    const double deviationAngleAmplitude = 30.*M_PI/180.;  // rad

    // matrices are prepared for root-to-tip axis in _negative_ x-direction
    //
    //               ^ y
    //    ___________|
    //   /           |
    //  <          --+------> x
    //   `  wing   / |
    // tip `______/   root

    // deviation rotation
    const double deviationAngle = deviationAngleAmplitude *
        .5 * ( 1. - std::cos( 2.*M_PI*time/(.5*strokePeriod) ) );
    Vec3 aY; aY.setZero( ); aY( 1 ) = 1.;
    const Mat3x3 rotY = corlib::rotationMatrix( aY, -deviationAngle );
    const Mat3x3 rot = rotY;
        
    // time-derivative of matrix
    if ( rotMatDeriv ) {
        const double deviationAngleD = deviationAngleAmplitude *
            .5 * std::sin( 2.*M_PI*time/(.5*strokePeriod) ) * 2.*M_PI/(.5*strokePeriod);
        Vec3 bY; bY.setZero( );
        const Mat3x3 derY = corlib::rotationMatrixDeriv( aY, -deviationAngle,
                                                         bY, -deviationAngleD );
        *rotMatDeriv = derY;
    }        

    // return matrix
    return rot;
}

//------------------------------------------------------------------------------
// prescribed displacements due to flapping motion
corlib::eigenX::VectorSd< 3 > app::flappingDisplacements(
    const eigenX::VectorSd< 3 > & x,
    const double & period,
    const double & time
    )
{
    namespace eigenX = corlib::eigenX;

    // rotation matrix
    const eigenX::MatrixSd< 3, 3 > rot =
        flappingRotMatrix( period, time, NULL );

    // perform rotation
    const eigenX::VectorSd< 3 > dis = ( rot * x ) - x;
    return dis;
}

//------------------------------------------------------------------------------
// prescribed velocities due to flapping motion
corlib::eigenX::VectorSd< 3 > app::flappingVelocities(
    const eigenX::VectorSd< 3 > & x,
    const double & period,
    const double & time
    )
{
    namespace eigenX = corlib::eigenX;

    // time-derivative of rotation matrix
    eigenX::MatrixSd< 3, 3 > rotDeriv;
    flappingRotMatrix( period, time, &rotDeriv );

    // perform rotation
    const eigenX::VectorSd< 3 > vel = ( rotDeriv * x );
    return vel;
}



#endif
