// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//
//                 Computational Structural Mechanics Lab
//                         University of Cambridge
//
// This software is copyrighted and all rights are retained by the University of
// Cambridge. All use, disclosure, and/or reproduction of any part not expressly
// authorised by the University of Cambridge is prohibited. (C) 2010
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

//! @file Loads.hpp

#ifndef beam_application_ring_loads_h
#define beam_application_ring_loads_h

//------------------------------------------------------------------------------
// headers
#include <Eigen/Core>
#include <corlib/eigenX.hpp>
#include <corlib/verify.hpp>

namespace eigenX = corlib::eigenX;

//------------------------------------------------------------------------------
// declarations
eigenX::VectorSd<2>
gravity( const eigenX::VectorSd<2> coorRef );

eigenX::VectorSd<2>
centralForce( const eigenX::VectorSd<2> coorRef );

eigenX::VectorSd<2>
windMill( const eigenX::VectorSd<2> coorRef );

//------------------------------------------------------------------------------
// implementations

//------------------------------------------------------------------------------
// Body force function (has to be given the coordinate as argument)
eigenX::VectorSd<2>
gravity( const eigenX::VectorSd<2> coorRef )
{
    eigenX::VectorSd<2> weigDens;
    weigDens( 0 ) = 0.0;  // N/m^3
    weigDens( 1 ) = - 7.85e3 * 9.81;  // Kg/m^3 * N/Kg
    return weigDens;
}

//------------------------------------------------------------------------------
// central force
eigenX::VectorSd<2>
centralForce( const eigenX::VectorSd<2> coorRef )
{
    FTL_VERIFY_DESCRIPTIVE( coorRef.norm( )>1.0e-6,
			    "Hit singular point" );
    const double scale = 0.1/0.02;
    eigenX::VectorSd<2> forceDens;
    forceDens = scale*coorRef/( coorRef.norm( ) );
    return forceDens;
}

//------------------------------------------------------------------------------
// normal traction follows windmill pattern,
//       
//        |Y
//        |
//   *         *bisectrix
//    *    p  *
//     *     *
//      *   *
//       *|*
// p =>   +-- <= p --------> X
//       * *
//      *   *
//     *  ^  *
//    *   p   *bisectrix
// 
eigenX::VectorSd<2>
windMill( const eigenX::VectorSd<2> coorRef )
{
    FTL_VERIFY_DESCRIPTIVE( coorRef.norm( )>1.0e-6,
                            "Hit singular point" );
    // reference co-ordinates (convenience)
    const double X = coorRef( 0 );
    const double Y = coorRef( 1 );

    // body load
    const double scale = 0.0000001/0.02;

    // body load
    eigenX::VectorSd<2> forceDens =
        Eigen::VectorXd::Zero( 2 );
    if ( Y <= X ) {
        if ( Y >= -X ) {
            forceDens( 0 ) = -scale;
        }
        else {
            forceDens( 1 ) = scale;
        }
    }
    else {
        if ( Y >= -X ) {
            forceDens( 1 ) = -scale;
        }
        else {
            forceDens( 0 ) = scale;
        }
    }
    return forceDens;
}

#endif
