// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//
//                 Computational Structural Mechanics Lab
//                         University of Cambridge
//
// This software is copyrighted and all rights are retained by the University of
// Cambridge. All use, disclosure, and/or reproduction of any part not expressly
// authorised by the University of Cambridge is prohibited. (C) 2010
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

//! @file helpers.hpp

#ifndef solid_apps_newmark_helpers_h
#define soild_apps_newmark_helpers_h

//------------------------------------------------------------------------------
// headers
#include <corlib/verify.hpp>
#include <corlib/fuzzyEqual.hpp>

//==============================================================================
// declarations
namespace app{

    namespace ublas = boost::numeric::ublas;

    ublas::bounded_vector< double, 2 > 
    gravity( const ublas::bounded_vector< double, 2 > coorRef );

    double scalePresribedDisplacement( const double time, const double dt );

}

//==============================================================================
// implementations

//------------------------------------------------------------------------------
// Body force function (has to be given the coordinate as argument)
boost::numeric::ublas::bounded_vector< double, 2 > 
app::gravity( const ublas::bounded_vector< double, 2 > coorRef )
{
    ublas::bounded_vector< double, 2 > weigDens;
    weigDens( 0 ) =  0.1;
    weigDens( 1 ) =  0.0;
    return weigDens;
}

//------------------------------------------------------------------------------
// Function to scale prescribed displacements depending on time
double app::scalePresribedDisplacement( const double time, const double dt )
{
    const double factor = corlib::fuzzyEqual( time, dt ) ? 1. : 0.;
    return factor;
}

#endif
