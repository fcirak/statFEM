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

#ifndef solid_application_newmark_loads_h
#define solid_application_newmark_loads_h

//------------------------------------------------------------------------------
// headers
#include <corlib/verify.hpp>

//------------------------------------------------------------------------------
// declarations
boost::numeric::ublas::bounded_vector< double, 2 > 
gravity( const boost::numeric::ublas::bounded_vector< double, 2 > coorRef );

//------------------------------------------------------------------------------
// implementations

//------------------------------------------------------------------------------
// Body force function (has to be given the coordinate as argument)
boost::numeric::ublas::bounded_vector< double, 2 > 
gravity( const boost::numeric::ublas::bounded_vector< double, 2 > coorRef )
{
    boost::numeric::ublas::bounded_vector< double, 2 > weigDens;
    weigDens( 0 ) =  0.1;
    weigDens( 1 ) =  0.0;
    return weigDens;
}

#endif
