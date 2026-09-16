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

#ifndef beam_application_cantilever_loads_h
#define beam_application_cantilever_loads_h

//------------------------------------------------------------------------------
// headers
#include <Eigen/Core>
#include <corlib/eigenX.hpp>

namespace eigenX = corlib::eigenX;

//------------------------------------------------------------------------------
// Body force function (has to be given the coordinate as argument)
eigenX::VectorSd< 2 >
gravity( const eigenX::VectorSd< 2 > coorRef )
{
    eigenX::VectorSd< 2 > weigDens;
    weigDens( 0 ) = 0.0;
    weigDens( 1 ) = -0.01;
    return weigDens;
}


#endif
