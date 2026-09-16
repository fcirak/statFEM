// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//
//                 Computational Structural Mechanics Lab
//                         University of Cambridge
//
// This software is copyrighted and all rights are retained by the University of
// Cambridge. All use, disclosure, and/or reproduction of any part not expressly
// authorised by the University of Cambridge is prohibited. (C) 2010
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

//! @file PotentialQuadrature.hpp

#ifndef del2_fem_potentialquadrature_h
#define del2_fem_potentialquadrature_h

//------------------------------------------------------------------------------
// headers
#include <corlib/Shape.hpp>

//------------------------------------------------------------------------------
// forward declaration
namespace del2 {
    namespace fem {

        //======================================================================
        //! Provide number of default/optimal Gauss points
        //! required for integration
        //!
        //! These traits allow to select the number of quadrature (or Gauss)
        //! points based on the element shape and the polynomial order
        //! of the shape functions.
        //!
        //! \tparam SHAPE    basic shape of element
        //! \tparam NUMFUNC  number of shape functions
        template< corlib::shape SHAPE, unsigned NUMFUNC >
        struct PotentialQuadratureTraits;

        //----------------------------------------------------------------------
        // specialisations
        template< >
        struct PotentialQuadratureTraits< corlib::TRIANGLE, 3 >
        {
            static const unsigned numPoints = 1;
        };

        template< >
        struct PotentialQuadratureTraits< corlib::TRIANGLE, 6 >
        {
            static const unsigned numPoints = 3;
        };

        template< >
        struct PotentialQuadratureTraits< corlib::QUADRILATERAL, 4 >
        {
            static const unsigned numPoints = 4;
        };

        template< >
        struct PotentialQuadratureTraits< corlib::QUADRILATERAL, 8 >
        {
            static const unsigned numPoints = 9;
        };

        template< >
        struct PotentialQuadratureTraits< corlib::QUADRILATERAL, 9 >
        {
            static const unsigned numPoints = 9;
        };

    }
}

#endif
