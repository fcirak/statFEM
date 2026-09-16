// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//
//                 Computational Structural Mechanics Lab
//                         University of Cambridge
//
// This software is copyrighted and all rights are retained by the University of
// Cambridge. All use, disclosure, and/or reproduction of any part not expressly
// authorised by the University of Cambridge is prohibited. (C) 2010
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

//! @file QuadratureTraitsStatic.hpp

#ifndef solid_fem_quadraturetraitsstatic_h
#define solid_fem_quadraturetraitsstatic_h

namespace solid {
    namespace fem {

        //======================================================================
        //! This traits allows to deduce the (total) number of Gauss
        //! points depending on simplex shape (triangle, quadrilateral, hexahedron, ...)
        //! and the number of used Lagrangian shape functions.
        //!
        //! \tparam  SHAPE    Simplex geometry
        //! \tparam  NFUN     Number of shape functions
        template< enum corlib::shape SHAPE, unsigned NFUN >
        struct QuadratureTraitsStatic;

        template< >
        struct QuadratureTraitsStatic< corlib::TRIANGLE, 3 >
        {
            static const unsigned numPoints = 1;
        };
        
        template< >
        struct QuadratureTraitsStatic< corlib::TRIANGLE, 6 >
        {
            static const unsigned numPoints = 6;
        };
        
        template< >
        struct QuadratureTraitsStatic< corlib::QUADRILATERAL, 4 >
        {
            static const unsigned numPoints = 4;
        };
        
        template< >
        struct QuadratureTraitsStatic< corlib::QUADRILATERAL, 9 >
        {
            static const unsigned numPoints = 9;
        };

        template< >
        struct QuadratureTraitsStatic< corlib::HEXAHEDRON, 8 >
        {
            static const unsigned numPoints = 8;
        };

        template< >
        struct QuadratureTraitsStatic< corlib::HEXAHEDRON, 27 >
        {
            static const unsigned numPoints = 27;
        };

    }
}

#endif
