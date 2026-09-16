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

#ifndef subdiv_surf_shape_h
#define subdiv_surf_shape_h

#include <corlib/Shape.hpp>

#include <Eigen/Core>
#include <corlib/eigenX.hpp>

#include <subdiv/surf/ShapeDim.hpp>
#include <subdiv/surf/MeshTags.hpp>

namespace subdiv{
    namespace surf{

    	namespace eigenX = corlib::eigenX;

        //======================================================================
        /// Compute minimal (characteristic) distance of parametric
        /// co-ordinate to parametric domain boundary
        ///
        /// \tparam      SHAPE  Simplex shape type 
        ///
        /// \param[in]   xi     Parametric co-ordinate
        /// \return             Minimal distance
       	template< corlib::shape SHAPE >
        double shapeMinDistanceToBoundary(
            const eigenX::VectorSd< corlib::ShapeTraits< SHAPE >::dim > & xi
            );

        template< >
        double shapeMinDistanceToBoundary< corlib::TRIANGLE >(
            const eigenX::VectorSd< corlib::ShapeTraits< corlib::TRIANGLE >::dim > & xi
            );

        template< >
        double shapeMinDistanceToBoundary< corlib::QUADRILATERAL >(
            const eigenX::VectorSd< corlib::ShapeTraits< corlib::QUADRILATERAL >::dim > & xi
            );


        //======================================================================
        ///
        /// \tparam   SHAPE   Simplex shape, i.e. triangle or quadrilateral
        /// \tparam   FACET      Facet (tree) type
        template< corlib::shape SHAPE, typename FACET >
        class ShapeChildCoordinate;

        //----------------------------------------------------------------------
        template< typename FACET >
        class ShapeChildCoordinate< corlib::TRIANGLE, FACET >
        {
        public:
            static const corlib::shape myShape = corlib::TRIANGLE;
            static const unsigned lDim         = corlib::ShapeTraits< myShape >::dim;

            typedef FACET                           Facet;
            typedef eigenX::VectorSd<lDim>          VecLDim;
            typedef eigenX::MatrixSd<lDim,lDim>     MatLDimLDim;
            
        public:
            /// Determine which sub-facet contains evaluation point xi
            /// and transform xi to co-ordinate system for sub-facet
            ///
            /// \param[in]       facet   Parent facet which gets children
            /// \param[in,out]   xi      Parametric evaluation point
            /// \param[in,out]   jac     Jacobian at parametric evaluation point
            /// \param[out]      child   Index of sub-facet
            void operator()( Facet * facet,
                             VecLDim & xi,
                             MatLDimLDim & jac,
                             unsigned & child );
        };

        //----------------------------------------------------------------------
        template< typename FACET >
        class ShapeChildCoordinate< corlib::QUADRILATERAL, FACET >
        {
        public:
            static const corlib::shape myShape = corlib::QUADRILATERAL;
            static const unsigned numVertices  = corlib::ShapeTraits< myShape >::numVertices;
            static const unsigned lDim         = corlib::ShapeTraits< myShape >::dim;

            typedef FACET                           Facet;
            typedef eigenX::VectorSd<lDim>          VecLDim;
            typedef eigenX::MatrixSd<lDim,lDim>     MatLDimLDim;
            
        public:
            /// Determine which sub-facet contains evaluation point xi
            /// and transform xi to co-ordinate system for sub-facet
            ///
            /// \param[in]       facet   Parent facet which gets children
            /// \param[in,out]   xi      Parametric evaluation point
            /// \param[in,out]   jac     Jacobian at parametric evaluation point
            /// \param[out]      child   Index of sub-facet
            void operator()( Facet * facet,
                             VecLDim & xi,
                             MatLDimLDim & jac,
                             unsigned & child );
        };

    }
}

//------------------------------------------------------------------------------
#include <subdiv/surf/ShapeFunHelpers.ipp>

#endif
