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

#ifndef subdiv_surf_buildchildren_h
#define subdiv_surf_buildchildren_h

#include<Eigen/Core>

namespace subdiv {
    namespace surf {

        //======================================================================
        /// Make children of facet
        ///
        /// \tparam   SHAPE   Simplex shape, i.e. triangle or quadrilateral
        /// \tparam   FT      Facet (tree) type
        template< corlib::shape SHAPE, typename FT >
        class BuildChildren;

        //----------------------------------------------------------------------
        template< typename FT >
        class BuildChildren< corlib::TRIANGLE, FT >
        {
        public:
            static const corlib::shape myShape = corlib::TRIANGLE;

            typedef FT                                           FTree;
            typedef typename FTree::Vertex                       Vertex;
            typedef typename FTree::VecVPtrNV                    VecVPtrNV;
            typedef typename FTree::VecVPtr                      VecVPtr;

        public:
            /// Make children
            ///
            /// \param[in,out]  ft   Facet which gets children
            ///                      (not needed only for interface compatibility)
            /// \param[in]      fv   Refined vertex at facet centre
            /// \param[in]      ev   Refined vertices on edges
            /// \param[in,out]  nv   Newly allocated vertices
            void operator()( FTree * ft, Vertex * fv, const VecVPtrNV & ev, VecVPtr & nv );

        private:
            /// Permutation: Gives local index of next node depending on current node index
            static const std::array< int, corlib::ShapeTraits< myShape >::numVertices >       next_;
            /// Permutation: Gives local index of previous node depending on current node index
            static const std::array< int, corlib::ShapeTraits< myShape >::numVertices >       previous_;
        };

        //----------------------------------------------------------------------
        template< typename FT >
        class BuildChildren< corlib::QUADRILATERAL, FT >
        {
        public:
            static const corlib::shape myShape = corlib::QUADRILATERAL;

            typedef FT                                           FTree;
            typedef typename FTree::Vertex                       Vertex;
            typedef typename FTree::VecVPtrNV                    VecVPtrNV;
            typedef typename FTree::VecVPtr                      VecVPtr;

        public:
            /// Make children
            ///
            /// \param[in,out]  ft   Facet which gets children
            /// \param[in]      fv   Refined vertex at facet centre
            /// \param[in]      ev   Refined vertices on edges
            /// \param[in,out]  nv   Newly allocated vertices
            void operator()( FTree * ft, Vertex * fv, const VecVPtrNV & ev, VecVPtr & nv );

        private:
             /// Permutation: Gives local index of next node depending on current node index
             static const std::array< int, corlib::ShapeTraits< myShape >::numVertices >       next_;
             /// Permutation: Gives local index of previous node depending on current node index
             static const std::array< int, corlib::ShapeTraits< myShape >::numVertices >       previous_;
        };

    }
}

//------------------------------------------------------------------------------
#include "BuildChildren.ipp"
#endif
