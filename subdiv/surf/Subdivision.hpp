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

#ifndef subdiv_surf_subdivision_h
#define subdiv_surf_subdivision_h

#include <iostream>
#include <string>
#include <map>
#include <array>

#include <corlib/verify.hpp>
#include <corlib/Shape.hpp>

#include <corlib/TensorSpline.hpp>
#include <subdiv/surf/BoxSpline.hpp>
#include <subdiv/surf/ShapeDim.hpp>
#include <subdiv/surf/MeshTags.hpp>

//------------------------------------------------------------------------------
namespace subdiv{
    namespace surf{

        //----------------------------------------------------------------------
        /// Subdivision method enumerator
        enum method {
            CATMULL_CLARK,   ///< Catmull-Clark's scheme for quadrilateral simplexes
            LOOP             ///< Loop's scheme scheme for triangular simplexes
        };

        //======================================================================
        /// give subdivision shape
        ///
        /// \tparam  METHOD     Subdivision method, e.g. Loop's scheme
        template< method METHOD > 
        struct SubdivisionShape;

        template<>
        struct SubdivisionShape< CATMULL_CLARK >
        {
            static const corlib::shape myShape = corlib::QUADRILATERAL;
        };

        template<>
        struct SubdivisionShape< LOOP >
        {
            static const corlib::shape myShape = corlib::TRIANGLE;
        };

        //======================================================================
        /// Subdivision rules 
        ///
        /// \tparam  METHOD     Subdivision method, e.g. Loop's scheme
        /// \tparam  FT         Facet tree type (simplex type sampling surface)
        /// \tparam  DEGREE     Spline degree (for NURBS)
        template< method METHOD, typename FT, unsigned DEGREE = 0 > 
        class Subdivision;

        //----------------------------------------------------------------------
        /// Subdivision rule specialisation for Loop's scheme
        ///
        /// Reference:
        /// - C Loop, Smooth subdivision surfaces based on triangles, PhD thesis, 1987.
        ///
        /// \tparam  FT          Facet tree type (triangles sampling surface)
        template< typename FT >
        class Subdivision< subdiv::surf::LOOP, FT >
        {
        public:
            static const method        myMethod          = subdiv::surf::LOOP;
            static const corlib::shape myShape           = corlib::TRIANGLE;
            static const unsigned      degree            = 4;
            static const unsigned      numRings          = 1;
            static const unsigned      numExtraRings     = 0;
            static const unsigned      numVertices       = 3;
            static const unsigned      numFunctions      = 12;
            
            typedef subdiv::surf::BoxSpline                 Spline;
            typedef typename FT::Vertex                     Vertex; 
            typedef typename Vertex::Point                  Point;

            typedef std::map< Vertex *, double >            MapVPtrD;
            typedef typename MapVPtrD::iterator             MapVPtrDIter;

        private:
            typedef std::vector< Vertex * >                 VecVPtr_;

        public:            
            /// empty constructor
            Subdivision( ) { }
            
            /// Refine facet in given mesh: Face point
            ///
            /// Actually, does nothing, because the facets are not refined
            /// in Loop's subdivision scheme.
            ///
            /// \param[in]     ft      Facet tree entry which is refined
            /// \param[in]     level   Base refinement level from which is refined
            /// \param[out]    posNew  New co-ordinate for refined vertex at centre
            ///                        of facet.  Must be initialised to zero.
            void refineFacet( FT * ft, const int level, Point & posNew );

            /// Refine edge in given mesh: Edge point
            ///
            /// Computes co-ordinate for a new vertex in the middle of the edge
            ///
            /// \param[in]     ft      Facet tree entry which is refined
            /// \param[in]     iedge   Local edge index
            /// \param[in]     level   Base refinement level from which is refined
            /// \param[out]    posNew  New co-ordinate for refined vertex at middle
            ///                        of edge.  Must be initialised to zero.
            void refineEdge( FT * ft, const int iedge, const int level, Point & posNew );
            
            /// Refine vertex in given mesh
            ///
            /// Computes co-ordinate of a vertex which already exists (update)
            ///
            /// \param[in]     ft      Facet tree entry which is refined
            /// \param[in]     ivtx    Local index of vertex
            /// \param[in]     level   Base refinement level from which is refined
            /// \param[out]    posNew  New co-ordinate for refined vertex at middle
            ///                        of edge.  Must be initialised to zero.
            void refineVertex( FT * ft, const int ivtx, const int level, Point & posNew );

            /// Return tangent vectors to limit surface
            ///
            /// Computes the tangents to the limit surface at vertex of facet.
            /// This takes the formulas contained in the appendix of
            /// - H Biermann, A Levin, D Zorin, Piecewise smooth subdivision surfaces
            ///       with normal control, SIGGRAPH 2000 conference proceedings, p 113-120, 2000.
            ///
            /// \param[in]     ft      Facet tree entry which is refined
            /// \param[in]     ivtx    Local index of vertex
            /// \param[in]     level   Base refinement level from which is refined
            /// \param[out]    coeff0  Coefficients of first tangent
            /// \param[out]    coeff1  Coefficients of second tangent
            void computeTangents( FT * ft, const int ivtx, const int level,
                                  MapVPtrD & coeff0, MapVPtrD & coeff1 );

            /// Return position on limit surface
            ///
            /// Computes the position on the limit surface at vertex of facet.
            /// This takes the formulas contained in the appendix of
            /// - H Biermann, A Levin, D Zorin, Piecewise smooth subdivision surfaces
            ///       with normal control, SIGGRAPH 2000 conference proceedings, p 113-120, 2000.
            /// - C Loop, Smooth subdivision surfaces based on triangles, PhD thesis, 1987.
            /// - F Cheng, F Fan, S Lai, C Huang, J Wang, J Yong, Loop Subdivision Surface
            ///       based Progressive Interpolation, J of Comp Sc and Tech 24(1), p 39–46, 2009.
            ///
            /// \param[in]     ft      Facet tree entry which is refined
            /// \param[in]     ivtx    Local index of vertex
            /// \param[in]     level   Base refinement level from which is refined
            /// \param[out]    coeff   Coefficients of limit point
            void computeLimit( FT * ft, const int ivtx, const int level,
                               MapVPtrD & coeff );

        private:
            static const std::array< int, numVertices >       next_;
            static const std::array< int, numVertices >       previous_;

        };

        //----------------------------------------------------------------------
        /// Catmull-Clark subdivision scheme
        ///
        /// \tparam  FT          Facet tree type (quadrilaterals sampling surface)
        template< typename FT >
        class Subdivision< subdiv::surf::CATMULL_CLARK, FT >
        {
        public:
            static const method        myMethod          = subdiv::surf::CATMULL_CLARK;
            static const corlib::shape myShape           = corlib::QUADRILATERAL;
            /// Spline degree: Here cubic
            static const unsigned      degree            = 3;
            static const unsigned      numRings          = 1;
            static const unsigned      numExtraRings     = 0;
            /// Number of corner vertices of simplex
            static const unsigned      numVertices       = 4;
            /// Number of shape functions whose supports intersects with facet domain
            static const unsigned      numFunctions      = 16; 
            
            typedef corlib::TensorSpline< degree, 2 >       Spline;
            typedef typename FT::Vertex                     Vertex; 
            typedef typename Vertex::Point                  Point;

            typedef std::map< Vertex *, double >            MapVPtrD;
            typedef typename MapVPtrD::iterator             MapVPtrDIter;

        private:
            typedef std::vector< Vertex * >                 VecVPtr_;

        public:            
            /// Empty constructor
            Subdivision( ) { }
            
            /// Refine facet in given mesh: Face point
            ///
            /// Computes co-ordinate for a new vertex in the centre of the facet
            ///
            /// \param[in]     ft      Facet tree entry which is refined
            /// \param[in]     level   Base refinement level from which is refined
            /// \param[out]    posNew  New co-ordinate for refined vertex at centre
            ///                        of facet.  Must be initialised to zero.
            void refineFacet( FT * ft, const int level, Point & posNew );

            /// Refine edge in given mesh: Edge point
            ///
            /// Computes co-ordinate for a new vertex in the middle of the edge
            ///
            /// \param[in]     ft      Facet tree entry which is refined
            /// \param[in]     iedge   Local edge index
            /// \param[in]     level   Base refinement level from which is refined
            /// \param[out]    posNew  New co-ordinate for refined vertex at middle
            ///                        of edge.  Must be initialised to zero.
            void refineEdge( FT * ft, const int iedge, const int level, Point & posNew );
            
            /// Refine vertex in given mesh
            ///
            /// Computes co-ordinate of a vertex which already exists (update)
            ///
            /// \param[in]     ft      Facet tree entry which is refined
            /// \param[in]     ivtx    Local index of vertex
            /// \param[in]     level   Base refinement level from which is refined
            /// \param[out]    posNew  New co-ordinate for refined vertex at middle
            ///                        of edge.  Must be initialised to zero.
            void refineVertex( FT * ft, const int ivtx, const int level, Point & posNew );

            /// Return tangent vectors to limit surface
            ///
            /// Computes the tangents to the limit surface at vertex of facet.
            /// This takes the formulas contained in the appendix of
            /// H Biermann, A Levin, D Zorin, Piecewise smooth subdivision surfaces
            /// with normal control, SIGGRAPH 2000 conference proceedings, p 113-120,
            /// July 2000.
            ///
            /// \param[in]     ft      Facet tree entry which is refined
            /// \param[in]     ivtx    Local index of vertex
            /// \param[in]     level   Base refinement level from which is refined
            /// \param[out]    coeff0  Coefficients of first tangent
            /// \param[out]    coeff1  Coefficients of second tangent
            void computeTangents( FT * ft, const int ivtx, const int level,
                                  MapVPtrD & coeff0, MapVPtrD & coeff1 );

            /// Return position on limit surface
            ///
            /// Computes the position on the limit surface at vertex of facet.
            /// This takes the formulas contained in the appendix of
            /// - H Biermann, A Levin, D Zorin, Piecewise smooth subdivision surfaces
            ///       with normal control, SIGGRAPH 2000 conference proceedings, p 113-120,
            ///       July 2000.
            ///
            /// \param[in]     ft      Facet tree entry which is refined
            /// \param[in]     ivtx    Local index of vertex
            /// \param[in]     level   Base refinement level from which is refined
            /// \param[out]    coeff   Coefficients of limit point
            void computeLimit( FT * ft, const int ivtx, const int level,
                               MapVPtrD & coeff );
            
        private:
            /// Permutation: Local index of next vertex
            static const std::array< int, numVertices >       next_;
            /// Permutation: Local index of previous vertex
            static const std::array< int, numVertices >       previous_;

        };

    }
}


//==============================================================================
#include "SubdivisionLoop.ipp"
#include "SubdivisionCatmullClark.ipp"

#endif
