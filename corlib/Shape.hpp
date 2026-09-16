// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//
//                 Computational Structural Mechanics Lab
//                         University of Cambridge
//
// This software is copyrighted and all rights are retained by the University of
// Cambridge. All use, disclosure, and/or reproduction of any part not expressly
// authorised by the University of Cambridge is prohibited. (C) 2010
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

//! @file Shape.hpp
//! @todo Do mesh size computations for ALL element types (paper of Hughes?)

#ifndef corlib_shape_h
#define corlib_shape_h
//------------------------------------------------------------------------------

#include <iostream>
#include <string>
#include <array>

#include <boost/algorithm/string.hpp>

#include <Eigen/Core>

#include <corlib/linalg.hpp>
#include <corlib/verify.hpp>
#include <corlib/fuzzyEqual.hpp>
#include <corlib/eigenX.hpp>

//------------------------------------------------------------------------------
namespace corlib{

    //--------------------------------------------------------------------------
    //! \brief Enumerator for description of the geometry of an element
    enum shape {
        UNDEFINED     = -1,  //!< undefined shape (a default)
        POINT,               //!< a point
        LINE,                //!< a line element
        TRIANGLE,            //!< a triangle element
        QUADRILATERAL,       //!< a quadrilateral element
        TETRAHEDRON,         //!< a tetrehedral element
        HEXAHEDRON,          //!< a hexahedral element
        PIXEL,               //!< a tensor-product 2D cell
        VOXEL                //!< a tensor-product 3D cell
    };

    //--------------------------------------------------------------------------
    /** \struct corlib::ShapeTraits 
     *  \brief Group of shape-specfic traits
     *  \details These traits allow access the natural dimension of a geometric
     *   shape, its number of corner points (vertices), its number of faces,
     *   the centroid in reference coordinates and the size of the corresponding
     *   reference element.
     */
    template<shape SHAPE> struct ShapeTraits;
    
    // -----------------------------------------
    //! \cond SKIPDOX
    template<> struct ShapeTraits<POINT>
    {
        static const unsigned dim          = 0;
        static const unsigned numVertices  = 1;
    };
    
    //--------------------------------------------------------------------------
    template<> struct ShapeTraits<LINE>
    {
        static const unsigned dim          = 1;
        static const unsigned numVertices  = 2;
        static const unsigned numFaces     = 2;
        static const shape    faceShape    = POINT;
        static eigenX::VectorSd<dim> centroid() {
            return 1./2. * eigenX::VectorSd<dim>::Ones();
        }
        static double refSize() { return 1.; }
        static bool isInside( const eigenX::VectorSd<dim> xi ) {
            return ( fuzzyGreaterEqual( xi(0), 0.0 ) and
                        fuzzyLessEqual( xi(0), 1.0 ) ); }
    };

    // -----------------------------------------
    template<> struct ShapeTraits<TRIANGLE>
    {
        static const unsigned dim          = 2;
        static const unsigned numVertices  = 3;
        static const unsigned numFaces     = 3;
        static const shape    faceShape    = LINE;
        static eigenX::VectorSd<dim> centroid() {
            return 1./3. * eigenX::VectorSd<dim>::Ones();
        } 
        static double refSize() { return 0.5; }
        static bool isInside( const eigenX::VectorSd<dim> xi ) {
            return ( fuzzyGreaterEqual( xi(0), 0.0       ) and
                        fuzzyLessEqual( xi(0), 1.0       ) and 
                     fuzzyGreaterEqual( xi(1), 0.0       ) and
                        fuzzyLessEqual( xi(1), 1.0-xi(0) ) );
        }
    };

    // -----------------------------------------
    template<> struct ShapeTraits<QUADRILATERAL>
    {
        static const unsigned dim          = 2;
        static const unsigned numVertices  = 4;
        static const unsigned numFaces     = 4;
        static const shape    faceShape    = LINE;
        static eigenX::VectorSd<dim> centroid() {
            return 1./2. * eigenX::VectorSd<dim>::Ones();
        }
        static double refSize() { return 1.; }
        static bool isInside( const eigenX::VectorSd<dim> xi ) {
            return ( fuzzyGreaterEqual( xi(0), 0.0 ) and
                        fuzzyLessEqual( xi(0), 1.0 ) and 
                     fuzzyGreaterEqual( xi(1), 0.0 ) and
                        fuzzyLessEqual( xi(1), 1.0 ) );
        }
    };

    // -----------------------------------------
    template<> struct ShapeTraits<TETRAHEDRON>
    {
        static const unsigned dim          = 3;
        static const unsigned numVertices  = 4;
        static const unsigned numFaces     = 4;
        static const shape    faceShape    = TRIANGLE;
        static eigenX::VectorSd<dim> centroid() {
            return 1./4. * eigenX::VectorSd<dim>::Ones();
        } 
        static double refSize() { return 1./6.; }
        static bool isInside( const eigenX::VectorSd<dim> xi ) {
            return ( fuzzyGreaterEqual( xi(0), 0.0             ) and
                        fuzzyLessEqual( xi(0), 1.0             ) and 
                     fuzzyGreaterEqual( xi(1), 0.0             ) and
                        fuzzyLessEqual( xi(1), 1.0-xi(0)       ) and 
                     fuzzyGreaterEqual( xi(2), 0.0             ) and 
                        fuzzyLessEqual( xi(2), 1.0-xi(0)-xi(1) ) );
        }
    };

    // -----------------------------------------
    template<> struct ShapeTraits<HEXAHEDRON>
    {
        static const unsigned dim          = 3;
        static const unsigned numVertices  = 8;
        static const unsigned numFaces     = 6;
        static const shape    faceShape    = QUADRILATERAL;
        static eigenX::VectorSd<dim> centroid() {
            return 1./2. * eigenX::VectorSd<dim>::Ones();
        }
        static double refSize() { return 1.; }
        static bool isInside( const eigenX::VectorSd<dim> xi ) {
            return ( fuzzyGreaterEqual( xi(0), 0.0 ) and
                        fuzzyLessEqual( xi(0), 1.0 ) and 
                     fuzzyGreaterEqual( xi(1), 0.0 ) and
                        fuzzyLessEqual( xi(1), 1.0 ) and 
                     fuzzyGreaterEqual( xi(2), 0.0 ) and 
                        fuzzyLessEqual( xi(2), 1.0 ) );
        }
    };

    // -----------------------------------------
    template<> struct ShapeTraits<PIXEL>
    {
        static const unsigned dim          = 2;
        static const unsigned numVertices  = 4;
        static const unsigned numFaces     = 4;
        static const shape    faceShape    = LINE;
        static eigenX::VectorSd<dim> centroid() {
            return 1./2. * eigenX::VectorSd<dim>::Ones();
        }
        static double refSize() { return 1.; }
    };

    // -----------------------------------------
    template<> struct ShapeTraits<VOXEL>
    {
        static const unsigned dim          = 3;
        static const unsigned numVertices  = 8;
        static const unsigned numFaces     = 6;
        static const shape    faceShape    = PIXEL;
        static eigenX::VectorSd<dim> centroid() {
            return 1./2. * eigenX::VectorSd<dim>::Ones();
        }
        static double refSize() { return 1.; }
    };
    //! \endcond

    //------------------------------------------------------------------------------
    //! Assign the shape of a Simplex to a space dimension
    template<unsigned DIM> struct Simplex;
    //!\cond SKIPDOX
    template<> struct Simplex<0> { static const shape theShape = POINT;       };
    template<> struct Simplex<1> { static const shape theShape = LINE;        };
    template<> struct Simplex<2> { static const shape theShape = TRIANGLE;    };
    template<> struct Simplex<3> { static const shape theShape = TETRAHEDRON; };
    //! \endcond

    //------------------------------------------------------------------------------
    //! Assign the shape of a tensor-product element
    template<unsigned DIM> struct HyperCube;
    //! \cond SKIPDOX
    template<> struct HyperCube<1> { static const shape theShape = LINE; };
    template<> struct HyperCube<2> { static const shape theShape = PIXEL; };
    template<> struct HyperCube<3> { static const shape theShape = VOXEL; };
    //! \endcond

    //------------------------------------------------------------------------------
    //! compute a critical mesh parameter given the vertices and element shape
    template<shape SHAPE>
    double elementSize( const std::array<Eigen::VectorXd,
                                         ShapeTraits<SHAPE>::numVertices> & X );

    /** Compute \f$ h = | x_1 - x_0 | \f$ of a line element
     *  \param[in] X  Matrix where each column is a vertex     */
    template<> 
    double elementSize<LINE>( const std::array< Eigen::VectorXd, 2> & X )
    {
        return ( X[1] - X[0] ).norm( );
    }

    /** Compute the radius of the circumcircle of a triangle by using 
     *  Heron's formula:
     *  \f[
     *      r = \frac{(a+b+c)(a+b-c)(a+c-b)(b+c-a)}{4}
     *  \f]
     *  where a,b, and c are the side lenghts of the triangle.
     *  \param[in] X  Matrix where each column is a vertex of the triangle   */
    template<> 
    double elementSize<TRIANGLE>( const std::array< Eigen::VectorXd, 3> & X )
    {
        // lenghts of the sides
        const double a( ( X[1] - X[0] ).norm() );
        const double b( ( X[2] - X[1] ).norm() );
        const double c( ( X[0] - X[2] ).norm() );
        // Heron's formula
        const double area = std::sqrt( (a+b+c)*(a+b-c)*(a+c-b)*(b+c-a) ) / 4.;
        // return radius of circumcircle
        return ( a * b * c ) / 4. / area;
    }

    /** Radius of the circumsphere of a tetrahedron by means of the formula:
     *  \f[
     *   r = \frac{| <a,a>(b \times c) + <b,b>(c \times a) + <c,c> (a \times b) |}
     *       { ( 12  V ) }
     *  \f]
     *
     *   with a,b,c the vectors spanning the tetrahedron (starting from one vertex)
     *   and \f$ V = | <a, (b \times c) > | / 6 \f$, the volume of the element   
     *  \param[in] X  Matrix where each column is a vertex     */
    template<> 
    double elementSize<TETRAHEDRON>( const std::array< Eigen::VectorXd, 4> & X )
    {
        // vectors which span the tetrahedron and their lenghts squared
    	Eigen::VectorXd a( X[1] - X[0] ); const double a2 = a.squaredNorm();
    	Eigen::VectorXd b( X[2] - X[0] ); const double b2 = b.squaredNorm();
    	Eigen::VectorXd c( X[3] - X[0] ); const double c2 = c.squaredNorm();
        // 12 times the volume of the tetrahedron
        const double V12 = 2. * a.dot( corlib::cross_prod( b, c ) );

        // return radius of sircumcircle;
        const eigenX::VectorSd<3> t = ( a2 * corlib::cross_prod( b, c ) +
                                        b2 * corlib::cross_prod( c, a ) +
                                        c2 * corlib::cross_prod( a, b ) );
        return t.norm( ) / V12;
    }

    template<> 
    double elementSize<QUADRILATERAL>( const std::array< Eigen::VectorXd, 4> & X )
    {
        //diagonals
        const double el1 = ( X[2] - X[0] ).norm( );
        const double el2 = ( X[3] - X[1] ).norm( );
        if( el1 > el2 ) return el1;
        return el2;
    }

    /****************
     * Still missing: the hexehedral mesh size parameters 
     * (circumthings??)
     */

    //------------------------------------------------------------------------------
    //! Convert (input) shape string to enum
    enum shape convertShapeStringToEnum( const std::string & ss )
    {
        enum shape se = UNDEFINED;
        const std::string ssL = boost::to_lower_copy( ss );
    
        if      ( ssL == "undefined"     ) se = UNDEFINED;
        else if ( ssL == "point"         ) se = POINT;
        else if ( ssL == "line"          ) se = LINE;
        else if ( ssL == "triangle"      ) se = TRIANGLE;
        else if ( ssL == "quadrilateral" ) se = QUADRILATERAL;
        else if ( ssL == "tetrahedron"   ) se = TETRAHEDRON;
        else if ( ssL == "hexahedron"    ) se = HEXAHEDRON;
        else if ( ssL == "pixel"         ) se = PIXEL;
        else if ( ssL == "voxel"         ) se = VOXEL;
        else {
            std::cerr << "Found unrecognized shape name: " << ss << std::endl;
            exit(0);
        }
        return se;
    }

    //! Convert shape enum to string
    std::string convertShapeEnumToString( const enum shape se )
    {
        std::string ss = "undefined";
        
        if      ( se == UNDEFINED     ) ss = "undefined";
        else if ( se == POINT         ) ss = "point";
        else if ( se == LINE          ) ss = "line";
        else if ( se == TRIANGLE      ) ss = "triangle";
        else if ( se == QUADRILATERAL ) ss = "quadrilateral";
        else if ( se == TETRAHEDRON   ) ss = "tetrahedron";
        else if ( se == HEXAHEDRON    ) ss = "hexahedron";
        else if ( se == PIXEL         ) ss = "pixel";
        else if ( se == VOXEL         ) ss = "voxel";
        else FTL_VERIFY_DESCRIPTIVE( false, "Cannot handle shape enum" );

        return ss;
    }

 
    //! Give the dimension of shape. This function is for cases when static 
    //! compilation is not sufficient (i.e. dynamic input).
    unsigned dimensionOfShape( const enum shape& shape)
    {
        // initialise output
        unsigned dim = 0;

        switch ( shape ) {
        case LINE:
            dim = ShapeTraits< LINE >::dim;
            break;
        case TRIANGLE :
            dim = ShapeTraits< TRIANGLE >::dim;
            break;
        case QUADRILATERAL :
            dim = ShapeTraits< QUADRILATERAL >::dim;
            break;
        case TETRAHEDRON : 
            dim = ShapeTraits< TETRAHEDRON >::dim;
            break;
        case HEXAHEDRON :
            dim = ShapeTraits< HEXAHEDRON >::dim;
            break;
        case PIXEL :
            dim = ShapeTraits< PIXEL >::dim;
            break;
        case VOXEL :
            dim = ShapeTraits< VOXEL >::dim;
            break;
        default:
            FTL_VERIFY_DESCRIPTIVE( false, "shape type not known." );
        }

        return dim;
    }
    
} 

#endif
