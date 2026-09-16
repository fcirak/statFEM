// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//
//                 Computational Structural Mechanics Lab
//                         University of Cambridge
//
// This software is copyrighted and all rights are retained by the University of
// Cambridge. All use, disclosure, and/or reproduction of any part not expressly
// authorised by the University of Cambridge is prohibited. (C) 2010
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

//! @file ShapefunTraits.hpp

#ifndef corlib_shapefuntraits_h
#define corlib_shapefuntraits_h

#include <corlib/Shape.hpp>
#include <corlib/verify.hpp>
#include <corlib/misc.hpp>

// forward declarations
namespace corlib{

    template<corlib::shape SHAPE, unsigned NNODS>
    class Shapefun;

    template<unsigned DEGREE, unsigned DIM>
    class TensorSpline;
}

//------------------------------------------------------------------------------
namespace corlib{

    //! Descriptor for the type of shape functions
    enum SfunClass {
        NOSHAPEFUN,      //!< Undefined
        SIMPLEX,         //!< d-Simplex shape function
        TENSORSPLINE,    //!< Tensor-product of B-Splines
        TENSORLAGRANGE,  //!< Tensor-product of Lagrangians
        SERENDIPITY      //!< Serendipity shape functions
    };

    template<corlib::shape SHAPE, unsigned numNodes> struct ShapeFunClassifier;

    //--------------------------------------------------------------------------
    //! \cond SKIPDOX
    namespace detail_{

        //======================================================================
        // Number of nodes on a simplex element for given degree
        // The number of nodes of a p-degree d-simplex are
        // N = (1/d!)  prod_{i=1}^d (p+i)
        // The value of N is here implemented recursively
        template<unsigned DEGREE, unsigned DIM>
        struct NumSimplexNodes
        {
            static const unsigned value = 
                (DEGREE + DIM ) * NumSimplexNodes<DEGREE,DIM-1>::value / DIM;
        };
    
        template<unsigned DEGREE>
        struct NumSimplexNodes<DEGREE,0>
        {
            static const unsigned value = 1;
        };

        //======================================================================
        // Degree of a simplex element, given its number of nodes
        template<unsigned NNODS, unsigned DIM>
        class SimplexDegree
        {
            // observation shows that for D > 1 the number of nodes is
            // N = ( (p-1)^D + l.o.t + D!) / D!
            // ignoring the lower-order terms, this yields
            // p = floor[ (N * D! - D!)^{1/D} ] - 1
        private:
            static const unsigned facD = corlib::factorial<DIM>::value;
            static const unsigned rad  = NNODS * facD - facD;
            static const unsigned root = corlib::NthRoot<rad,DIM>::value;
        public:
            static const unsigned value = root - 1;
        };

        // obvious result for 1D: p = N - 1
        template<unsigned NNODS>
        class SimplexDegree<NNODS,1>
        {
        public:
            static const unsigned value = NNODS - 1;
        };

        //======================================================================
        // Tensor-product Lagrange element shapes
        template<unsigned DIM> struct HyperLagrange;

        template<>
        struct HyperLagrange<1> 
        { 
            static const corlib::shape value = corlib::LINE; 
        };

        template<>
        struct HyperLagrange<2> 
        { 
            static const corlib::shape value = corlib::QUADRILATERAL; 
        };

        template<>
        struct HyperLagrange<3> 
        { 
            static const corlib::shape value = corlib::HEXAHEDRON; 
        };

        //======================================================================
        //  Lagrangian shape functions on a simplex 
        //  SfunClass = SIMPLEX
        //  Example DIM=2, DEGREE=2, NNODS=6:
        //                                                   
        //       2                     
        //       |`                     
        //       |  `                  
        //       5    `4              
        //       |      `             
        //       |        `           
        //       0-----3----1 
        //    

        //----------------------------------------------------------------------
        // Type of simplex shape function for given degree and dimension
        template<unsigned DEGREE, unsigned DIM>
        struct SimplexShapeFunction
        {
            static const corlib::shape theShape = corlib::Simplex<DIM>::theShape;
            static const unsigned      numNodes = NumSimplexNodes<DEGREE,DIM>::value;
            
            typedef typename corlib::Shapefun<theShape,numNodes> type;
        };

        //----------------------------------------------------------------------
        // Degree of simplex shape function for given number of nodes and dim
        template<unsigned NNODS, unsigned DIM>
        struct SimplexShapeFunctionDegree
        {
            static const unsigned value = SimplexDegree<NNODS,DIM>::value;
        };


        //======================================================================
        // Tensor spline shape function on a hyper-cube
        // SfunClass = TENSORSPLINE
        // Example  DIM=2, DEGREE=2, NNODS=9:
        //
        //    6        7        8
        //
        //        +---------+
        //        |         |
        //    3   |    4    |   5
        //        |         |
        //        +---------+
        //
        //    0        1        2        
        //

        //----------------------------------------------------------------------
        // Type of Tensor-product spline shape function for given degree and dim
        template<unsigned DEGREE, unsigned DIM>
        struct TensorSplineShapeFunction
        {
            typedef typename corlib::TensorSpline<DEGREE,DIM> type;
        };

        //----------------------------------------------------------------------
        // Degree of Tensor-product spline shape funs given the number of nodes
        template<unsigned NNODS, unsigned DIM>
        struct TensorSplineShapeFunctionDegree
        {
            static const unsigned value = corlib::NthRoot<NNODS,DIM>::value-1;
        };

        //======================================================================
        // Serendipity shape function on  quad/hexes
        // SfunClass = SERENDIPITY
        // Example DIM=2, DEGREE=2, NNODS=8:
        //
        //        3----6----2
        //        |         |
        //        7         5    
        //        |         |
        //        0----4----1

        //----------------------------------------------------------------------
        // Type of Serendipity shape functions given the degree and dimension
        template<unsigned DEGREE, unsigned DIM>
        struct SerendipityShapeFunction
        {
            // number of nodes is numVertices + (p-1)* numEdges =
            //                     2^d + (p-1) (2^d * d /2) = 2^{d-1}*(2+d*(p-1))
            static const unsigned      numNodes = 
                corlib::twoToTheN<DIM-1>::value * (2 + (DEGREE-1) * DIM);
            static const corlib::shape theShape = HyperLagrange<DIM>::value;

            typedef typename corlib::Shapefun<theShape,numNodes> type;
        };

        //----------------------------------------------------------------------
        // Degree of Serendipity shape functions given the number of nodes
        template<unsigned NNODS, unsigned DIM>
        struct SerendipityShapeFunctionDegree
        {
            // given above formula for number of nodes, it can be reversed to
            // p = (N/2^{d-1} - 2)/d + 1
            static const unsigned aux = 
                (NNODS/corlib::twoToTheN<DIM-1>::value) - 2;
            static const unsigned value = aux / DIM + 1;
        };

        //======================================================================
        // Tensor-produc Lagrangian on quads/hexes
        // SfunClass = TENSORLAGRANGE
        // Example DIM=2, DEGREE=2, NNODS=9
        //
        //        3----6----2
        //        |         |
        //        7    8    5    
        //        |         |
        //        0----4----1

        //----------------------------------------------------------------------
        // Type of Tensor-product Lagrange shape functions given the degree
        template<unsigned DEGREE, unsigned DIM>
        struct TensorLagrangeShapeFunction
        {
            static const corlib::shape theShape = HyperLagrange<DIM>::value;
            static const unsigned numNodes = corlib::mToTheN<DEGREE+1,DIM>::value;
            typedef typename corlib::Shapefun<theShape,numNodes> type;
        };

        //----------------------------------------------------------------------
        // Degree of tensor-product Lagrangian shape functions
        template<unsigned NNODS, unsigned DIM>
        struct TensorLagrangeShapeFunctionDegree
        {
            static const unsigned value = corlib::NthRoot<NNODS,DIM>::value - 1;
        };

    }//! \endcond SKIPDOX
}


//------------------------------------------------------------------------------
/** \brief Deduce shape function type from shape and number of nodes
 *  \details By checking the values of shape and number of nodes, 
 *  this object classifies the shape function type as one  of the 
 *  corlib::SfunClass values. 
 *  Note that this classification fails for 1D, since in this case all shapes
 *  are a line, the degree is the number of nodes minus one. It is thus 
 *  impossible to distinguish between the Lagrangian and Tensor-Spline types. 
 *  \tparam SHAPE
 *  \tparam NNODS
 */
template<corlib::shape SHAPE, unsigned NNODS>
struct corlib::ShapeFunClassifier
{
    // deduce dimension from the shape
    static const unsigned dim      = corlib::ShapeTraits<SHAPE>::dim;
    static const unsigned numNodes = NNODS;

    // check for tensor-product type: dim-th power of dim-th integer root 
    static const unsigned rootNumNodes =  corlib::NthRoot<numNodes,dim>::value;
    static const bool isTensor         = (corlib::mToTheN<rootNumNodes,dim>::value == numNodes);

    // check if shape is a hyperquad (QUADRILATERAL|HEXAHEDRON)
    static const bool isHyperLagrange  = (SHAPE == detail_::HyperLagrange<dim>::value);

    // check if shape is a simplex (TRIANGLE|TETRAHEDRON)
    static const bool isSimplex        = (SHAPE == corlib::Simplex<dim>::theShape);

    // check if shape is a hypercube (PIXEL|VOXEL)
    static const bool isTensorSpline   = (SHAPE == corlib::HyperCube<dim>::theShape);

    // hyperquad with tensor-product
    static const bool isTensorLagrange = (isHyperLagrange and isTensor);
    // hyperquad without tensor-product
    static const bool isSerendipity    = (isHyperLagrange and not isTensor);

    // abort for dim==1, as it would not allow to deduce the shape function class
    FTL_STATIC_ASSERT_MSG( (dim != 1), 
                           "Do not use this object for one-D, ambiguous outcome" );
    
    //logical check for tensor-products
    FTL_STATIC_ASSERT_MSG( (isTensor == (isTensorSpline or isTensorLagrange) ),
                           "No proper tensor-product shape function");

    
    // set shape function classifier according to above values
    static const corlib::SfunClass value = ( isSimplex ? SIMPLEX :
                                             ( isTensorSpline ? TENSORSPLINE :
                                               ( isTensorLagrange ? TENSORLAGRANGE :
                                                 ( isSerendipity ? SERENDIPITY : NOSHAPEFUN )
                                                 ) 
                                               ) 
                                             );
};
    
//------------------------------------------------------------------------------
namespace corlib{ 

    //--------------------------------------------------------------------------
    /** \brief Define the type of a shape function using degree and class
     *  \details By using the above traits, this object defines the type of
     *  a shape function if its class, degree and shape dimension are given.
     *  \tparam DEGREE Polynomial degree of the functions
     *  \tparam DIM    Spatial dimension of the underlying shape
     *  \tparam SFC    Class of shape function
     */
    template<unsigned DEGREE, unsigned DIM, corlib::SfunClass SFC> 
    struct ShapeFunType;

    //! \cond SKIPDOX
    template<unsigned DEGREE, unsigned DIM>
    struct ShapeFunType<DEGREE,DIM,SIMPLEX>
    {
        typedef typename 
        corlib::detail_::SimplexShapeFunction<DEGREE,DIM>::type type;
    };

    template<unsigned DEGREE, unsigned DIM>
    struct ShapeFunType<DEGREE,DIM,TENSORSPLINE>
    {
        typedef typename 
        corlib::detail_::TensorSplineShapeFunction<DEGREE,DIM>::type type;
    };

    template<unsigned DEGREE, unsigned DIM>
    struct ShapeFunType<DEGREE,DIM,TENSORLAGRANGE>
    {
        typedef typename 
        corlib::detail_::TensorLagrangeShapeFunction<DEGREE,DIM>::type type;
    };

    template<unsigned DEGREE, unsigned DIM>
    struct ShapeFunType<DEGREE,DIM,SERENDIPITY>
    {
        typedef typename 
        corlib::detail_::SerendipityShapeFunction<DEGREE,DIM>::type type;
    };
    //! \endcond SKIPDOX
}

//------------------------------------------------------------------------------
namespace corlib{

    //--------------------------------------------------------------------------
    /** \brief Define the polynomial degree of a shape function
     *  \details Given the classifier corlib::SfunClass, the number of nodes and
     *  the dimension, this object defines the signficant polynomial degree of
     *  the shape function using above traits structs.
     *  \tparam NNODS  Number of nodes
     *  \tparam DIM    Dimension of the shape
     *  \tparam SFC    Shape function class
     */
    template<unsigned NNODS, unsigned DIM, corlib::SfunClass SFC>
    struct ShapeFunDegree;

    //! \cond SKIPDOX
    template<unsigned NNODS, unsigned DIM>
    struct ShapeFunDegree<NNODS,DIM,SIMPLEX>
    {
        static const unsigned value = 
            corlib::detail_::SimplexShapeFunctionDegree<NNODS,DIM>::value;
    };

    template<unsigned NNODS, unsigned DIM>
    struct ShapeFunDegree<NNODS,DIM,TENSORSPLINE>
    {
        static const unsigned value = 
            corlib::detail_::TensorSplineShapeFunctionDegree<NNODS,DIM>::value;
    };

    template<unsigned NNODS, unsigned DIM>
    struct ShapeFunDegree<NNODS,DIM,TENSORLAGRANGE>
    {
        static const unsigned value = 
            corlib::detail_::TensorLagrangeShapeFunctionDegree<NNODS,DIM>::value;
    };

    template<unsigned NNODS, unsigned DIM>
    struct ShapeFunDegree<NNODS,DIM,SERENDIPITY>
    {
        static const unsigned value = 
            corlib::detail_::SerendipityShapeFunctionDegree<NNODS,DIM>::value;
    };
    //! \endcond SKIPDOX
}


#endif
