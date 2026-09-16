// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//
//                 Computational Structural Mechanics Lab
//                         University of Cambridge
//
// This software is copyrighted and all rights are retained by the University of
// Cambridge. All use, disclosure, and/or reproduction of any part not expressly
// authorised by the University of Cambridge is prohibited. (C) 2010
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

//! @file Shapefun.hpp

#ifndef corlib_shapefun_h
#define corlib_shapefun_h
//------------------------------------------------------------------------------
#include <Eigen/Core>
#include <corlib/Shape.hpp>
#include <corlib/ShapefunTraits.hpp>
#include <corlib/eigenX.hpp>

namespace corlib{
    template< corlib::shape SHAPE, unsigned NFUN > class Shapefun;

    //! \cond SKIPDOX
    namespace detail_{
        
        // use corlib::ShapeFunClassifier for defining the shape fun class
        template<corlib::shape SHAPE, unsigned NFUN, unsigned DIM>
        struct ShapeFunClassTraits
        {
            static const corlib::SfunClass value =
                corlib::ShapeFunClassifier<SHAPE,NFUN>::value;
        };

        // default to simplex shape function for dim = 1
        template<corlib::shape SHAPE, unsigned NFUN>
        struct ShapeFunClassTraits<SHAPE,NFUN,1>
        {
            static const corlib::SfunClass value = corlib::SIMPLEX;
        };

        // default to simplex shape function for dim = 0
        template<corlib::shape SHAPE, unsigned NFUN>
        struct ShapeFunClassTraits<SHAPE,NFUN,0>
        {
            static const corlib::SfunClass value = corlib::SIMPLEX;
        };

        //----------------------------------------------------------------------
        // Workaround to catch the DIM=0 case
        template<unsigned NFUN, corlib::SfunClass SFC, unsigned DIM>
        struct ShapeFunDegree
        {
            static const unsigned value = 
                corlib::ShapeFunDegree<NFUN,DIM,SFC>::value;
        };

        template<unsigned NFUN, corlib::SfunClass SFC>
        struct ShapeFunDegree<NFUN,SFC,0>
        {
            static const unsigned value = 0; 
        };

    }
    //! \endcond SKIPDOX
}

//------------------------------------------------------------------------------
/** \brief Finite element shape functions
 *  \details Evaluates the FE shape functions and their gradient
 *  \tparam SHAPE shape of the element (as defined in shape enumerator)
 *  \tparam NFUN  number of shape functions (usually equal to number of nodes)
 */
template< corlib::shape SHAPE, unsigned NFUN >
class corlib::Shapefun
{
public:
    static const corlib::shape myShape      = SHAPE;
    static const unsigned      numFunctions = NFUN;
    static const unsigned      dim          = corlib::ShapeTraits<SHAPE>::dim;
    // Classification of this shape function object 
    // {SIMPLEX|TENSORLAGRANGE|SERENDIPITY}
    static const corlib::SfunClass sfc      = 
        detail_::ShapeFunClassTraits<myShape,numFunctions,dim>::value;
    // (Signficant) polynomial degree of this shape function object 
    static const unsigned      degree       = 
        detail_::ShapeFunDegree<numFunctions,sfc,dim>::value;

    typedef eigenX::VectorSd<dim>                  VecDim;
    typedef eigenX::VectorSd<numFunctions>         VecNFun;
    typedef eigenX::MatrixSd<dim,numFunctions>     MatDimNFun;
    typedef Eigen::Matrix<VecNFun,dim,dim>         MatVecNFunDimDim;

public:
 
    /** @name Parameterized functions */
    //@{
    //------------------------------------------------------------------------------
    /** \brief \f$ \varphi^k(\xi) \f$ \details
     * Evaluation of the shape functions at given local coordinate xi
     *  \param[in]  xi  Local coordinate at which the functions are evaluated
     *  \param[out] phi Vector with result: phi[k] = \f$ \varphi^k(\xi) \f$
     */
    void evaluate( const VecDim & xi, VecNFun & phi ) const;

    //------------------------------------------------------------------------------
    /** \brief \f$ d \varphi^k / d \xi_i \f$ \details
     * Evaluation of the shape functions derivatives at given local coordinate xi
     *  \param[in]  xi   Local coordinate at which the derivatives are evaluated
     *  \param[out] dphi Result matrix: dphi[i,k] = \f$ (d\varphi^k/ d\xi_i)(\xi) \f$
     */
    void evaluateGradient( const VecDim & xi, MatDimNFun & dphi ) const;

    //------------------------------------------------------------------------------
    /** \brief \f$ d^2 \varphi^k / ( d \xi_i  d \xi_j ) \f$ \details
     * Evaluation of the shape functions second derivatives at given local coordinate xi
     *  \param[in]  xi    Local coordinate at which the derivatives are evaluated
     *  \param[out] ddphi Result matrix: ddphi(i,j)(k) = \f$ (d^2\varphi^k/(d\xi_i d\xi_j))(\xi) \f$
     */
    void evaluateHessian( const VecDim & xi, MatVecNFunDimDim & ddphi ) const;

    //------------------------------------------------------------------------------
    /** \brief Interpolation points \details
     * Return the location of the interpolation points in local coordinates
     *  \param[out] ipoints Interpolation points of this shape functions
     */
    static void giveInterpolationPoints( MatDimNFun & ipoints );
    //@}

};

#include "Shapefun.ipp"

//------------------------------------------------------------------------------
#endif 
