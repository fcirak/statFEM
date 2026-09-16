// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//
//                 Computational Structural Mechanics Lab
//                         University of Cambridge
//
// This software is copyrighted and all rights are retained by the University of
// Cambridge. All use, disclosure, and/or reproduction of any part not expressly
// authorised by the University of Cambridge is prohibited. (C) 2010
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

//! @file BSpline.hpp

#ifndef corlib_bspline_h
#define corlib_bspline_h
//------------------------------------------------------------------------------
#include <Eigen/Core>

#include <corlib/eigenX.hpp>
#include <corlib/misc.hpp>

//------------------------------------------------------------------------------
namespace corlib{
    template<unsigned DEGREE> class BSpline;
    namespace eigenX = corlib::eigenX;
}

//------------------------------------------------------------------------------
/** \brief   Evaluation of B-Spline FE shape functions and derivatives
 *  \details Using [0,1] as the reference interval, n+1 B-Splines of degree n
 *  (cardinal, uniform splines ) are non-zero on this interval. The splines are
 *  numbered from left to right. This object provides the functions to evaluate
 *  the spline functions and their first and second derivatives on [0,1].
 *  \tparam DEGREE Polynomial degree of the B-Splines                       */
template<unsigned DEGREE> 
class corlib::BSpline
{
public:
    static const unsigned degree       = DEGREE;
    static const unsigned numSplines   = degree+1;
    
    typedef eigenX::VectorSd<numSplines>         VecNspl;

public:
    /** Evaluation of B-Spline FE shape functions of order DEGREE
     *  \param[in]  xi  Local coordinate from [0,1]
     *  \param[out] phi Result \f$\phi_i(\xi)\f$      */
    void evaluate( const double & xi, VecNspl & phi ) const;

    /** Evaluate the first derivative of the B-Spline FE functions
     *  \param[in]  xi   Local coordinate from [0,1]
     *  \param[out] dphi Result \f$ \phi^\prime_i(\xi) \f$     */
    void evaluateGradient( const double & xi, VecNspl & dphi ) const;

    /** Evaluate the second derivative of the B-Spline FE functions
     *  \param[in]  xi    Local coordinate from [0,1]
     *  \param[out] ddphi Result \f$ \phi^{\prime\prime}_i(\xi) \f$ */
    void evaluateHessian( const double & xi, VecNspl & ddphi ) const;
};


#include "BSpline.ipp"
//------------------------------------------------------------------------------
#endif
