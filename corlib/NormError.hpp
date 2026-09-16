// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//
//                 Computational Structural Mechanics Lab
//                         University of Cambridge
//
// This software is copyrighted and all rights are retained by the University of
// Cambridge. All use, disclosure, and/or reproduction of any part not expressly
// authorised by the University of Cambridge is prohibited. (C) 2010
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

//! @file NormError.hpp

#ifndef corlib_normerror_h
#define corlib_normerror_h
//------------------------------------------------------------------------------
//! boost  includes
#include <boost/function.hpp>
#include <boost/type_traits.hpp>
//! corlib includes
#include <corlib/verify.hpp>


//------------------------------------------------------------------------------
namespace corlib{
    template<typename EXACT, typename APPROX, typename ELEMENT> 
    class NormError;

}

//------------------------------------------------------------------------------
/** \brief Computation of the approximation error using the analytical solution
 *  \details Denote \f$u\f$ as the analytical solution and \f$u^h\f$ the 
 *  corresponding approximation. Then, this functor computes
 *  \f[
 *            e_i(x) = ( D[ u(x(\xi)) ] - D[ u^h(x(\xi)) ] )_i^2
 *  \f]
 *  where \f$\xi\f$ is a point in reference coordinates and \f$x\f$ the
 *  corresponding physical coordinate provided by the element's geometry 
 *  function. Normally, \f$D\f$ is either identity or the gradient operation.
 *  The analytical solution is provided by a functor of type 'Result (x)'
 *  and the approximation by a functor of type 'Result (Cell*, xi)'.
 *  The error result \f$e(x)\f$ is weighted by some integration weight and
 *  added to a reference to some storage variable.
 *  If called by an integrator, this functor can be used to compute for
 *  instance the error in the \f$L_2\f$-norm or \f$H^1\f$-seminorm.
 *  \tparam EXACT    Functor for analytical solution
 *  \tparam APPROX   Functor for the approximate solution
 *  \tparam ELEMENT  Type of object to integrate over (can be inferred from APPROX)
 */
template<typename EXACT, typename APPROX, 
         typename ELEMENT = 
         typename boost::remove_pointer<typename APPROX::first_argument_type>::type>
class corlib::NormError 
    : public boost::function<void (ELEMENT *, 
                                   const typename ELEMENT::VecLDim &, 
                                   const double & ) >
{
public:
    //! @name Template parameters
    //@{
    typedef EXACT                                         ExactSolution;
    typedef APPROX                                        ApproximateSolution;
    typedef ELEMENT                                       Element;
    //@}

    //! @name Derived types
    //@{
    typedef typename ExactSolution::result_type           ResultType;
    typedef typename Element::VecLDim                     VecLDim;
    typedef typename Element::VecDim                      VecDim;
    //!}
    
    //! Sanity check
    FTL_STATIC_ASSERT_MSG( (boost::is_same<typename ApproximateSolution::result_type,
                                           ResultType>::value),
                           "Approximate and exact solution functors do not match" );

    //--------------------------------------------------------------------------
    /** Cstor with both functors and the reference to the storage variable
     *  \param[in] exactSolution       Functor for the exact solution at a point
     *  \param[in] approximateSolution Functor for the element approximation
     *  \param[in] errorSquared        Reference to the square of the error
     */
    NormError( ExactSolution       exactSolution, 
               ApproximateSolution approximateSolution,
               ResultType        & errorSquared )
        : exactSolution_(       exactSolution ), 
          approximateSolution_( approximateSolution ),
          errorSquared_(        errorSquared ) { }

    //--------------------------------------------------------------------------
    /** Function call on cell pointer given the integration point and weight
     *  \param[in] ep        Pointer to element
     *  \param[in] quadPoint Quadrature point
     *  \param[in] weight    Quadrature weight
     */
    void operator()( const Element * ep, const VecLDim & quadPoint, 
                     const double & weight )
    {
        //! - get global point
        VecDim x; ep -> geometry( quadPoint, x );

        //! - get exact and apprximate solutions, compute the difference
        ResultType analytical  = exactSolution_( x );
        ResultType approximate = approximateSolution_( ep, quadPoint );
        ResultType difference  = analytical - approximate;

        //! - get jacobian from element
        const double jac = ep -> jacobian( quadPoint );

        //! - store the squared error components
        errorSquared_ += 
            difference.cwiseProduct( difference ) * jac * weight;

        return;
    }

private:
    //--------------------------------------------------------------------------
    ExactSolution        exactSolution_;       //!< Exact solution functor
    ApproximateSolution  approximateSolution_; //!< Approximate solution functor
    ResultType         & errorSquared_;        //!< Reference to result
};
//------------------------------------------------------------------------------
#endif
