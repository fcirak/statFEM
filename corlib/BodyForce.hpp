// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//
//                 Computational Structural Mechanics Lab
//                         University of Cambridge
//
// This software is copyrighted and all rights are retained by the University of
// Cambridge. All use, disclosure, and/or reproduction of any part not expressly
// authorised by the University of Cambridge is prohibited. (C) 2010
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

//! @file BodyForce.hpp

#ifndef corlib_bodyforce_h
#define corlib_bodyforce_h

#include <boost/bind.hpp>
#include <boost/function.hpp>

//------------------------------------------------------------------------------
namespace corlib{
    template<typename ELEMENT, typename FUNC> 
    class BodyForce;
}

//------------------------------------------------------------------------------
/** \brief Integrand wrapper for a prescribed body force / source terms
 *  \details The nodal forces due to applied volume forces normally is of the
 *  form
 *  \f[
 *         F[i] = \int_\Omega f(x) \phi_i(x) dx
 *  \f]
 *  where \f$\phi_i\f$ refers to the i-th nodal shape function and \f$f(x)\f$ is
 *  the prescribed force as a function of the global coordinate. The function
 *  for the force term is passed by the user in form of some functor which is
 *  passed to the element level function stored in #bfEval_. 
 *  \tparam ELEMENT  Type of element 
 *  \tparam FUNC     Type of body force functor f(x)
 */
template<typename ELEMENT, typename FUNC>
class corlib::BodyForce
    : public boost::function<void (ELEMENT*, const typename ELEMENT::VecLDim &,
                                   const double &)>
{
public:
    typedef ELEMENT             Element;
    static const unsigned dim = Element::dim;

private:
    typedef typename Element::VecLDim                           VecLDim_;
    typedef boost::function< void( Element *, const VecLDim_ &, 
                                   const double &, 
                                   FUNC, const double & ) >     BodyForceEvalFun;

public:
    //! Constructor given the function object representing the body force and
    //! (if wanted) a functor representing the element's evaluation function
    BodyForce( FUNC func , 
               BodyForceEvalFun bfEval = 
               boost::bind( &Element::template bodyForce<FUNC>, _1, _2, _3, _4, _5 )
        ) 
        : func_( func ), 
          bfEval_( bfEval ),
          factor_( 1. ) {}

    //! Basic operation for called by a quadrature 
    void operator()(ELEMENT * ep, const VecLDim_ & xi, const double & weight )
    {
        // let element compute forces by passing the force function
        bfEval_( ep, xi, weight, func_, factor_ );
        return;
    }

public:    
    void setFactor( const double & factor ) { factor_ = factor; }

private:
    FUNC func_;       //!< Function object of type vec = f(vec)

    //! Local function representing the integrand of the linear form
    BodyForceEvalFun bfEval_;

    double factor_;   //!< Multiplier for the forces
};

//------------------------------------------------------------------------------
#endif
