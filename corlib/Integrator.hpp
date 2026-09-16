// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//
//                 Computational Structural Mechanics Lab
//                         University of Cambridge
//
// This software is copyrighted and all rights are retained by the University of
// Cambridge. All use, disclosure, and/or reproduction of any part not expressly
// authorised by the University of Cambridge is prohibited. (C) 2010
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

//! @file Integrator.hpp

#ifndef corlib_integrator_h
#define corlib_integrator_h

#include <corlib/verify.hpp>

#include <boost/numeric/ublas/io.hpp>
#include <boost/static_assert.hpp>
#include <boost/type_traits.hpp>
#include <boost/function.hpp>

//------------------------------------------------------------------------------
namespace corlib{
    template<typename QUAD, typename INTEGRAND, typename ELEMENT, bool DUMMY> 
    class Integrator;
}

//------------------------------------------------------------------------------
/** \brief Generic numerical integration function object
 *  \details Functor to numerically integrate over an element by using a 
 *  Quadrature rule and letting the element evaluate its integrand function
 *  at the passed point.
 *  \tparam QUAD      Type of quadrature
 *  \tparam INTEGRAND Type of integrand function (functor)
 *  \tparam ELEMENT   Type of element to integrate over
 */
template<typename QUAD, typename INTEGRAND, 
         typename ELEMENT = typename boost::remove_pointer<typename INTEGRAND::arg1_type>::type,
         bool DUMMY = false>
class corlib::Integrator 
    : public boost::function< void( ELEMENT * ) >
{
    FTL_STATIC_ASSERT_MSG( ( (QUAD::theShape == ELEMENT::myShape) or DUMMY), 
                           "Geometric shape of quadrature and element must coincide" );
public:
    typedef QUAD      Quadrature;
    typedef INTEGRAND Integrand;

    //! Cstor with reference to integrand and quadrature objects
    Integrator( Quadrature & quadrature,
                Integrand  & integrand )
        : quadrature_( quadrature ),
          integrand_(  integrand )
    {}

   //--------------------------------------------------------------------------
   /** Perform the numerical integration: instantiate a quadrature object and 
    *  perform the quadrature.
    *  \param[in] ep Pointer to an element
    */
    void operator()( ELEMENT*  ep ) 
    {
        typename Quadrature::QuadIter qiter = quadrature_.begin( );
        typename Quadrature::QuadIter qend  = quadrature_.end( );
        
        for ( ; qiter != qend; qiter ++ ) {
            // obtain quadrature weight and point
            const                      double  weight = qiter -> first;
            const typename Quadrature::VecDim  point  = qiter -> second;
            // call integrand functor 
            integrand_( ep, point, weight);
        }
    }

protected:
    Quadrature & quadrature_; //!< Functor representing the quadrature rule
    Integrand  & integrand_;  //!< Functor representing the integrand
};
//------------------------------------------------------------------------------
#endif
