// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//
//                 Computational Structural Mechanics Lab
//                         University of Cambridge
//
// This software is copyrighted and all rights are retained by the University of
// Cambridge. All use, disclosure, and/or reproduction of any part not expressly
// authorised by the University of Cambridge is prohibited. (C) 2010
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

//! @file ComputeElementMatrix.hpp

#ifndef corlib_computeelementmatrix_h
#define corlib_computeelementmatrix_h

#include <boost/type_traits.hpp>
#include <boost/function.hpp>
#include <boost/bind.hpp>

#include <Eigen/Core>
#include <corlib/eigenX.hpp>

//------------------------------------------------------------------------------
namespace corlib{
    
    template<typename QUAD, typename INTEGRAND> 
    class ComputeElementMatrix;

}

//------------------------------------------------------------------------------
/** \brief Numerical integration of an element stiffness matrix
 *  \details
 *  \tparam QUAD      Type of quadrature to be used
 *  \tparam INTEGRAND Type of integrand to operate on 
 */
template<typename QUAD, typename INTEGRAND>
class corlib::ComputeElementMatrix
    : public boost::function< void( typename INTEGRAND::arg1_type, 
                                    Eigen::MatrixXd & ) >
{
public:
    //! @name Convenience typedefs
    //@{
    typedef          QUAD                                         Quadrature;
    typedef          INTEGRAND                                    Integrand;
    typedef typename Integrand::arg1_type                         ElementPtr;
    //@}

    //! Cstor with reference to integrand and quadrature objects
    ComputeElementMatrix( Quadrature & quadrature,
                          Integrand  & integrand )
        : quadrature_ ( quadrature ),
          integrand_(   integrand  )
    {}
    
    /** Function call operator applied to element ptr and matrix storage
     *  \attention <i> Storage matrix has to be initialised with zeros</i>
     *  \param[in]     ep      Pointer to element 
     *  \param[in,out] matrix  Result container
     */
    void operator()( ElementPtr ep, Eigen::MatrixXd & matrix  )
    {
        // - Access to quadrature 
        typename Quadrature::QuadIter qiter = quadrature_.begin( );
        typename Quadrature::QuadIter qend  = quadrature_.end( );
        // - Iterate over quadrature (weight,point)-pairs
        for ( ; qiter != qend; qiter ++ ) {
            // obtain quadrature weight and point
            const                      double  weight = qiter -> first;
            const typename Quadrature::VecDim  point  = qiter -> second;
            // call integrand functor 
            integrand_( ep, point, weight, matrix );
        }

        return;
    }

private:
    Quadrature & quadrature_; //!< Functor representing the quadrature rule
    Integrand  & integrand_;  //!< Functor representing the integrand 
};

#endif
