// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//
//                 Computational Structural Mechanics Lab
//                         University of Cambridge
//
// This software is copyrighted and all rights are retained by the University of
// Cambridge. All use, disclosure, and/or reproduction of any part not expressly
// authorised by the University of Cambridge is prohibited. (C) 2010
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

//! @file   EvaluateField.hpp

#ifndef corlib_evaluatefield_h
#define corlib_evaluatefield_h

#include <array>
#include <corlib/verify.hpp>
#include <corlib/Accessor.hpp>
#include <corlib/eigenX.hpp>

namespace corlib{

    template<typename GIVDAT, typename ELEMENT> 
    class EvaluateField;

    template<typename GIVDAT, typename ELEMENT> 
    class EvaluateFieldGradient;

    namespace detail_{
        template<typename GIVDAT, typename OUT> class GetNodalDatum;
    }

}

//------------------------------------------------------------------------------
/** \brief Functor for the evaluation of an approximated field
 *  \details Assume some abstract finite element approximation
 *  \f[
 *       u_h(\xi) = \sum_n \phi^n(\xi) u^n
 *  \f]
 *  where \f$ u \f$ denotes the field, \f$ \phi \f$ the element's shape 
 *  functions, \f$ u^n \f$ the nodal coefficients and the sum is carried out
 *  over all nodes. This object will evaluate this expression at a given
 *  local coordinate \f$ \xi \f$ and return the result.
 *
 *  Note that this writing does not require classical Lagrange
 *  interpolation but can also be applied if the coefficients of the 
 *  approximation are stored some node-like structure.
 *  This object relies on the function 'iterateOverNodes' that has to be 
 *  available in the element implementation.
 *  \tparam GIVDAT  Type of accessor to get the nodal coefficients
 *  \tparam ELEMENT Type of element to act on
 */
template<typename GIVDAT, typename ELEMENT> 
class corlib::EvaluateField
    : public boost::function<typename GIVDAT::result_type(ELEMENT *,
                                               const typename ELEMENT::VecLDim)>
{
public:
    typedef typename GIVDAT::result_type   result_type;
    typedef typename GIVDAT::argument_type NodePtr;
    static const unsigned numNodes = ELEMENT::numNodes;
    typedef eigenX::VectorSd<numNodes> VecNN;
    typedef boost::function<void( const ELEMENT*, const typename ELEMENT::VecLDim&,
                                  VecNN &) > SFun;

    /** Cstor with node accessor and shape function functor, 
     *  \param[in] op    Node accessor to get the nodal datum
     *  \param[in] sFun  Functor for the evaluation of the shape functions
     */
    EvaluateField( GIVDAT op, 
                   SFun sFun = &ELEMENT::sfun ) 
        : op_(   op ),
          sFun_( sFun ) { }

    /** Function call operator applied to element with local coordinate
     *  \param[in] ep Pointer to element
     *  \param[in] xi Local coordinate of evaluation
     */
    result_type operator()( ELEMENT * ep,
                            const typename ELEMENT::VecLDim & xi )
    {
        // evaluate element's shape functions
        VecNN phi;
        sFun_( ep, xi, phi );

        // set up some storage for the nodal data
        typedef std::array<result_type, numNodes> Storage;
        Storage storage;

        // Accessor functor with iterator to the storage
        detail_::GetNodalDatum<GIVDAT, typename Storage::iterator> 
            gnd( op_, storage.begin() );

        // get nodal data from element
        ep -> iterateOverNodes( gnd );
        
        // compute result as linear combination
        result_type result = corlib::detail_::AccessorTraits<result_type>::zero();
        for ( unsigned n = 0; n < numNodes; n ++ ) {
            result += 
                static_cast<typename corlib::detail_::AccessorTraits<result_type>::ValueType>(phi[n]) 
                * storage[n];
        }
        
        // pass back the result
        return result;
    }

private:
    GIVDAT   op_;   //! Accessor acting on the element's nodes
    SFun sFun_; //! Shape function evaluation functor
};

//------------------------------------------------------------------------------
/** \brief Functor for the evaluation of a field's gradient at some local point
 *  \details Let a field \f$ u \f$ have the approximation as in 
 *  corlib::EvaluateField, then its gradient is written as
 *  \f[
 *        \nabla_x u_h = \sum_n \nabla_x \phi^n(\xi) \otimes u^n
 *  \f]
 *  where all derivatives are taken with respect to the physical coordinate 
 *  \f$ x \f$.
 *  \tparam GIVDAT  Type of accessor to get the nodal coefficients
 *  \tparam ELEMENT Type of element to act on
 */
template<typename GIVDAT, typename ELEMENT> 
class corlib::EvaluateFieldGradient
    : public boost::function<eigenX::MatrixSd<ELEMENT::dim,
                                              GIVDAT::result_type::RowsAtCompileTime>(ELEMENT *,
                                              const typename ELEMENT::VecLDim)>
{
public:
    typedef typename GIVDAT::result_type   datum_type;
    typedef typename GIVDAT::argument_type NodePtr;
    static const unsigned numNodes = ELEMENT::numNodes;
    static const unsigned dof      = datum_type::RowsAtCompileTime;
    static const unsigned dim      = ELEMENT::dim;

    typedef eigenX::MatrixSd<dim,dof>                  MatDimDof;
    typedef eigenX::MatrixSd<dim,numNodes>             MatDimNN;

    typedef boost::function<void( const ELEMENT*, const typename ELEMENT::VecLDim&,
                                  MatDimNN &) > SFunGradX;

    /** Cstor with node accessor and shape function functor, 
     *  \param[in] op        Node accessor to get the nodal datum
     *  \param[in] sFunGrad  Functor for the evaluation of the shape 
     *                       derivatives
     */
    EvaluateFieldGradient( GIVDAT op, 
                           SFunGradX sFunGrad = &ELEMENT::globalDerivatives ) 
        : op_(       op ),
          sFunGrad_( sFunGrad ) { }

    /** Function call operator applied to element with local coordinate
     *  \param[in] ep Pointer to element
     *  \param[in] xi Local coordinate of evaluation
     */
    MatDimDof operator()( ELEMENT * ep,
                          const typename ELEMENT::VecLDim & xi )
    {
        // evaluate element's shape functions
        MatDimNN dPhiDx;
        sFunGrad_( ep, xi, dPhiDx );

        // set up some storage for the nodal data
        typedef std::array<datum_type, numNodes> Storage;
        Storage storage;

        // Accessor functor with iterator to the storage
        detail_::GetNodalDatum<GIVDAT, typename Storage::iterator> 
            gnd( op_, storage.begin() );

        // get nodal data from element
        ep -> iterateOverNodes( gnd ); 
        
        // compute result as linear combination
        MatDimDof result; result.setZero();
        for ( unsigned n = 0; n < numNodes; n ++ ) {
            result += ( dPhiDx.col( n ) * storage[n].transpose() );
        }
        
        // pass back the result
        return result;
    }

private:
    GIVDAT    op_;       //! Accessor acting on the element's nodes
    SFunGradX sFunGrad_; //! Shape function evaluation functor
};

//------------------------------------------------------------------------------
/** \brief Functor which allows to externally store the result of some operator
 *  \details When evaluating an approximated field, it is required to access
 *  a datum stored in the element's nodes. This functor applies an accessor 
 *  operator to a node and store the result in a de-referenced iterator.
 *  \tparam GIVDAT  Type of accessor applied to the nodes
 *  \tparam OUT Type of output iterator for the external storage
 */
template<typename GIVDAT, typename OUT>
class corlib::detail_::GetNodalDatum
    : public boost::function<void( typename GIVDAT::argument_type, GIVDAT &)>
{
public:
    /** Cstor with storage iterator, 
     *  \param[in] out Storage iterator
     *  \param[in] op  Access operator
     */
    GetNodalDatum( GIVDAT op, OUT out ) : op_( op ), out_( out ) {  }

    /** Function call operator to apply accessor to node
     *   \param[in] np Node pointer to act on
     */
    void operator()( typename GIVDAT::argument_type np )
    {
        *out_++ = op_( np );
        return;
    }

private:
    GIVDAT op_;  //!< Operator to access the nodal datum
    OUT    out_; //!< Iterator to some external storage
};

#endif
