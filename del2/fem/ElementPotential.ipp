// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//
//                 Computational Structural Mechanics Lab
//                         University of Cambridge
//
// This software is copyrighted and all rights are retained by the University of
// Cambridge. All use, disclosure, and/or reproduction of any part not expressly
// authorised by the University of Cambridge is prohibited. (C) 2010
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

//! @file ElementPotential.ipp

#include <functional>
#include <Eigen/Core>
#include <corlib/linalg.hpp>
#include <corlib/EvaluateField.hpp>
#include <corlib/eigenX.hpp>

namespace eigenX = corlib::eigenX;

//------------------------------------------------------------------------------
/** Accessor to the indices of the nodal degrees of freedom
 *  \param[in] dofIndices Vector for the storage of the dof indices
 */
template<typename BELEMENT>
void del2::fem::ElementPotential<BELEMENT>::
getDofIndices( std::vector<unsigned> & dofIndices ) const
{
    // resize the vector (important for assembly!)
    dofIndices.resize( numNodes * dof );
    std::vector<unsigned>::iterator iter = dofIndices.begin();

    typename BasisElement::NodeConstIterator first = BasisElement::nodesBegin();
    typename BasisElement::NodeConstIterator last  = BasisElement::nodesEnd();
    for ( ; first != last; ++first ) {
        std::vector<unsigned> aux;
        (*first) -> copyDofArray( aux );
        iter = std::copy( aux.begin(), aux.end(), iter );
    }

    return;
}

//------------------------------------------------------------------------------
/**  Contribution to the element stiffness matrix in a quadrature rule
 *   The element stiffness matrix for the Laplace operator reads
 *   \f[
 *         K[ m, n ] = \int_\Omega \mu \nabla phi^n \cdot \nabla phi^m dx
 *   \f]
 *   This object adds the weighted integrand evaluated at a local coordinate
 *   \f$\xi\f$ to a provided storage.
 *   \param[in]  xi       Local coordinate: quadrature point
 *   \param[in]  weight   Weight corresponding to the quadrature point
 *   \param[out] result   Result storage
 */
template<typename BELEMENT>
void del2::fem::ElementPotential<BELEMENT>::
stiffnessIntegrand( const VecLDim & xi, const double & weight,
                    Eigen::MatrixXd & result ) const
{
    // debug check
    assert( (result.rows() == (numNodes * dof)) and
            (result.cols() == (numNodes * dof)) );

    // get kinematic quantities: jacobian, Dphi/DX
    MatDimNN dphiDx; 
    const double detG = BasisElement::globalDerivatives( xi, dphiDx );

    // contribute to element stiffness
    result += conductivity_ * detG * weight * 
        (dphiDx.transpose() * dphiDx);

    return;
}

//------------------------------------------------------------------------------
/** Computation of the term 
 *  \f[
 *      T[ 0, m * 1 + 0] =  \mu (\nabla_x \phi^m \cdot n )
 *  \f]
 *  which is part of the surface integral in Nitsche's method.
 *  \param[in]  xi      Local coordinate of evaluation
 *  \param[in]  normal  Outward unit normal vector
 *  \param[out] result  Result storage
 */
template<typename BELEMENT>
void del2::fem::ElementPotential<BELEMENT>::
nitscheTraction( const VecLDim & xi, const VecDim & normal, 
                 eigenX::MatrixSd<dof, numNodes*dof> & result ) const
{
    // shape function derivatives w.r.t global x
    MatDimNN dphiDx; 
    BasisElement::globalDerivatives( xi, dphiDx );
    // dphi / dn
    const VecNN dphiDn = conductivity_ * (dphiDx.transpose() * normal); 

    // store in result array
    result.row(0) = dphiDn;
}

//------------------------------------------------------------------------------
//! Evaluate the field gradient times conductivity at give coordinate #xi
template<typename BELEMENT>
eigenX::VectorSd<3>
del2::fem::ElementPotential<BELEMENT>::internalFlux( const VecLDim & xi ) const
{
    // use external class for field evaluation
    typedef boost::function< VecDof( const Node* ) > GetPot;
    corlib::EvaluateFieldGradient<GetPot,const SelfType_> 
        fluxEval( boost::bind( &Node::getPotential, _1 ) );
    eigenX::MatrixSd <dim, dof> gradP = fluxEval( this, xi );

    // apply material law
    gradP *= conductivity_;  

    // convert to Vec3
    Vec3 flux; flux.setZero();
    flux.head(dim) = gradP.col(0);
    return flux;
}

//------------------------------------------------------------------------------
//! Evaluate the potential at some point #xi inside the element
template<typename BELEMENT>
eigenX::VectorSd<BELEMENT::Node::dof>
del2::fem::ElementPotential<BELEMENT>::potential( const VecLDim & xi ) const
{
    // use external object for the field evaluation
    typedef boost::function< VecDof( const Node* ) > GetPot;
    corlib::EvaluateField<GetPot,const SelfType_> 
        potEval( boost::bind( &Node::getPotential, _1 ) );
    return potEval( this, xi );
}

//------------------------------------------------------------------------------
/** Compute nodal forces due to an applied body force. This function is given
 *  the coordinates and weight of a quadrature point. Moreover, a function object
 *  is passed with the operator of type VecDof = func( VecDim), where the 
 *  argument refers to the reference coordinate of the integration point. 
 *  \tparam FUNC       Type of function object
 *  \param[in] xi      Local coordinate of quadrature point
 *  \param[in] weight  Weight of quadrature point
 *  \param[in] func    The force function object
 *  \param[in] factor  Scalar multiplier
 */
template<typename BELEMENT>
template<typename FUNC>
void del2::fem::ElementPotential<BELEMENT>::
bodyForce( const VecLDim & xi, const double & weight, FUNC func, 
           const double & factor )
{
    // get global coordinate
    VecDim x; BasisElement::geometry( xi, x );

    // evaluate shape functions
    VecNN phi; BasisElement::sfun( xi, phi );

    // compute value of body force function
    const VecDof fofX = func( x );

    // get jacobian of this element
    const double detG = this -> jacobian( xi );

    // Compute the nodal forces 
    const MatDofNN F = weight * detG * factor * (fofX * phi.transpose());

    // pass to nodes
    typename BasisElement::NodeIterator  first = BasisElement::nodesBegin();
    typename BasisElement::NodeIterator  last  = BasisElement::nodesEnd();
    for ( unsigned n = 0; first != last; ++first, n++ ) {
        (*first) -> addToForce( F.col(n) );
    }

    return;
}

//------------------------------------------------------------------------------
template<typename BELEMENT>
template<typename FUNC>
typename BELEMENT::Node::VecDof
del2::fem::ElementPotential<BELEMENT>::
equationResidual( const VecLDim & xi, FUNC func ) const
{
    // force term
    VecDof f; f.clear();

    if ( func ) {
        VecDim x; BasisElement::geometry( xi, x );

        f = func( x );
    }

    // get second derivatives of shape functions
    typename BasisElement::MatVecNNDimDim ddphiDdx;
    BasisElement::globalSecondDerivatives( xi, ddphiDdx );

    // get nodal solutions
    MatDofNN U;
    typename BasisElement::NodeConstIterator  first = BasisElement::nodesBegin();
    typename BasisElement::NodeConstIterator  last  = BasisElement::nodesEnd();
    for ( unsigned n = 0; first != last; ++first, n++ ) {
        U.col(n) = (*first)->getPotential(); 
    }

    // compute residual of Poisson's equation
    VecDof res = - f;

    for ( unsigned n = 0; n < numNodes; n ++ ) {
        // Laplacian applied to nodal shape function
        double laplPhi = 0;
        for ( unsigned d = 0; d < dim; d ++ ) 
            laplPhi += ddphiDdx( d, d )( n );

        // subtract from residual
        res -= laplPhi * U.col(n);
    }
    res *= conductivity_;

    return res;
}

//------------------------------------------------------------------------------
template <typename BELEMENT>
void del2::fem::ElementPotential<BELEMENT>::
stiffnessDerivIntegrand ( const VecLDim & xi, const double & weight, 
                          Eigen::MatrixXd & result ) const
{
    // get kinematic quantities: jacobian, Dphi/DX
    MatDimNN dphiDx;
    const double detG = BasisElement::globalDerivatives( xi, dphiDx );
    
    // contribution to the element stiffness derivative with respect to conductivity
    result += 1.0 * detG * weight
            * ( dphiDx.transpose( ) * dphiDx );
 
    return;
}
