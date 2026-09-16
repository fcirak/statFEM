// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//
//                 Computational Structural Mechanics Lab
//                         University of Cambridge
//
// This software is copyrighted and all rights are retained by the University of
// Cambridge. All use, disclosure, and/or reproduction of any part not expressly
// authorised by the University of Cambridge is prohibited. (C) 2010 --
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

//! @file ElementMatern.ipp

#include <functional>
#include <corlib/verify.hpp>

//------------------------------------------------------------------------------
/** Accessor to the indices of the nodal degrees of freedom
 *  \param[in] dofIndices Vector for the storage of the dof indices
 */
template<typename BELEMENT>
void del2::fem::ElementMatern<BELEMENT>::
getDofIndices( std::vector<unsigned> & dofIndices ) const
{
    // resize the vector (important for assembly!)
    dofIndices.resize( numNodes * dof );
    std::vector<unsigned>::iterator iter = dofIndices.begin( );

    typename BasisElement::NodeConstIterator first = BasisElement::nodesBegin( );
    typename BasisElement::NodeConstIterator last  = BasisElement::nodesEnd( );
    for ( ; first != last; ++ first ) {
        std::vector<unsigned> aux;
        ( *first )->copyDofArray( aux );
        iter = std::copy( aux.begin( ), aux.end( ), iter );
    }

    return;
}

//------------------------------------------------------------------------------
/**  Contribution to the element Matern matrix in a quadrature rule.
 *   The element matrix for the Matern operator reads
 *   \f[
 *         K[ m, n ] = \int_\Omega ( H \nabla \phi^n ) \cdot \nabla \phi^m dx
 *                     + \int_\Omega \kappa^2 \phi^n \phi^m dx
 *   \f]
 *   For the fractional exponent case, either lumped or consistent mass matrix
 *   should be applied throughout the algorithm for accurate approximation. To
 *   ensure that the precision matrix is sparse, therefore, the mass integral is
 *   lumped to the diagonal entries of the element matrix. This object adds the
 *   weighted integrand evaluated at a local coordinate \f$\xi\f$ to a provided
 *   storage.
 *   \param[in]  xi       Local coordinate: quadrature point
 *   \param[in]  weight   Weight corresponding to the quadrature point
 *   \param[out] result   Result storage
 */
// todo Generalise the integrand for the anistropic case on a manifold.
template<typename BELEMENT>
void del2::fem::ElementMatern<BELEMENT>::
maternIntegrand( const VecLDim & xi, const double & weight,
                 Eigen::MatrixXd & result ) const
{
    // debug check
    FTL_VERIFY( ( result.rows( ) == ( numNodes * dof ) ) and
                ( result.cols( ) == ( numNodes * dof ) ) );

    // get kinematic quantities: jacobian, Dphi/DX
    MatDimNN dphiDx;
    const double detJ = BasisElement::globalDerivatives( xi, dphiDx );

    // contribute Laplacian component to element matrix
    result += detJ * weight
            * ( ( diffusivity_ * dphiDx ).transpose( ) * dphiDx );

    // evaluate shape functions
    VecNN phi;
    BasisElement::sfun( xi, phi );

    // contribute lumped mass component to element matrix
    result.diagonal( ) += detJ * weight * ( kappa_ * kappa_ * phi );

    return;
}

//--------------------------------------------------------------------------
/** Contribution of unit lumped mass to nodal force vector in a quadrature
 *  rule. The element unit lumped mass vector reads
 *  \f[
 *        f[ m ] = \int_\Omega \phi^m dx
 *  \f]
 *  This object adds the weighted integrand evaluated at a local coordinate
 *  \f$\xi\f$ and stores in the force vector to avoid defining additional
 *  storage.
 *  \param[in]  xi      Local coordinate: quadrature point
 *  \param[in]  weight  Weight corresponding to the quadrature point
 */
template<typename BELEMENT>
void del2::fem::ElementMatern<BELEMENT>::
lumpedMassIntegrand( const VecLDim & xi, const double & weight )
{
    // get Jacobian
    const double detJ = BasisElement::jacobian( xi );

    // evaluate shape functions
    VecNN phi;
    BasisElement::sfun( xi, phi );

    // compute unit lumped mass
    const MatDofNN massMat = weight * detJ * ( phi.transpose( ) );

    // pass to nodes
    typename BasisElement::NodeIterator first = BasisElement::nodesBegin( );
    typename BasisElement::NodeIterator last  = BasisElement::nodesEnd( );
    for ( unsigned n = 0; first != last; ++ first, n ++ )
    {
        ( *first )->addToForce( massMat.col( n ) );
    }

    return;
}
