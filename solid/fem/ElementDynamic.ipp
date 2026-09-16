// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//
//                 Computational Structural Mechanics Lab
//                         University of Cambridge
//
// This software is copyrighted and all rights are retained by the University of
// Cambridge. All use, disclosure, and/or reproduction of any part not expressly
// authorised by the University of Cambridge is prohibited. (C) 2010
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

//! @file ElementDynamic.ipp

//------------------------------------------------------------------------------
/** \param[out] V Matrix with the nodal velocities as columns                 */
template<typename BELEMENT, typename MAT>
void solid::fem::ElementDynamic<BELEMENT,MAT>::nodalVelocities( MatDofNN & V ) const
{
    this -> nodalQuantity_( V, boost::bind( &Node::giveVelocities, _1 ) );
    return;
}

//------------------------------------------------------------------------------
/** \param[out] A Matrix with the nodal accelerations as columns              */
template<typename BELEMENT, typename MAT>
void solid::fem::ElementDynamic<BELEMENT,MAT>::nodalAccelerations( MatDofNN & A ) const
{
    this -> nodalQuantity_( A, boost::bind( &Node::giveAccelerations, _1 ) );
    return;
}

//------------------------------------------------------------------------------
/** Evaluate the integral kernel function for the mass matrix at the given local
 *  coordinate xi and sum the result into a matrix.
 *  \f[
 *      M[m d + i, n d + i ] = \int \rho \varphi^m \varphi^n \det G d X
 *  \f]
 *  where \f$ d \f$ refers to the number of dofs per node.
 *  \param[in]  xi      Local coordinate at which this function is evaluated
 *  \param[in]  weight  The corresponding quadrature weight
 *  \param[out] result  Result storage
 */
template<typename BELEMENT, typename MAT>
void solid::fem::ElementDynamic<BELEMENT,MAT>::
massIntegrand( const VecLDim & xi, const double  & weight, 
               ublas::matrix<double> & result ) const
{
    // debug check
    assert( (result.size1() == numNodes * dof) and
            (result.size2() == numNodes * dof) );

    // material mass density
    const double rho = ElementStatic::material_ -> density( );

    // get jacobian from element
    const double detG = BasisElement::jacobian( xi );

    // evalute the shape functions
    typename BasisElement::VecNN phi; 
    BasisElement::sfun( xi, phi );

    // identity
    ublas::identity_matrix<double> delta_ij(dof);

    for ( unsigned k = 0; k < numNodes; k ++ ) {
        for( unsigned el = 0; el < numNodes; el ++ ) {
            // Block: M[k,el] = rho * phi(k) * phi(el) * delta_ij
            ublas::subrange( result, k*dof, (k+1)*dof, el*dof, (el+1)*dof) +=
                rho * phi(k) * phi(el) * detG * weight * delta_ij;
        }
    }
    
    return;
}

