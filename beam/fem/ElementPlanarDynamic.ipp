// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//
//                 Computational Structural Mechanics Lab
//                         University of Cambridge
//
// This software is copyrighted and all rights are retained by the University of
// Cambridge. All use, disclosure, and/or reproduction of any part not expressly
// authorised by the University of Cambridge is prohibited. (C) 2010
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

//! @file ElementPlanarDynamic.ipp

//------------------------------------------------------------------------------
/** Evaluate the integral kernel function for the mass matrix at the given local
 *  coordinate xi and sum the result into a matrix.
 *  \f[
 *      M[m d +i, n d +i] = \int \rho \varphi^m \varphi^n \det G dX
 *  \f]
 *  where \f$ d \f$ refers to the number of dofs per node.
 *  \param[in] xi      Local coordinate at which this function is evaluated
 *  \param[in] weight  The corresponding quadrature weight
 */
template< typename NODE, typename SFUN, typename MAT >
void beam::fem::ElementPlanarDynamic<NODE,SFUN,MAT>::massIntegrand( const VecLDim& xi, 
                                                                    const double& weight )
{
    // shape functions
    VecNF shpFct;  this -> sFun_( xi, shpFct );

    // material mass density
    const double rho = material_ -> density();

    // get jacobian from element
    const double det = this -> jac_( xi, false /* in reference config */ );

    // compute element mass matrix
    //
    // HINT: Only the translational inertia is computed.
    //       Including the rotational inertia leads to a non-linear
    //       inertia force. In this case, the inertia forces are not simply
    //       mass matrix times nodal acceleration. 
    for ( unsigned iDof = 0; iDof < dof; ++iDof ) {
        for ( unsigned kNod = 0; kNod < numNodesSN; ++kNod ) {
            const unsigned ikDofNod = kNod*dof + iDof;

            for ( unsigned jDof = 0; jDof < dof; ++jDof ) {
                for ( unsigned lNod = 0; lNod < numNodesSN; ++lNod ) {
                    const unsigned jlDofNod = lNod*dof + jDof;

                    if ( jDof == iDof )
                        elemMass_(ikDofNod, jlDofNod) += 
                            shpFct[ kNod ] * shpFct[ lNod ] * rho * area_ * det * weight;

                }
            }

        }
    }

    return;
}

//------------------------------------------------------------------------------
/** Return nodal velocities
 *
 *  \param[out]  vel  Matrix with the nodal velocities as columns
 */
template< typename NODE, typename SFUN, typename MAT >
void beam::fem::ElementPlanarDynamic<NODE,SFUN,MAT>::nodalVelocities( MatDofNF & vel ) const
{
    this -> nodalQuantity_( vel, std::bind( &NODE::giveVelocities, std::placeholders::_1 ) );
    return;
}

//------------------------------------------------------------------------------
/** \param[out] acc Matrix with the nodal accelerations as columns              */
template< typename NODE, typename SFUN, typename MAT >
void beam::fem::ElementPlanarDynamic<NODE,SFUN,MAT>::nodalAccelerations( MatDofNF & acc ) const
{
    this -> nodalQuantity_( acc, std::bind( &NODE::giveAccelerations, std::placeholders::_1 ) );
    return;
}

//------------------------------------------------------------------------------
template< typename NODE, typename SFUN, typename MAT >
typename NODE::VecDof
beam::fem::ElementPlanarDynamic<NODE,SFUN,MAT>::giveVelocity( const VecLDim& xi ) const
{
    // nodal coordinates
    MatDofNF velN;  this -> nodalVelocities( velN );
    // shape function derivatives
    VecNF phi;      this -> sFun_( xi, phi );
    // return interpolated result
    return velN * phi;
}

