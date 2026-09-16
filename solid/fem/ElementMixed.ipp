// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//
//                 Computational Structural Mechanics Lab
//                         University of Cambridge
//
// This software is copyrighted and all rights are retained by the University of
// Cambridge. All use, disclosure, and/or reproduction of any part not expressly
// authorised by the University of Cambridge is prohibited. (C) 2010
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

//! @file ElementMixed.ipp

//------------------------------------------------------------------------------
//! so, far only for original Stokes (-Delta u) + Grad p = f
template<typename NODE, typename SFUNU, typename SFUNP, typename MAT>
void solid::fem::ElementMixed<NODE,SFUNU,SFUNP,MAT>::
stiffnessIntegrand( const VecLDim & xi, const double & weight )
{
    // get kinematic quantities: jacobian, Dphi/DX, F
    MatDimNN   Dphi; 
    Mat3x3     F;
    const double detG = this -> kinematics_( xi, Dphi, F );

    // compute second derivative of strain energy ( d^2 W(F) ) / ( dF dF )
    //Mat3x3x3x3_ Cmat( material_ -> elasticityTensor( F ) );

    // bilinear form with Grad u * Grad v
    for (unsigned m = 0; m < numNodes; m++ ) {
        for (unsigned n = 0; n < numNodes; n++ ) {
            double Kentry = 0;
            for ( unsigned j = 0; j < dof; j ++ ) {
                Kentry += 0.001 * Dphi( j, m ) * Dphi( j, n ) * detG;
            }
            for (unsigned i = 0; i < dof; i++ ) {
                elemStiff_( m * dof + i, n * dof + i ) +=  Kentry * weight;
            }
        }
    }
    return;
}

//------------------------------------------------------------------------------
//! compute kernel of b(u,p)
template< typename NODE, typename SFUNU, typename SFUNP, typename MAT>
void solid::fem::ElementMixed<NODE,SFUNU,SFUNP,MAT>::
mixedFormIntegrand( const VecLDim & xi, const double & weight )
{
    // compute global derivatives of shape functions //! (d phi) / (dX)
    MatDimNN  Dphi;
    const double detG = this -> globalDerivatives_( xi, Dphi );

    // evaluate the pressure shape function
    VecNPN_ psi;
    this -> sfunPress_( xi, psi );

    // product phiP * div( DphiU )
    for ( unsigned m = 0; m < numNodes; m ++ ) {
        for ( unsigned n = 0; n < numPressNodes; n ++ ) {
            for ( unsigned i = 0; i < dof; i ++ ) {
                mixedMat_( m * dof + i, n ) += Dphi( i, m ) * psi( n ) * detG * weight;
            }
        }
    }
    return;
}

//------------------------------------------------------------------------------
template< typename NODE, typename SFUNU, typename SFUNP, typename MAT>
void solid::fem::ElementMixed<NODE,SFUNU,SFUNP,MAT>::
meanPressIntegrand( const VecLDim & xi, const double & weight )
{
    // get jacobian from element
    const double detG = this -> jacobian_( xi );

    // evaluate the pressure shape function
    VecNPN_ psi;
    this -> sfunPress_( xi, psi );

    // assemble outer product 
    stabilMat_ += detG * outer_prod( psi, psi ) * weight;

    return;
}


//------------------------------------------------------------------------------
//! interpolate pressures to non-pressure nodes
template< typename NODE, typename SFUNU, typename SFUNP, typename MAT>
void solid::fem::ElementMixed<NODE,SFUNU,SFUNP,MAT>::interpolatePressure()
{
    typedef boost::numeric::ublas::matrix_column<MatDimNN> column;
        
    SFUNU sfunu_;

    // get local coordinates of the interpolation nodes of U
    MatLDimNN ipointsU; 
    sfunu_.giveInterpolationPoints( ipointsU );

    // get pressure field
    VecNPN_ pressureValues;
    for ( unsigned pv = 0; pv < numPressNodes; pv ++ ) {
        pressureValues( pv ) = ElementBasic_::nodes_[ pv ] -> givePressure();
    }

    // go through all nodes and pass interpolated pressure to them
    for ( unsigned uv = 0; uv < numNodes; uv ++ ) {
        const bool wp = ElementBasic_::nodes_[ uv ] -> givePressureFlag( );
        if ( !wp ) {
            // get local coordinate of the node
            VecDim xi( column( ipointsU, uv ) );
            // compute interpolation function at that node
            VecNPN_ psi;
            sfunp_.evaluate( xi, psi );
            // compute interpolated pressure
            const double pressure = inner_prod(  psi, pressureValues );
            typename NODE::Vec1 aux; aux[0] = pressure;
            // pass value to node
            ElementBasic_::nodes_[ uv ] -> storePressure( aux );
            // set flag to true (??)
            ElementBasic_::nodes_[ uv ] -> setPressureFlag( );
        }
    }
}

