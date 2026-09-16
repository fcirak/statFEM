//
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//
//                 Computational Structural Mechanics Lab
//                         University of Cambridge
//
// This software is copyrighted and all rights are retained by the CSMLab. All
// use, disclosure, and/or reproduction of any part not expressly authorized by
// F Cirak is prohibited. (C) 2010.
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//

//! @author Burkhard Bornemann, Kosala Bandara
//! @date   2011

//------------------------------------------------------------------------------
template< typename BELEMENT, typename MAT >
void gshell::fem::ElementKLDynamic<BELEMENT,MAT>::massIntegrand(
    const VecLDim & xi, 
    const double & weight,
    Eigen::MatrixXd & elemMass
    ) const
{
    // number of shape functions
    const unsigned numFunctions = this->numFunctions( );
    const unsigned matSize = dof * numFunctions;
    FTL_VERIFY( elemMass.rows() == matSize and
                elemMass.cols() == matSize );

    // shape functions
    VecNF phi;          this->BasisElement::sfun( xi, phi );
    // 1st derivatives
    MatLDimNF dPhiDXi;  this->BasisElement::sfunGrad( xi, dPhiDXi );

    // Collect nodal quantities
    MatDimNF xNF;  this->BasisElement::supportNodeCoordinates( xNF );
    // Compute base
    double aRef[dim][dim];
    
    MatDimNF xNFT = xNF.transpose();
    MatLDimNF dPhiDXiT = dPhiDXi.transpose();
    const double detA = SurfaceGeometer::covariantBase( numFunctions, &(dPhiDXiT(0,0)),
                                                        &(xNFT(0,0)), aRef );

    // material mass density
    const double rho = material_->density( );  // kg/m^3

    // compute element mass matrix
    //
    // HINT: Only the translational inertia is computed.
    for ( unsigned l=0; l<numFunctions; ++l ) {
        for ( unsigned k=0; k<numFunctions; ++k ) {
            const double mass_dlck = thickness_ * rho * phi[ l ] * phi[ k ];
            for ( unsigned d=0; d<dof; ++d ) {
                const unsigned c = d;
                elemMass( l*dof+d, k*dof+c ) += mass_dlck * detA * weight;  // positive LHS
            }
        }
    }

    return;
}

//------------------------------------------------------------------------------
template< typename BELEMENT, typename MAT >
void gshell::fem::ElementKLDynamic<BELEMENT,MAT>::lumpedMassIntegrand(
    const VecLDim & xi, 
    const double & weight,
    Eigen::MatrixXd & lumpedMass
    ) const
{
    // number of shape functions
    const unsigned numFunctions = this->BasisElement::numFunctions( );
    const unsigned matSize = dof * numFunctions;
    FTL_VERIFY( lumpedMass.rows() == matSize and
                lumpedMass.cols() == 1 );

    // shape functions
    VecNF phi;          this->BasisElement::sfun( xi, phi );
    // 1st derivatives
    MatLDimNF dPhiDXi;  this->BasisElement::sfunGrad( xi, dPhiDXi );

    // Collect nodal quantities
    MatDimNF xNF;  this->BasisElement::supportNodeCoordinates( xNF );
    // Compute base
    double aRef[dim][dim];
    MatDimNF xNFT = xNF.transpose();
    MatLDimNF dPhiDXiT = dPhiDXi.transpose();
    const double detA = SurfaceGeometer::covariantBase( numFunctions, &(dPhiDXiT(0,0)),
                                                        &(xNFT(0,0)), aRef );

    // material mass density per volume(!)
    const double rho = material_->density( );  // kg/m^3

    // compute element mass matrix
    //
    // HINT: Only the translational inertia is computed.
    for ( unsigned l=0; l<numFunctions; ++l ) {
        const double mass_dl = thickness_ * rho * phi[ l ];// * phi[ k ];
        for ( unsigned d=0; d<dof; ++d ) {
            lumpedMass( l*dof+d, 0 ) += mass_dl * detA * weight;  // positive LHS
        }
    }

    return;
}

//------------------------------------------------------------------------------
template< typename BELEMENT, typename MAT >
void gshell::fem::ElementKLDynamic<BELEMENT,MAT>::nodalVelocities(
    MatDimNF & vel
    ) const
{
    this->nodalQuantity_( vel, std::mem_fun( &Node::giveVelocities ) );
    return;
}

//------------------------------------------------------------------------------
template< typename BELEMENT, typename MAT >
void gshell::fem::ElementKLDynamic<BELEMENT,MAT>::nodalAccelerations(
    MatDimNF & acc
    ) const
{
    this->nodalQuantity_( acc, std::mem_fun( &Node::giveAccelerations ) );
    return;
}

