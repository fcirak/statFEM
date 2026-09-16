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

//! @file   PlaneStressDecorator.ipp
//! @author Thomas Rueberg
//! @date   2010

//------------------------------------------------------------------------------
//! iterate until state of plane stress is reached
void solid::material::PlaneStressDecorator::
iterateUntilPlaneStress( Mat3x3_ & F ) const
{
    // tolerance
    const double tol = 1.e-8; //? more sophisticed criterion needed

    while ( true ) {
        // obtain P
        Mat3x3_ P;  
        material_->FPKstress( F, P );

        // check criterion
        if ( std::fabs( P(2,2) ) < tol ) break;

        // obtain elasticity tensor according to passed deformation gradient
        Mat3x3x3x3_ C;
        material_->elasticityTensor( F, C );

        // update F_33
        F(2,2) = F(2,2) - P(2,2)/C(2,2)(2,2); // Newton iterate
    }
}


//------------------------------------------------------------------------------
// functions which simply mask the material behaviour

//! strain energy function
double solid::material::PlaneStressDecorator::strainEnergy( const Mat3x3_ & F ) const
{
    // iterate to state of plane stress
    Mat3x3_ planeF = F;
    this -> iterateUntilPlaneStress( planeF );
    
    // pass back the energy
    return material_->strainEnergy( planeF );
}

//! First Piola-Kirchhoff stress tensor
void solid::material::PlaneStressDecorator::FPKstress( const Mat3x3_ & F, Mat3x3_ & P ) const
{
    // iterate to state of plane stress
    Mat3x3_ planeF = F;
    this -> iterateUntilPlaneStress( planeF );

    // pass back the energy
    return material_->FPKstress( planeF, P );
}

//! elasticity tensor
void solid::material::PlaneStressDecorator::elasticityTensor( const Mat3x3_ & F, Mat3x3x3x3_ & C ) const
{
    // iterate to state of plane stress
    Mat3x3_ planeF = F;
    this -> iterateUntilPlaneStress( planeF );
    
    // pass back the energy
    return material_->elasticityTensor( planeF, C );
}

