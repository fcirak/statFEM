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

//! @file   Neohookean.ipp
//! @author Thomas Rueberg
//! @date   2010

#include <corlib/linalg.hpp>

//------------------------------------------------------------------------------
solid::material::Neohookean :: Neohookean( )
    : MaterialBase( numberInternalVariables_, true)
{
    registerVariable( "youngModulus", youngModulus_ );
    registerVariable( "poissonRatio", poissonRatio_ );
    registerVariable( "massDensity",  massDensity_ );
}

solid::material::Neohookean :: ~Neohookean( ) { }

//------------------------------------------------------------------------------
void solid::material::Neohookean :: computeLame( double & lambda, double & mu ) const
{
    lambda = youngModulus_ * poissonRatio_ / 
        (1.0 + poissonRatio_) / (1.0 - 2.0 * poissonRatio_); 
    mu = youngModulus_ / 2.0 / (1.0 + poissonRatio_);
    return;
}


//------------------------------------------------------------------------------
/** Compute value of the Strain Energy function
 *  \f[
 *         W = \frac{\lambda}{2} (\log J)^2 -\mu \log J +\frac{\mu}{2} (tr(C) - 3)
 *  \f]
 *  with \f$ C = F^T F \f$ and \f$ J = det F \f$.
 *  \param[in]  F Deformation Gradient
 *  \retval     W the elastic strain energy
 */
double solid::material::Neohookean::strainEnergy( const Mat3x3_ & F ) const
{
    // Lame parameters
    double lambda, mu;
    this -> computeLame( lambda, mu );
    // compute Right Cauchy-Green tensor
    Mat3x3_ C = F.transpose( ) * F;
    // compute log( det( F ) )
    const double logJ = log( corlib::determinant( F ) );
    // trace of C
    const double trC = corlib::trace<3>( C );
    
    // return energy value
    return lambda / 2. * logJ * logJ - mu * logJ + mu/2. * (trC -3.);
}

//------------------------------------------------------------------------------
/** Compute the first Piola-Kirchhoff stress tensor
 *  \f[
 *       P_{iJ} = F_{iK} [ (\lambda \log J - \mu) C^{-1}_{KJ} + \mu \delta_{KJ}]
 *  \f]
 *  with \f$ C_{IJ} = F_{kI} F_{kJ} \f$ and \f$ J = det F \f$.
 *  \param[in]  F Deformation gradient
 *  \param[out] P First Piola-Kirchhoff stress tensor
 */
void solid::material::Neohookean :: FPKstress( const Mat3x3_ & F, Mat3x3_ & P ) const
{
    // Lame parameters
    double lambda, mu;
    this -> computeLame( lambda, mu );
    // compute inverse Right Cauchy Green 
    Mat3x3_ Ci = F.transpose( ) * F;
    corlib::inverse( Ci );
    // compute log( det( F ) )
    const double logJ = log( corlib::determinant( F ) );
    // compute 2nd PK stress
    const double factor = lambda * logJ - mu;
    Eigen::MatrixXd eye_ = Eigen::MatrixXd::Identity( 3, 3 );
    P = factor * Ci + mu * eye_;
    // 1st PK stress:
    P = F * P;
    return;
}

//------------------------------------------------------------------------------
/** Compute the elasticity tensor
 *  \f[
 *       C_{iJkL} = S_{JL} \delta_{ik} + F_{iM} \hat{C}_{MJNL} F_{kN}
 *       \qquad 
 *       \hat{C}_{MJNL} = \lambda C^{-1}_{MJ} C^{-1}_{NL} + (\mu - \lambda \log J)
 *                        (C^{-1}_{MN} C^{-1}_{JL} + C^{-1}_{ML} C^{-1}_{JN})
 *  \f]
 *  with \f$ S_{JL} \f$ such that \f$ P_{iL} = F_{iJ} S_{JL} \f$ and P as in
 *  Neohookean::FPKstress. Don't confuse \f$ C_{iJkL} \f$ (fourth order elasticity 
 *  tensor) with \f$ C_{IJ} \f$ (second order Right Cauchy Green stretch tensor)!!
 *  \param[in]  F  Deformation gradient
 *  \param[out] C  Elasticity tensor
 */
void solid::material::Neohookean::elasticityTensor( const Mat3x3_ & F, Mat3x3x3x3_ & C ) const
{
    // Lame parameters
    double lambda, mu;
    this -> computeLame( lambda, mu );
    // compute inverse Right Cauchy Green 
    Mat3x3_ Ci = F.transpose( ) * F;
    corlib::inverse( Ci );
    // compute log( det( F ) )
    const double logJ = log( corlib::determinant( F ) );
    // compute 2nd PK stress
    const double factor = lambda * logJ - mu;
    Eigen::MatrixXd eye_ = Eigen::MatrixXd::Identity( 3, 3 );
    Mat3x3_ SPK = factor * Ci + mu * eye_;
    // compute (d^2 W)/(d C^2)
    Mat3x3x3x3_ Ctmp;
    for (int L = 0; L < 3; ++L ) {
        for (int  N = 0; N < 3; ++N ) {
            for (int J = 0; J < 3; ++J ) {
                for (int M = 0; M < 3; ++M ) {
                    Ctmp(M, J)(N, L) = lambda * Ci(M, J) * Ci(N, L)
                        - factor * ( Ci(M, N) * Ci(J, L) + Ci(M, L) * Ci(J, N) );
                }
            }
        }	
    }
    
    // change to (d^2 W)/(d F^2)
    for (int L = 0; L < 3; ++L ) {
        for (int k = 0; k < 3; ++k ) {
            for (int J = 0; J < 3; ++J ) {
                for (int i = 0; i < 3; ++i ) {
                    double CiJkL = 0.0;
                    for (int N = 0; N < 3; ++N )  {
                        for (int M = 0; M < 3; ++M ) {
                            CiJkL += F(i,M) * Ctmp(M,J)(N,L) * F(k,N);
                        }
                    }
                    if ( i == k ) {
                        CiJkL += SPK( J, L );
                    }
                    C(i,J)(k,L) = CiJkL;
                }
            }
        }
    }
    return;
}

