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

//! @file   SVenant.ipp
//! @author Thomas Rueberg
//! @date   2010

#include <corlib/linalg.hpp>
#include <Eigen/Core>


//------------------------------------------------------------------------------
solid::material::SVenant::SVenant( )
    : MaterialBase( numberInternalVariables_, true)
{
    registerVariable( "youngModulus", youngModulus_ );
    registerVariable( "poissonRatio", poissonRatio_ );
    registerVariable( "massDensity",  massDensity_ );
}

solid::material::SVenant::~SVenant( ) { }

//------------------------------------------------------------------------------
void solid::material::SVenant::computeLame( double & lambda, double & mu ) const
{
    lambda = youngModulus_ * poissonRatio_ / 
        (1.0 + poissonRatio_) / (1.0 - 2.0 * poissonRatio_); 
    mu = youngModulus_ / 2.0 / (1.0 + poissonRatio_);
}


//------------------------------------------------------------------------------
/** Compute value of the Strain Energy function
 *  \f[
 *       W = \frac{\lambda}{2} (tr E)^2 + \mu tr (E E)
 *  \f]
 *  with \f$ E = 1/2 ( F^T F - I ) \f$.
 *  \param[in]  F Deformation gradient
 *  \retval     W strain energy
 */
double solid::material::SVenant::strainEnergy( const Mat3x3_ & F ) const
{
    // Lame parameters
    double lambda, mu;
    this -> computeLame( lambda, mu );
    // compute Euler-Lagrange strain tensor
    Eigen::MatrixXd eye_ = Eigen::MatrixXd::Identity( 3, 3 );
    
    Mat3x3_ E;
    E = 0.5 * ( ( F.transpose( ) * F ) - eye_ );
    
    // compute traces
    const double trE  = corlib::trace<3>( E );
    const Mat3x3_ E2 = E * E;
    const double trE2 = corlib::trace< 3 >( E2 );
    // return energy value
    return lambda / 2. * trE * trE + mu * trE2;
}

//------------------------------------------------------------------------------
/** Compute First Piola-Kirchhoff stress
 *  \f[
 *      P = F ( \lambda tr(E) I + 2 \mu E )
 *  \f]
 *  with the Euler-Lagrange strain tensor E.
 *  \param[in]  F Deformation gradient
 *  \param[out] P First Piola-Kirchhoff stress tensor
 */
void solid::material::SVenant::FPKstress( const Mat3x3_ & F, Mat3x3_ & P ) const
{
    // Lame parameters
    double lambda, mu;
    this -> computeLame( lambda, mu );
    // compute Euler-Lagrange strain tensor
    const Eigen::MatrixXd eye_ = Eigen::MatrixXd::Identity( 3, 3 );
    
    const Mat3x3_ E = 0.5 * ( ( F.transpose( ) * F ) - eye_ );
    
    // compute trace
    const double trE = corlib::trace< 3 >( E );
    // compute SPK stress
    P = lambda * trE * eye_ + 2. * mu * E;
    // convert to FPK stress
    P = F * P;
     return;
}

//------------------------------------------------------------------------------
/** Compute the elasticity tensor
 *  \f[
 *       C_{iJkL} = \lambda F_{iJ} F_{kL} + \mu F_{iL} F_{kJ}
 *                  + \mu (F F^T)_{ik} \delta_{JL} + S_{JL} \delta_{ik}
 *  \f]
 *  with \f$ S_{JL} \f$ such that \f$ P_{iL} = F_{iJ} S_{JL} \f$ and P as in
 *  SVenant::FPKstress. 
 *  \param[in]  F  Deformation gradient
 *  \param[out] C  Elasticity tensor
 */
void solid::material::SVenant::elasticityTensor( const Mat3x3_ & F, Mat3x3x3x3_ & C ) const
{
    // Lame parameters
    double lambda, mu;
    this -> computeLame( lambda, mu );
    // compupte Euler-Lagrange strain tensor
    const Eigen::MatrixXd eye_ = Eigen::MatrixXd::Identity( 3, 3 );

    const Mat3x3_ E = 0.5 * ( ( F.transpose( ) * F ) - eye_ );
    
    // compute SPK stress
    const double trE = corlib::trace< 3 >( E );
    const Mat3x3_ S = lambda * trE * eye_ + 2. * mu * E;
    // F times F^T
    const Mat3x3_ FFt = F * F.transpose( );
    
    // set up elasticity tensor
    double C_iJkL;
    for ( unsigned i = 0; i < 3; i ++ ) {
        for ( unsigned J = 0; J < 3; J ++ ) {
            for ( unsigned k = 0; k < 3; k ++ ) {
                for ( unsigned L = 0; L < 3; L ++ ) {
                    C_iJkL = lambda * F( i,J ) * F( k, L );
                    C_iJkL += mu * F(i,L) * F( k, J);
                    if ( J == L ) C_iJkL += mu * FFt( i, k );
                    if ( i == k ) C_iJkL += S( J, L);
                    
                    C( i, J)( k, L ) = C_iJkL;
                }
            }
        }
    }
    return;
}

