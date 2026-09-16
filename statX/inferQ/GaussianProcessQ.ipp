// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//
//                 Computational Structural Mechanics Lab
//                         University of Cambridge
//
// This software is copyrighted and all rights are retained by the University of
// Cambridge. All use, disclosure, and/or reproduction of any part not expressly
// authorised by the University of Cambridge is prohibited. (C) 2010 --
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

//! @file GaussianProcessQ.ipp

//------------------------------------------------------------------------------
void statX::inferQ::GaussianProcessQ::
setHypIdPre( const std::string & preamble )
{
    hypIdPre_ = preamble;

    return;
}

//------------------------------------------------------------------------------
/**
 * - Step 1a: solve \f$ ( \sigma_e^2 Q_v + n_o P_v^\trans P_v ) X =
 *            P_v^\trans \left( \sum\limits_i^{n_o} y_i \right) \f$
 *            and compute \f$ P_v X \f$.
 * - Step 1b: compute \f$ - \frac{1}{2} \sum\limits_j^{n_o} \left( y_j ^ \trans
 *            \left( \sigma_e^{-2} ( y_j - P_v X ) \right) \right) \f$.
 * - Step 2a: compute \f$ \log \vert \sigma_e^2 Q_v + n_o P_v^\trans P_v \vert \f$.
 * - Step 2b: compute \f$ \log \vert  Q_v  \vert \f$.
 * - Step 2c: compute \f$ - \frac{1}{2} \log \left( ( \sigma_e^2 )^{n_y n_o - n_u}
 *            \frac{ \vert \sigma_e^2 Q_v + n_o P_v^\trans P_v \vert }
 *            { \vert Q_v \vert } \right) \f$.
 *
 * - Step 3:  compute \f$ - \frac{n_y n_o}{2} \log 2 \pi \f$.
 */
double statX::inferQ::GaussianProcessQ::
giveLogLikelihood( const HpsValueMap & hParams,
                   const Eigen::SparseMatrix<double> & qMat,
                   const Eigen::SparseMatrix<double> & prMat )
{
    // size check
    FTL_VERIFY_DESCRIPTIVE( pMat_.cols( ) == qMat.rows( ),
                            "Number of mesh points not matched.\n" );

    // get observation error
    const double eps = this->extractEps_( hParams );

    // number of data points and number of observations
    const unsigned numSensors  = obsValues_.rows( ); // n_y
    const unsigned numReadings = obsValues_.cols( ); // n_o

    // get effective projection matrix
    const Eigen::SparseMatrix<double> pMatV = pMat_ * prMat;

    // step 1a
    solver_.compute( eps * eps * qMat +
                     numReadings * pMatV.transpose( ) * pMatV );
    FTL_VERIFY( solver_.info( ) == Eigen::Success );
    Eigen::VectorXd sumObsValues = Eigen::VectorXd::Zero( numSensors );
    for ( unsigned i = 0; i < numReadings; ++i )
        sumObsValues += obsValues_.col( i );
    const Eigen::VectorXd fixedVec =
                  pMatV * solver_.solve( pMatV.transpose( ) * sumObsValues );
    FTL_VERIFY( solver_.info( ) == Eigen::Success );

    // step 1b
    double tmp1 = 0.;
    for ( unsigned i = 0; i < numReadings; ++i ) {
        tmp1 += ( obsValues_.col( i ) ).dot(
                ( 1. / ( eps * eps ) ) * ( obsValues_.col( i ) - fixedVec ) );
    }
    tmp1 *= -0.5;

    // step 2a
    Eigen::MatrixXd lMat = solver_.matrixL( );
    Eigen::ArrayXd lMatDiag = lMat.diagonal( );
    const double tmp2a = 2. * lMatDiag.log( ).sum( );

    // step 2b
    solver_.compute( qMat );
    FTL_VERIFY( solver_.info( ) == Eigen::Success );
    lMat = solver_.matrixL( );
    lMatDiag = lMat.diagonal( );
    const double tmp2b = 2. * lMatDiag.log( ).sum( );

    // step 2c
    const int diff = numSensors * numReadings - qMat.rows( ); // can be +ve or -ve
    const double tmp2 = -0.5 * ( diff * std::log( eps * eps ) + tmp2a - tmp2b );

    // step 3
    const double tmp3 = -0.5 * numReadings * numSensors * std::log( 2. * M_PI );

    return tmp1 + tmp2 + tmp3;
}

//------------------------------------------------------------------------------
/**
 * The mean vector of the posterior \f$ p( u \vert Y ) \f$ is
 * \f{align}{
 *        u_{\vert Y} &= P_r v_{\vert Y} \quad \text{where} \\
 *        v_{\vert Y} &= ( \sigma_e^2 Q_v + n_o P_v^\trans P_v )^{-1}
 *                       P_v^\trans \left( \sum\limits_i^{n_o} y_i \right)
 * \f}
 */
void statX::inferQ::GaussianProcessQ::
predictMean( const HpsValueMap & hParams,
             const Eigen::SparseMatrix<double> & qMat,
             const Eigen::SparseMatrix<double> & prMat,
             Eigen::VectorXd & uVec )
{
    // size check
    FTL_VERIFY_DESCRIPTIVE( pMat_.cols( ) == qMat.rows( ),
                            "Number of mesh points not matched.\n" );

    // get observation error
    const double eps = this->extractEps_( hParams );

    // get effective projection matrix
    const Eigen::SparseMatrix<double> pMatV = pMat_ * prMat;

    // number of data points and number of observations
    const unsigned numSensors  = obsValues_.rows( ); // n_y
    const unsigned numReadings = obsValues_.cols( ); // n_o

    // posterior mean of auxiliary variable
    solver_.compute( eps * eps * qMat +
                     numReadings * pMatV.transpose( ) * pMatV );
    FTL_VERIFY( solver_.info( ) == Eigen::Success );
    Eigen::VectorXd sumObsValues = Eigen::VectorXd::Zero( numSensors );
    for ( unsigned i = 0; i < numReadings; ++i )
        sumObsValues += obsValues_.col( i );
    const Eigen::VectorXd vVec =
            solver_.solve( pMatV.transpose( ) * sumObsValues );
    FTL_VERIFY( solver_.info( ) == Eigen::Success );

    // project to the mean vector of inferred variable
    uVec.resize( pMat_.cols( ) );
    uVec = prMat * vVec;

    return;
}

//------------------------------------------------------------------------------
/**
 * The covariance matrix of the posterior \f$ p( u \vert Y ) \f$ is
 * \f{align}{
 *        C_{u \vert Y} &= P_r C_{v \vert Y} P_r^\trans \quad \text{where} \\
 *        C_{v \vert Y} &= \sigma_e^2 ( \sigma_e^2 Q_v + n_o P_v^\trans P_v )^{-1}
 * \f}
 */
void statX::inferQ::GaussianProcessQ::
predictCovariance( const HpsValueMap & hParams,
                   const Eigen::SparseMatrix<double> & qMat,
                   const Eigen::SparseMatrix<double> & prMat,
                   Eigen::MatrixXd & uCovMat )
{
    // size check
    FTL_VERIFY_DESCRIPTIVE( pMat_.cols( ) == qMat.rows( ),
                            "Number of mesh points not matched.\n" );

    // get observation error
    const double eps = this->extractEps_( hParams );

    // get effective projection matrix
    const Eigen::SparseMatrix<double> pMatV = pMat_ * prMat;

    // number of observations
    const unsigned numReadings = obsValues_.cols( ); // n_o

    // compute the inverse
    solver_.compute( eps * eps * qMat +
                     numReadings * pMatV.transpose( ) * pMatV );
    FTL_VERIFY( solver_.info( ) == Eigen::Success );
    Eigen::SparseMatrix<double> idenMat( pMat_.cols( ), pMat_.cols( ) );
    idenMat.setIdentity( );
    const Eigen::MatrixXd invMat = solver_.solve( idenMat );
    FTL_VERIFY( solver_.info( ) == Eigen::Success );

    // scale by sigma_e ^ 2 and project to covariance matrix of actual variable
    uCovMat.resize( pMat_.cols( ), pMat_.cols( ) );
    uCovMat = prMat * ( eps * eps * invMat ) * prMat.transpose( );

    return;
}

//------------------------------------------------------------------------------
void statX::inferQ::GaussianProcessQ::
predictVariance( const HpsValueMap & hParams,
                 const Eigen::SparseMatrix<double> & qMat,
                 const Eigen::SparseMatrix<double> & prMat,
                 Eigen::VectorXd & uVarVec )
{
    Eigen::MatrixXd uCovMat;
    this->predictCovariance( hParams, qMat, prMat, uCovMat );
    uVarVec = uCovMat.diagonal( );
    return;
}

//------------------------------------------------------------------------------
double statX::inferQ::GaussianProcessQ::
extractEps_( const HpsValueMap & hParams ) const
{
    const auto it = hParams.find( hypIdPre_ + hypIds_[ 0 ] );
    FTL_VERIFY( it != hParams.end( ) );
    return it->second;
}
