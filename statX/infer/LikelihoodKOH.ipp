// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//
//                 Computational Structural Mechanics Lab
//                         University of Cambridge
//
// This software is copyrighted and all rights are retained by the University of
// Cambridge. All use, disclosure, and/or reproduction of any part not expressly
// authorised by the University of Cambridge is prohibited. (C) 2010 --
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

//! @file LikelihoodKOH.ipp
#include <vector>

#include <corlib/verify.hpp>
#include <statX/infer/NormalSampler.hpp>

//------------------------------------------------------------------------------
template <typename FORWARDGAUSS, typename MISMATCHCOV>
void statX::infer::LikelihoodKOH<FORWARDGAUSS, MISMATCHCOV>::computeTotalCovariance_ (
        const HpsValueMap & hyperParam, Eigen::MatrixXd & covTotal )
{    
    //! compute the covariance of the FE solution
    Eigen::MatrixXd covU;
    forwardGauss_ -> giveProjectedCovarU( covU );

    //! get rho
    const auto it = hyperParam.find( hypIdPre_ + hypId_ );
    FTL_VERIFY( it != hyperParam.end( ) );
    const double rho = it->second;

    mismatchCov_ -> giveCovariance( hyperParam, covTotal );
    
    FTL_VERIFY( covU.rows( ) == covTotal.rows( ) );
    FTL_VERIFY( covU.cols( ) == covTotal.cols( ) );

    covTotal += rho * rho * covU;

    return;
}

//------------------------------------------------------------------------------
template <typename FORWARDGAUSS, typename MISMATCHCOV>
void statX::infer::LikelihoodKOH<FORWARDGAUSS, MISMATCHCOV>::giveTypicalReadings (
        const HpsValueMap & hyperParam, unsigned numReadings,
        Eigen::MatrixXd & yTypicalMat, unsigned seed )
{
    //! get solution mean
    Eigen::VectorXd feMean;
    forwardGauss_ -> giveProjectedMeanU( feMean );

    //! compute total covariance
    Eigen::MatrixXd covTotal;
    computeTotalCovariance_( hyperParam, covTotal );

    //! instantiate a normal sampler
    statX::infer::NormalSampler normalSampler( seed );

    //! get rho
    const auto it = hyperParam.find( hypIdPre_ + hypId_ );
    FTL_VERIFY( it != hyperParam.end( ) );
    const double rho = it->second;

    yTypicalMat.resize( feMean.rows( ), numReadings );
    for ( unsigned i = 0; i < numReadings; ++ i ) {
        Eigen::VectorXd yTypicalVec;
        normalSampler.giveSample( rho * feMean, covTotal, yTypicalVec );
        yTypicalMat.col( i ) = yTypicalVec;
    }

    return;
}

//------------------------------------------------------------------------------
template <typename FORWARDGAUSS, typename MISMATCHCOV>
double statX::infer::LikelihoodKOH<FORWARDGAUSS, MISMATCHCOV>::giveLogLikelihood (
        const HpsValueMap & hyperParam )
{
    //! compute total covariance
    Eigen::MatrixXd covTotal;
    computeTotalCovariance_( hyperParam, covTotal );

    //! get solution mean
    Eigen::VectorXd feMean;
    forwardGauss_ -> giveProjectedMeanU( feMean );

    const unsigned numSensors = sensorValues_.rows( );
    FTL_VERIFY( numSensors == feMean.rows( ) );

    const double tmp1 = -0.5 * numSensors * log( 2. * M_PI );

    // second term (determinant) - depends on the number of sensors
    Eigen::LLT < Eigen::MatrixXd > lltOfMk( covTotal );
    FTL_VERIFY( lltOfMk.info( ) == Eigen::Success );

    const Eigen::MatrixXd cholMca = lltOfMk.matrixL( );
    Eigen::ArrayXd cholMcaDiag = cholMca.diagonal( );

    const double tmp2 = -0.5 * ( 2. * cholMcaDiag.log( ).sum( ) );

    // third term - depends on the number of sensors and readings
    const unsigned numReadings = sensorValues_.cols( );

    //! get rho
    const auto it = hyperParam.find( hypIdPre_ + hypId_ );
    FTL_VERIFY( it != hyperParam.end( ) );
    const double rho = it->second;

    double tmp3 = 0.0;

    for ( auto i = 0; i < numReadings; ++ i ) {
        const Eigen::VectorXd diff = sensorValues_.col( i ) - rho * feMean;
        const double tmp3Comp = -0.5 * diff.dot( lltOfMk.solve( diff ) );
        FTL_VERIFY( tmp3Comp < 0. );
        tmp3 += tmp3Comp;
    }

    double value = numReadings * ( tmp1 + tmp2 ) + tmp3;

    return value;
}

//------------------------------------------------------------------------------
template <typename FORWARDGAUSS, typename MISMATCHCOV>
double statX::infer::LikelihoodKOH<FORWARDGAUSS, MISMATCHCOV>::giveLogLikelihoodTimeSeries (
        const HpsValueMap &hyperParam )
{
    // number of readings
    const unsigned numReadings = sensorValues_.cols( );
    
    //! first term
    const unsigned numSensors = sensorValues_.rows( );   
    const double tmp1 = -0.5 * numSensors * log( 2. * M_PI );
    
    //! get \f$ \sigma_e \f$ and compute \f$ C_e \f$
    auto itE = hyperParam.find( "mismatch.sensorEps" );
    FTL_VERIFY( itE != hyperParam.end( ) );
    const double eps = itE -> second;
    
    Eigen::MatrixXd covE;
    covE.resize( numSensors, numSensors );
    covE = eps * eps * Eigen::MatrixXd::Identity( numSensors, numSensors );

    //! get unscaled \f$ C_d \f$
    Eigen::MatrixXd covD;
    mismatchCov_->giveCovariance( hyperParam, covD );
    
    //! covD contains \f$ C_e \f$, remove it
    covD -= covE;
    
    //! get projected \f$ C_u \f$
    Eigen::MatrixXd covU;
    forwardGauss_->giveProjectedCovarU( covU );
    
    //! get \f$ \rho \f$
    const auto it = hyperParam.find( hypIdPre_ + hypId_ );
    FTL_VERIFY( it != hyperParam.end( ) );
    const double rho = it->second;

    //! second and third terms
    double tmp2 = 0.0;
    double tmp3 = 0.0;
    FTL_VERIFY( mismatchScaling_.rows() == numReadings );
    for ( auto i = 0; i < numReadings; ++ i ) {
        
        //! get solution mean
        Eigen::VectorXd feMean;
        forwardGauss_->giveProjectedMeanU( feMean, i );
        FTL_VERIFY( numSensors == feMean.rows( ) );
        
        //! compute total covariance
        Eigen::MatrixXd covTotal;
        covTotal = ( rho * rho * covU ) + ( mismatchScaling_( i ) * mismatchScaling_( i ) * covD ) + covE;
        
        // second term (determinant) - depends on the number of sensors
        Eigen::LLT < Eigen::MatrixXd > lltOfMk( covTotal );
        FTL_VERIFY( lltOfMk.info( ) == Eigen::Success );

        const Eigen::MatrixXd cholMca = lltOfMk.matrixL( );
        Eigen::ArrayXd cholMcaDiag = cholMca.diagonal( );

        tmp2 += -0.5 * ( 2. * cholMcaDiag.log( ).sum( ) );

        const Eigen::VectorXd diff = sensorValues_.col( i ) - rho * feMean;
        const double tmp3Comp = -0.5 * diff.dot( lltOfMk.solve( diff ) );
        FTL_VERIFY( tmp3Comp < 0. );
        tmp3 += tmp3Comp;

    }
    
    double value = ( numReadings * tmp1 ) + tmp2 + tmp3;

    return value;
}

//------------------------------------------------------------------------------
template <typename FORWARDGAUSS, typename MISMATCHCOV>
void statX::infer::LikelihoodKOH<FORWARDGAUSS, MISMATCHCOV>::giveLogLikelihoodGrad (
        const HpsValueMap & hyperParam, HpsValueMap & mllDeriv )
{
    //! compute total covariance
    Eigen::MatrixXd covTotal;
    computeTotalCovariance_( hyperParam, covTotal );

    //! compute gradient of log marginal likelihood with respect to rho
    const auto it = hyperParam.find( hypIdPre_ + hypId_ );
    FTL_VERIFY( it != hyperParam.end( ) );
    const double rho = it->second;

    //! get the mean of the FE solution
    Eigen::VectorXd feMean;
    forwardGauss_ -> giveProjectedMeanU( feMean );
    
    const Eigen::VectorXd diff = sensorValues_ - rho * feMean;

    //! get the covariance of the FE solution
    Eigen::MatrixXd covU;
    forwardGauss_ -> giveProjectedCovarU( covU );
    Eigen::MatrixXd covTotalDrho = 2. * rho * covU;

    Eigen::MatrixXd covTotalInv = covTotal.inverse( );

    mllDeriv[ hypIdPre_ + hypId_ ] = -0.5
            * ( covTotalInv * covTotalDrho ).trace( )
            + feMean.dot( covTotalInv * diff )
            + 0.5 * diff.dot( covTotalInv * covTotalDrho * covTotalInv * diff );

    //! compute gradient of log likelihood with respect to sigma
    HpsMatMap hpsMatMap;
    Eigen::MatrixXd mismatchDmag, mismatchDl;
    mismatchCov_ -> giveCovarianceDerv( hyperParam, hpsMatMap );

    std::vector < std::string > hypIds;
    for ( auto it = hpsMatMap.begin( ); it != hpsMatMap.end( ); ++ it ) {
        hypIds.push_back( it->first );
    }

    mismatchDmag = hpsMatMap[ hypIdPre_ + hypIds[ 0 ] ];
    mismatchDl = hpsMatMap[ hypIdPre_ + hypIds[ 1 ] ];

    mllDeriv[ hypIdPre_ + hypIds[ 0 ] ] = -0.5
            * ( covTotalInv * mismatchDmag ).trace( )
            + 0.5 * diff.dot( covTotalInv * mismatchDmag * covTotalInv * diff );

    mllDeriv[ hypIdPre_ + hypIds[ 1 ] ] = -0.5
            * ( covTotalInv * mismatchDl ).trace( )
            + 0.5 * diff.dot( covTotalInv * mismatchDl * covTotalInv * diff );

    return;
}
    
//------------------------------------------------------------------------------
template <typename FORWARDGAUSS, typename MISMATCHCOV>
void statX::infer::LikelihoodKOH<FORWARDGAUSS, MISMATCHCOV>::giveFEconditionedData (
        const HpsValueMap & hyperParam, Eigen::VectorXd & feMeanCondData,
        Eigen::MatrixXd & feCovarCondData ) const
{
    
    //! compute \f$ K_s = C_d + C_e \f$
    Eigen::MatrixXd matKs;
    mismatchCov_ -> giveCovariance( hyperParam, matKs );

    //! compute \f$ K_s = ( C_d + C_e )^{-1} \f$
    Eigen::FullPivLU < Eigen::MatrixXd > luKs( matKs );
    FTL_VERIFY( luKs.isInvertible( ) );
    Eigen::MatrixXd matKsInv = luKs.inverse( );

    //! get \f$ C_u = I A^{-1} C_f  A^{-T} I^T \f$
    Eigen::MatrixXd covUNotProjected;
    Eigen::MatrixXd matIdentity = Eigen::MatrixXd::Identity( pMat_.cols(), pMat_.cols() );
    forwardGauss_ -> setProjectionMatrix( matIdentity );
    forwardGauss_ -> updateSolution( hyperParam );
    forwardGauss_ -> giveProjectedCovarU( covUNotProjected );

    //! compute \$ C_u^{-1} \$
    Eigen::MatrixXd matCuInv;
    Eigen::FullPivLU < Eigen::MatrixXd > luCu( covUNotProjected );
    FTL_VERIFY( luCu.isInvertible( ) );
    matCuInv = luCu.inverse( );

    //! get rho and numReadings
    const auto it = hyperParam.find( hypIdPre_ + hypId_ );
    FTL_VERIFY( it != hyperParam.end( ) );
    const double rho = it->second;
    const unsigned numReadings = sensorValues_.cols( );

    //! compute matrix \f$ C_{u|y} = ( \rho^2 n_o K_s^{-1} + ( P C_u P^T )^{-1} )^{-1} \f$
    const Eigen::MatrixXd covSum = ( ( rho * rho * numReadings ) * 
                                       pMat_.transpose() * matKsInv * pMat_ ) +
                                       matCuInv;
    Eigen::FullPivLU < Eigen::MatrixXd > luSum( covSum );
    FTL_VERIFY( luSum.isInvertible( ) );
    feCovarCondData = luSum.inverse( );

    //! compute \f$ \sum_k y_k \f$
    Eigen::VectorXd sumSensorValues = sensorValues_.col( 0 );
    for ( unsigned i = 1; i < numReadings; ++ i ) {
        sumSensorValues += sensorValues_.col( i );
    }

    //! compute the FE mean conditioned on data \f$ \overline u_{u|y} \f$
    Eigen::VectorXd feMean;
    forwardGauss_ -> giveProjectedMeanU( feMean );
    feMeanCondData = feCovarCondData
            * ( rho * pMat_.transpose( ) * matKsInv * sumSensorValues + matCuInv * feMean );
    
    return;
}

//------------------------------------------------------------------------------
template <typename FORWARDGAUSS, typename MISMATCHCOV>
void statX::infer::LikelihoodKOH<FORWARDGAUSS, MISMATCHCOV>::giveProjFEconditionedData (
        const HpsValueMap & hyperParam, Eigen::VectorXd & feMeanCondData,
        Eigen::MatrixXd & feCovarCondData ) const
{
    
    //! compute \f$ K_s = C_d + C_e \f$
    Eigen::MatrixXd matKs;
    mismatchCov_ -> giveCovariance( hyperParam, matKs );

    //! compute \f$ K_s = ( C_d + C_e )^{-1} \f$
    Eigen::FullPivLU < Eigen::MatrixXd > luKs( matKs );
    FTL_VERIFY( luKs.isInvertible( ) );
    Eigen::MatrixXd matKsInv = luKs.inverse( );

    //! get \f$ P C_u P^T \f$
    Eigen::MatrixXd covU;
    forwardGauss_ -> giveProjectedCovarU( covU );

    //! compute \$ ( P C_u P^T )^{-1} \$
    Eigen::MatrixXd matCuInv;
    Eigen::FullPivLU < Eigen::MatrixXd > luCu( covU );
    FTL_VERIFY( luCu.isInvertible( ) );
    matCuInv = luCu.inverse( );

    //! get rho and numReadings
    const auto it = hyperParam.find( hypIdPre_ + hypId_ );
    FTL_VERIFY( it != hyperParam.end( ) );
    const double rho = it->second;
    const unsigned numReadings = sensorValues_.cols( );

    //! compute matrix \f$ C_{u|y} = ( \rho^2 n_o K_s^{-1} + ( P C_u P^T )^{-1} )^{-1} \f$
    const Eigen::MatrixXd covSum = ( ( rho * rho * numReadings ) * matKsInv )
                                 + matCuInv;
    Eigen::FullPivLU < Eigen::MatrixXd > luSum( covSum );
    FTL_VERIFY( luSum.isInvertible( ) );
    feCovarCondData = luSum.inverse( );

    //! compute \f$ \sum_k y_k \f$
    Eigen::VectorXd sumSensorValues = sensorValues_.col( 0 );
    for ( unsigned i = 1; i < numReadings; ++ i ) {
        sumSensorValues += sensorValues_.col( i );
    }

    //! compute the FE mean conditioned on data \f$ \overline u_{u|y} \f$
    Eigen::VectorXd feMean;
    forwardGauss_ -> giveProjectedMeanU( feMean );
    feMeanCondData = feCovarCondData
            * ( rho * matKsInv * sumSensorValues + matCuInv * feMean );
    
    return;
}

//------------------------------------------------------------------------------
template <typename FORWARDGAUSS, typename MISMATCHCOV>
void statX::infer::LikelihoodKOH<FORWARDGAUSS, MISMATCHCOV>::giveProjFEconditionedDataTimeSeries (
        const HpsValueMap & hyperParam, Eigen::VectorXd & feMeanCondData,
        Eigen::MatrixXd & feCovarCondData, unsigned index ) const
{
    
    //! compute \f$ K_s = C_d + C_e \f$
    Eigen::MatrixXd matKs;
    mismatchCov_ -> giveCovariance( hyperParam, matKs );
    
    //! get eps and compute covE
    auto itE = hyperParam.find( "mismatch.sensorEps" );
    FTL_VERIFY( itE != hyperParam.end( ) );
    const double eps = itE -> second;
    Eigen::MatrixXd covE = eps * eps * Eigen::MatrixXd::Identity( matKs.rows( ), matKs.rows( ) );   
    
    //! consider mismatch scaling
    matKs -= covE;
    double mismatchScale = mismatchScaling_( index );
    matKs *= ( mismatchScale * mismatchScale );
    matKs += covE;

    //! compute \f$ K_s = ( C_d + C_e )^{-1} \f$
    Eigen::FullPivLU < Eigen::MatrixXd > luKs( matKs );
    FTL_VERIFY( luKs.isInvertible( ) );
    Eigen::MatrixXd matKsInv = luKs.inverse( );

    //! get \f$ P C_u P^T \f$
    Eigen::MatrixXd covU;
    forwardGauss_ -> giveProjectedCovarU( covU );

    //! compute \$ ( P C_u P^T )^{-1} \$
    Eigen::MatrixXd matCuInv;
    Eigen::FullPivLU < Eigen::MatrixXd > luCu( covU );
    FTL_VERIFY( luCu.isInvertible( ) );
    matCuInv = luCu.inverse( );

    //! get rho
    const auto it = hyperParam.find( hypIdPre_ + hypId_ );
    FTL_VERIFY( it != hyperParam.end( ) );
    const double rho = it->second;

    //! compute matrix \f$ C_{u|y} = ( \rho^2 K_s^{-1} + ( P C_u P^T )^{-1} )^{-1} \f$
    //! here n_o = 1, i.e. consider only data at one time instance
    const Eigen::MatrixXd covSum = ( ( rho * rho ) * matKsInv )
                                 + matCuInv;
    Eigen::FullPivLU < Eigen::MatrixXd > luSum( covSum );
    FTL_VERIFY( luSum.isInvertible( ) );
    feCovarCondData = luSum.inverse( );
    
    //! compute matrix \f$ S = ( (1 / \rho^2 n_o) K_s ) + P C_u P^T \f$
    Eigen::MatrixXd matSinv;
    {
        Eigen::MatrixXd matS = ( matKs / ( rho * rho ) ) + covU;
        Eigen::FullPivLU<Eigen::MatrixXd> luMatS( matS );
        FTL_VERIFY( luMatS.isInvertible( ) );
        matSinv = luMatS.inverse( );
    }
    
    //! Identity matrix
    Eigen::MatrixXd matI;
    matI = Eigen::MatrixXd::Identity( pMat_.rows( ), pMat_.rows( ) );
    
    //! get the projected mean \f$ Pu \f$
    Eigen::VectorXd projFeMean;
    forwardGauss_ -> giveProjectedMeanU( projFeMean, index );
    
    //! compute the FE mean conditioned on data \f$ \overline u_{u|y} \f$
    Eigen::VectorXd vecA, vecB;
    vecA = feCovarCondData * rho * matKsInv * sensorValues_.col( index );
    vecB = ( matI - covU * matSinv ) * projFeMean;
    feMeanCondData = vecA + vecB;
    
    return;
}

//------------------------------------------------------------------------------
template <typename FORWARDGAUSS, typename MISMATCHCOV>
void statX::infer::LikelihoodKOH<FORWARDGAUSS, MISMATCHCOV>::giveFEconditionedDataWoodbury (
        const HpsValueMap & hyperParam, Eigen::VectorXd & feMeanCondData,
        Eigen::MatrixXd & feCovarCondData ) const
{

    //! Woodbury expansion of p(u|y)
    //! compute \f$ K_s = C_d + C_e \f$
    Eigen::MatrixXd matKs;
    mismatchCov_ -> giveCovariance( hyperParam, matKs );

    //! compute \f$ ( C_d + C_e )^{-1} \f$
    Eigen::FullPivLU < Eigen::MatrixXd > luKs( matKs );
    FTL_VERIFY( luKs.isInvertible( ) );
    Eigen::MatrixXd matKsInv = luKs.inverse( );

    //! get \f$ C_u = P A^{-1} C_r  A^{-T} P^T \f$
    Eigen::MatrixXd covU;
    forwardGauss_ -> giveProjectedCovarU( covU );

    //! get rho and numReadings
    const auto it = hyperParam.find( hypIdPre_ + hypId_ );
    FTL_VERIFY( it != hyperParam.end( ) );
    const double rho = it->second;
    const unsigned numReadings = sensorValues_.cols( );

    //! compute matrix \f$ S = ( (1 / \rho^2 n_o) K_s ) + P C_u P^T \f$
    Eigen::MatrixXd matSinv;
    {
        Eigen::MatrixXd matS = ( matKs / ( rho * rho * numReadings ) ) + covU;
        Eigen::FullPivLU<Eigen::MatrixXd> luMatS( matS );
        FTL_VERIFY( luMatS.isInvertible( ) );
        matSinv = luMatS.inverse( );
    }

    //! compute unprojected \f$ C_u \f$ using identity projection matrix 
    Eigen::MatrixXd covUNotProjected;
    Eigen::MatrixXd matIdentity = Eigen::MatrixXd::Identity( pMat_.cols(), pMat_.cols() );
    forwardGauss_ -> setProjectionMatrix( matIdentity );
    forwardGauss_ -> updateSolution( hyperParam );
    forwardGauss_ -> giveProjectedCovarU( covUNotProjected );

    //! compute matrix \f$ D = I - ( C_u P^T S^{-1} P ) \f$
    Eigen::MatrixXd matD;
    matD = matIdentity - ( covUNotProjected * pMat_.transpose( ) * matSinv * pMat_ );
    
    //! compute matrix \f$ C_{u|y} = D C_u \f$
    feCovarCondData = matD * covUNotProjected;

    //! compute unprojected \f$ \overline u \f$ using identity projection matrix
    Eigen::VectorXd feMeanNotProjected;
    forwardGauss_ -> giveProjectedMeanU( feMeanNotProjected );

    //! compute \f$ \sum_k y_k \f$
    Eigen::VectorXd sumSensorValues = Eigen::VectorXd::Zero( sensorValues_.rows( ) );
    for ( auto i = 0; i < numReadings; ++i ) {
        sumSensorValues += sensorValues_.col( i );
    }

    //! compute the FE mean conditioned on data \f$ \overline u_{u|y} \f$
    const Eigen::VectorXd vecA = rho * covUNotProjected * pMat_.transpose( ) * matKsInv * sumSensorValues;
    const Eigen::VectorXd vecB = feMeanNotProjected;
    feMeanCondData =  matD * ( vecA + vecB );
    
    return;
}


