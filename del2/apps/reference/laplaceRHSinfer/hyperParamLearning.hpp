// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//
//                 Computational Structural Mechanics Lab
//                         University of Cambridge
//
// This software is copyrighted and all rights are retained by the University of
// Cambridge. All use, disclosure, and/or reproduction of any part not expressly
// authorised by the University of Cambridge is prohibited. (C) 2010 --
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

//! @file hyperParamLearning.hpp

#ifndef del2_apps_reference_laplaceRHSinfer_hyperparamlearning_h
#define del2_apps_reference_laplaceRHSinfer_hyperparamlearning_h

#include <vector>
#include <map>
#include <string>
#include <fstream>
#include <cmath>
#include <functional>
#include <utility>
#include <random>
#include <algorithm>

#include <corlib/PropertiesParser.hpp>
#include <corlib/verify.hpp>

#include <statX/infer/LikelihoodKOH.hpp>
#include <statX/infer/MCMC.hpp>

#include "ioHelpers.hpp"

namespace infer  = statX::infer;
namespace eigenX = corlib::eigenX;

//------------------------------------------------------------------------------
class NormalPDF
{
public:
    NormalPDF ( const double mu, const double sigma )
            : mu_( mu ), sigma_( sigma )
    { }

    double giveValue ( const double x ) const
    {
        const double val = 1. / ( std::sqrt( 2. * M_PI ) * sigma_ )
                * std::exp( -( x - mu_ ) * ( x - mu_ ) / ( 2. * sigma_ * sigma_ ) );
        return val;
    }

private:
    const double mu_;
    const double sigma_;

};

//------------------------------------------------------------------------------
template <unsigned dim>
class MultiVariateNormalPDF
{
public:
    typedef eigenX::VectorSd<dim> Vec;
    typedef eigenX::MatrixSd<dim> Mat;

    MultiVariateNormalPDF ( const Vec & mean, const Mat & covariance )
            : mean_( mean ), cov_( covariance )
    { }

    double giveValue ( const Vec & x ) const
    {

        const unsigned d = mean_.rows( );
        const Eigen::VectorXd diff = x - mean_;

        Eigen::FullPivLU < Eigen::MatrixXd > covLU( cov_ );
        FTL_VERIFY( covLU.isInvertible( ) );
        const double tmp = diff.transpose( ) * covLU.inverse( ) * diff;
        const double val = 1. / ( std::sqrt( std::pow( 2. * M_PI, d ) * covLU.determinant( ) ) ) *
                            std::exp( -0.5 * tmp );
        return val;
    }

private:
    const Vec mean_;
    const Mat cov_;
};

//------------------------------------------------------------------------------
template <typename LikelihoodKOH, typename HpsValueMap>
void hyperParamLearning ( LikelihoodKOH & likelihoodKOH,
                          const std::string & inferenceInputProps,
                          HpsValueMap & hParamMap )
{

    //! reading inference parameters
    corlib::PropertiesParser *inferenceProps = new corlib::PropertiesParser;

    //! number of samples and readings
    unsigned numSamples, numReadings, enrichWithPrior;
    inferenceProps -> registerPropertiesVar( "numSamples",       numSamples     );
    inferenceProps -> registerPropertiesVar( "numReadings",      numReadings    );
    inferenceProps -> registerPropertiesVar( "enrichWithPrior",  enrichWithPrior);
    //! output file names
    std::string hpsMeanFile, sampleMCMCFile;
    inferenceProps -> registerPropertiesVar( "sampleMCMCFile",   sampleMCMCFile );
    inferenceProps -> registerPropertiesVar( "hpsMeanFile",      hpsMeanFile    );
    //! MCMC input file names
    std::string priorMeanFile, priorCovarFile, initHPFile, varianceIncFile;
    inferenceProps -> registerPropertiesVar( "priorMeanFile",    priorMeanFile  );
    inferenceProps -> registerPropertiesVar( "priorCovarFile",   priorCovarFile );
    inferenceProps -> registerPropertiesVar( "initHPFile",       initHPFile     );
    inferenceProps -> registerPropertiesVar( "varianceIncFile",  varianceIncFile);

    //! read variables from the inferenceInput.dat file
    std::ifstream inferenceInputFile( inferenceInputProps );
    FTL_VERIFY( inferenceInputFile.is_open( ) );
    inferenceProps -> readValues( inferenceInputFile );

    //! delete properties parser and close the input file
    delete inferenceProps;
    inferenceInputFile.close( );

    //! the free hyperparameters to be sampled
    std::array<std::string, 3> freeHParamIds = { "rho",
                                                 "mismatch.sigma",
                                                 "mismatch.length" };
    const unsigned numAllHParams = hParamMap.size( );
    const unsigned numFreeHParams = freeHParamIds.size( );

    //! define the posterior density p(\theta | y)
    typedef infer::MCMC<numFreeHParams>          MCMC;

    //! sampling the hyperparameters in the log space to ensure non-negativity
    Eigen::VectorXd priorMean( numFreeHParams );
    Eigen::MatrixXd priorCovar;
    priorCovar = Eigen::MatrixXd::Zero( numFreeHParams, numFreeHParams );

    //! get priorMean and priorCovar from parameter map
    for ( unsigned i = 0; i < numFreeHParams; ++ i ) {
        const auto itMean = hParamMap.find( "priorMean." + freeHParamIds[ i ] );
        FTL_VERIFY( itMean != hParamMap.end( ) );
        priorMean( i ) = std::log( itMean->second );

        const auto itCov = hParamMap.find( "priorCov." + freeHParamIds[ i ] );
        FTL_VERIFY( itCov != hParamMap.end( ) );
        priorCovar( i, i ) = itCov->second;
    }

    MultiVariateNormalPDF<numFreeHParams> priors( priorMean, priorCovar );

    typedef eigenX::VectorSd<numFreeHParams>     VecHps;

    //! get hyperparameters from the map and store it in a vector
    auto getHParams = [ & freeHParamIds ] ( HpsValueMap & hpMap, VecHps & freeHParam ) {
        for ( unsigned i = 0; i < freeHParamIds.size( ); ++i ) {
            const auto it = hpMap.find( freeHParamIds[ i ] );
            FTL_VERIFY( it != hpMap.end( ) );
            freeHParam( i ) = it->second;
        }
    };

    //! set hyperparameters from a vector to the map
    auto setHparams = [ & freeHParamIds ] ( HpsValueMap & hpMap, VecHps & freeHParam ) {
        for ( unsigned i = 0; i < freeHParamIds.size( ); ++i ) {
            const auto it = hpMap.find( freeHParamIds[ i ] );
            FTL_VERIFY( it != hpMap.end( ) );
            it->second = freeHParam( i );
        }
    };

    auto posteriorGivenData = [ & likelihoodKOH,  & priors, & hParamMap, & enrichWithPrior,
                                & setHparams ]
                                ( VecHps & hpVec ) {
        VecHps expHpVec;
        for( unsigned i = 0; i < hpVec.size( ); ++i ) {
            expHpVec( i ) = std::exp( hpVec( i ) );
        }
        setHparams( hParamMap, expHpVec );

        //! compute the log-likelihood considering the non-negativity of the hyperparameters
        double likelihood = likelihoodKOH.giveLogLikelihood( hParamMap ) + hpVec.sum( );

        //! add contribution from prior when enrichWithPrior == true
        if ( enrichWithPrior ) {
            likelihood += std::log( priors.giveValue( expHpVec ) );
        }

        return likelihood;
    };

    //! MCMC sampling of the posterior density of the hyperparameters p( w | y )
    unsigned seed( 0 );
    MCMC mcmc( posteriorGivenData, seed );

    //! get the initial sample and proposal proposed increment
    VecHps varianceInc, initHP;
    for ( unsigned i = 0; i < numFreeHParams; ++ i ) {
        varianceInc( i ) = hParamMap[ "propStDev." + freeHParamIds[ i ] ];
        initHP( i ) = std::log( hParamMap[ "init." + freeHParamIds[ i ] ] );
    }

    //mcmc.computePropVariance( propVariance, varianceInc, initHP );
    std::cout << "posterior of the initial value = "
              << posteriorGivenData( initHP ) << std::endl;

    mcmc.sample( numSamples, varianceInc, initHP );

    typedef typename MCMC::VecHparPair VecHPpair;
    std::vector<VecHPpair> samples = mcmc.giveSamples( );

    const unsigned sampleIdBegin = unsigned( samples.size( ) * 3.0 / 10.0 );
    Eigen::VectorXd sampleMean( numFreeHParams );
    sampleMean.setZero( );

    std::ofstream mlikelihoodStream( sampleMCMCFile.c_str( ) );
    FTL_VERIFY( mlikelihoodStream.is_open( ) );

    std::ofstream hpsMeanStream( hpsMeanFile.c_str( ) );
    FTL_VERIFY( hpsMeanStream.is_open( ) );

    for ( unsigned j = sampleIdBegin; j < samples.size( ); ++ j ) {
        const VecHps hparam = samples[ j ].first;
        const double value = samples[ j ].second;

        for ( unsigned k = 0; k < numFreeHParams; ++ k ) {
            sampleMean( k ) += std::exp( hparam( k ) );
            mlikelihoodStream << std::exp( hparam( k ) ) << "  ";
        }
        mlikelihoodStream << value << std::endl;

    }

    mlikelihoodStream.close( );

    for ( unsigned k = 0; k < numFreeHParams; ++ k ) {
        //! mean
        const double mean = sampleMean( k ) / ( samples.size( ) - sampleIdBegin );
        hpsMeanStream << freeHParamIds[ k ] << "  " << mean << std::endl;
        std::cout << mean << std::endl;
    }
    hpsMeanStream.close( );

    //! set the inferred values to hyperparameters map
    for ( unsigned k = 0; k < numFreeHParams; ++ k ) {
        hParamMap[ freeHParamIds[ k ] ] = sampleMean( k )
                / ( (double) samples.size( ) - sampleIdBegin );
    }

}

#endif
