// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//
//                 Computational Structural Mechanics Lab
//                         University of Cambridge
//
// This software is copyrighted and all rights are retained by the University of
// Cambridge. All use, disclosure, and/or reproduction of any part not expressly
// authorised by the University of Cambridge is prohibited. (C) 2010 --
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

//! @file LikelihoodKOH.hpp

#ifndef statX_infer_likelihoodkoh_h
#define statX_infer_likelihoodkoh_h

#include <map>
#include <random>
#include <string>

#include <Eigen/Core>

#include <corlib/eigenX.hpp>

#include <statX/infer/MismatchCovariance.hpp>
#include <statX/infer/ForwardGaussRHS.hpp>
#include <statX/infer/ForwardGaussLHS.hpp>

//------------------------------------------------------------------------------
namespace statX {
    namespace infer {
        template <typename FORWARDGAUSS, typename MISMATCHCOV>
        class LikelihoodKOH;

        namespace eigenX = corlib::eigenX;
    }
}

//------------------------------------------------------------------------------
/** \brief Kennedy and O'Hagan regression with uncertainty only on right hand side
 *  \details
 *
 *  Refer to paper Kennedy, M. & O'Hagan, A. (2011). Bayesian calibration
 *  of computer models. Journal of the Royal Statistical Society, 63(3): 425-464.
 *
 *  \tparam FORWARDGAUSS     forward problem solver providing \f$ p(u) \f$
 *  \tparam MISMATCHCOV     container providing the covariance \$ C_d + C_e\$
 *
 */
template <typename FORWARDGAUSS, typename MISMATCHCOV>
class statX::infer::LikelihoodKOH
{
public:
    //! template parameters
    typedef FORWARDGAUSS    ForwardGauss;
    typedef MISMATCHCOV     MismatchCov;

    //! public convenience typedefs and values
    static constexpr unsigned numHps = 1;
    typedef std::map<std::string, double>              HpsValueMap;
    typedef std::map<std::string, Eigen::MatrixXd>     HpsMatMap;

public:
    //! constructor
    LikelihoodKOH ( ForwardGauss * forwardGauss, MismatchCov * mismatchCov )
            : forwardGauss_( forwardGauss ), mismatchCov_( mismatchCov ),
              hypId_( "rho" )
    {
        // set preamble for sensor parameter
        mismatchCov_ -> setHypIdPre( "mismatch." );
    }

	//! set preamble for hyperparameter names
    void setHypIdPre ( const std::string & preamble )
    {
        hypIdPre_ = preamble;
    }

	//! set sensor readings
    void setSensorValues ( const Eigen::MatrixXd & sensorValues )
    {
        sensorValues_ = sensorValues;
        mismatchScaling_ = Eigen::VectorXd::Ones( sensorValues.cols( ) );
        return;
    }

    //! set scaling for mismatch covariance - for time series data
    void setMismatchScaling ( const Eigen::VectorXd & scaling )
    {
        FTL_VERIFY( mismatchScaling_.size( ) );
        FTL_VERIFY( mismatchScaling_.size( ) == scaling.size( ) );
        mismatchScaling_ = scaling;
    }

    //! set projection matrix
    void setProjectionMatrix ( Eigen::MatrixXd & pMat )
    {
        pMat_ = pMat;
    }

    //! give sample of sensor reading at observation points
    void giveTypicalReadings ( const HpsValueMap & hyperParam, unsigned numReadings,
                               Eigen::MatrixXd & yTypicalMat,
                               unsigned seed = std::random_device( )( ) );

    //! give log likelihood for given hyperparameters
    double giveLogLikelihood ( const HpsValueMap & hyperParam );

    //! give log likelihood for given hyperparameters using time-series data
    //! with mean of p(u) varies accross time
    double giveLogLikelihoodTimeSeries ( const HpsValueMap & hyperParam );


    //! give log likelihood gradient with respect to hyperparameters
    void giveLogLikelihoodGrad ( const HpsValueMap & hyperParam,
                                 HpsValueMap & mllDeriv );

    //! give FE mean and covariance, i.e. \f$ p(u|y) \f$ conditioned on sensor data
    void giveFEconditionedData ( const HpsValueMap & hyperParam,
                                 Eigen::VectorXd & feMeanCondData,
                                 Eigen::MatrixXd & feCovarCondData ) const;

    //! give projected FE mean and covariance, i.e. \f$ p(u|y) \f$ conditioned on sensor data
    void giveProjFEconditionedData ( const HpsValueMap & hyperParam,
                                 Eigen::VectorXd & feMeanCondData,
                                 Eigen::MatrixXd & feCovarCondData ) const;

    //! give projected FE mean and covariance, i.e. \f$ p(Pu|y) \f$ conditioned on time-series sensor data
    void giveProjFEconditionedDataTimeSeries ( const HpsValueMap & hyperParam,
                                          Eigen::VectorXd & feMeanCondData,
                                          Eigen::MatrixXd & feCovarCondData,
                                          unsigned index ) const;

    //! give FE mean and covariance, i.e. \f$ p(u|y) \f$ conditioned on sensor data
    //! with Woodbury matrix expansion
    void giveFEconditionedDataWoodbury ( const HpsValueMap & hyperParam,
                                         Eigen::VectorXd & feMeanCondData,
                                         Eigen::MatrixXd & feCovarCondData ) const;

private:
    //! method to compute the total covariance
    void computeTotalCovariance_ ( const HpsValueMap & hyperParam,
                                   Eigen::MatrixXd & covTotal );

private:
    //! local parameter names
    const std::string      hypId_;
    std::string            hypIdPre_;

    //! sensor values
    Eigen::MatrixXd        sensorValues_;

    //! mismatch covariance and forwardGauss
    MismatchCov          * const mismatchCov_;
    ForwardGauss         * const forwardGauss_;

    //! the projection matrix
    Eigen::MatrixXd        pMat_;

    //! mismatch scaling
    Eigen::VectorXd        mismatchScaling_;
};

//------------------------------------------------------------------------------
#include "LikelihoodKOH.ipp"

#endif
