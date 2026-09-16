// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//
//                 Computational Structural Mechanics Lab
//                         University of Cambridge
//
// This software is copyrighted and all rights are retained by the University of
// Cambridge. All use, disclosure, and/or reproduction of any part not expressly
// authorised by the University of Cambridge is prohibited. (C) 2010 --
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

//! @file GaussianProcess.hpp

#ifndef statX_infer_gaussianprocess_h
#define statX_infer_gaussianprocess_h

#include <vector>
#include <string>
#include <map>

#include <Eigen/Core>
#include <Eigen/Dense>

#include <corlib/verify.hpp>
#include <corlib/eigenX.hpp>

#include <statX/infer/ExpKernel.hpp>

//------------------------------------------------------------------------------
namespace statX {
    namespace infer {
        template <typename KERNEL, unsigned NUMPTS>
        class GaussianProcess;

        namespace eigenX = corlib::eigenX;
    }
}

//------------------------------------------------------------------------------
/** \brief Gaussian Process regression for given locations and observations
 *  \details
 *
 *   Implemented according to page 19 of the book 'Gaussian Processes for Machine
 *   Learning'
 *
 *   \tparam KERNEL    kernel type to be used for the covariance
 *   \tparam NUMPTS    number of observation points
 *
 */
template <typename KERNEL, unsigned NUMPTS>
class statX::infer::GaussianProcess
{
public:
    //! template parameters and values
    typedef KERNEL          Kernel;
    const unsigned numPts = NUMPTS;

    //! public convenience typedefs
    typedef typename Kernel::HpsValueMap HpsValueMap;

public:
    //! constructor
    GaussianProcess ()
    {
        //! hyperparameter names to be used in HpsValueMap
        //! regularisation parameter
        hypIds_.push_back( "eps" );

        //! create names of the observation variables w0, w1, w2, ...
        for ( unsigned i = 0; i < numPts; ++ i ) {
            hypIds_.push_back( "w" + std::to_string( i ) );
        }

        kernel_.setHypIdPre( "gaussianProcess." );
    }

    //! set preamble for hyperparameter names
    void setHypIdPre ( const std::string & preamble )
    {
        hypIdPre_ = preamble;
    }

    //! set the locations
    void createGaussianProcess ( HpsValueMap & hParams, const Eigen::MatrixXd & xMat )
    {
        //! store observation values in hParams_
        hParams_ = hParams;

        FTL_VERIFY( xMat.rows( ) == numPts );

        //! store observation points in xMat_
        xMat_ = xMat;

        //! create Gaussian process to be yet conditioned on observation values
        this->fit( );

        return;
    }

    //! set the observation values
    void setObservationValues ( const HpsValueMap & hParams )
    {
        //! extract the observation values from the map
        yVec_.resize( numPts );

        for ( unsigned i = 0; i < numPts; ++ i ) {
            // note: hypIds_[ 0 ] = "eps"
            const auto it = hParams.find( hypIdPre_ + hypIds_[ i + 1 ] );
            FTL_VERIFY( it != hParams.end( ) );

            yVec_( i ) = it->second;
        }

        return;
    }
    
    //! predict the function value mean at the given location
    double predictMean ( const Eigen::VectorXd & xStar ) const
    {
        //! check the compatibility of test location
        FTL_VERIFY( xStar.rows( ) == xMat_.cols( ) );

        //! compute the k* = k(x, x*), the covariance between input and test locations
        Eigen::VectorXd kStarVec( numPts );

        for ( unsigned i = 0; i < numPts; ++ i ) {
            kStarVec( i ) = kernel_.giveValue( hParams_, xMat_.row( i ), xStar );
        }

        FTL_VERIFY( yVec_.size( ) == numPts );

        //! compute the mean at the given location x*
        //  k*^T (k(x, x) + sigma_n^2 I)^-1 y
        const double predictedMean = kStarVec.dot( cyMatSolver_.solve( yVec_ ) );

        return predictedMean;
    }

    //! predict the function value variance at the given location
    double predictVariance ( const Eigen::VectorXd & xStar ) const
    {
        //! check the compatibility of test location
        FTL_VERIFY( xStar.rows( ) == xMat_.cols( ) );

        //! compute the k* = k(x, x*), the covariance between input and test locations
        Eigen::VectorXd kStarVec( numPts );

        for ( unsigned i = 0; i < numPts; ++ i ) {
            kStarVec( i ) = kernel_.giveValue( hParams_, xMat_.row( i ), xStar );
        }

        //! compute the variance at the given location x*
        //  k(x*, x*) - k*^T (k(x, x) + sigma_n^2 I)^{-1} k*
        const double predictedVar = kernel_.giveValue( hParams_, xStar, xStar )
                - kStarVec.dot( cyMatSolver_.solve( kStarVec ) );

        return predictedVar;
    }

    //! determine the covariance matrix for the given test locations xStars
    void predictCovariance ( const Eigen::MatrixXd & xStars, Eigen::MatrixXd & covMat )
    {
        //! check the compatibility of test location
        FTL_VERIFY( xStars.cols( ) == xMat_.cols( ) );

        //! number of test locations
        const unsigned numStars = xStars.rows( );

        //! compute the k* = k(x, x*), the covariance between input and test locations
        Eigen::MatrixXd kInputStar( numPts, numStars );
        for ( unsigned i = 0; i < numPts; ++ i ) {
            for ( unsigned j = 0; j < numStars; ++ j ) {
                kInputStar( i, j ) = kernel_.giveValue( hParams_, xMat_.row( i ), xStars.row( j ) );
            }
        }

        //! compute the kStarStar = k(x*, x*), the covariance computed at the test locations
        Eigen::MatrixXd kStarStar( numStars, numStars );
        for ( unsigned i = 0; i < numStars; ++ i ) {
            for ( unsigned j = 0; j < numStars; ++ j ) {
                kStarStar( i, j ) = kernel_.giveValue( hParams_, xStars.row( i ), xStars.row( j ) );
            }
        }

        //! compute the covariance conditioned on the observations
        //  k(x*, x*) - k*^T (k(x, x) + sigma_n^2 I)^{-1} k*
        covMat = kStarStar - kInputStar.transpose( ) * ( cyMatSolver_.solve( kInputStar ) );

        return;
    }

    //! give the marginal log likelihood \f$ \log p(y) \f$
    double giveLogLikelihood () const
    {
        FTL_VERIFY( yVec_.size( ) == numPts );

        //! tmp1 = -0.5 * y^T * (k(x, x) + sigma_n^2 I)^-1 y
        const double tmp1 = -0.5 * yVec_.dot( cyMatSolver_.solve( yVec_ ) );

        //! tmp2 = -log( det(k(x, x) + sigma_n^2 I) )
        const Eigen::MatrixXd yCovLMat = cyMatSolver_.matrixL( );
        const Eigen::ArrayXd yVarLMatDiag = yCovLMat.diagonal( );
        const double tmp2 = -yVarLMatDiag.log( ).sum( );

        //! tmp3 = -0.5 * n * log( 2 pi )
        const double tmp3 = -0.5 * xMat_.rows( ) * std::log( 2.0 * M_PI );

        //! log p(y|x) = tmp1 + tmp2 + tmp3
        return tmp1 + tmp2 + tmp3;
    }

private:
    //! Gaussian process regression
    void fit ()
    {
        //! compute k(x, x), the kernel matrix corresponding to the input locations
        Eigen::MatrixXd kMat( numPts, numPts );
        for ( unsigned i = 0; i < numPts; ++ i ) {
            kMat( i, i ) = kernel_.giveValue( hParams_, xMat_.row( i ),
                                              xMat_.row( i ) );
            for ( unsigned j = i + 1; j < numPts; ++ j ) {
                kMat( i, j ) = kernel_.giveValue( hParams_, xMat_.row( i ),
                                                  xMat_.row( j ) );
                kMat( j, i ) = kMat( i, j );
            }
        }

        //! compute k(x, x) + sigma_n^2 I
        const auto it = hParams_.find( hypIdPre_ + "eps" );
        FTL_VERIFY( it != hParams_.end( ) );
        const double eps = it->second;

        Eigen::MatrixXd yVarMat = kMat
                + eps * eps * Eigen::MatrixXd::Identity( numPts, numPts );

        //! compute the factorised y covariance matrix
        cyMatSolver_.compute( yVarMat );

        //! if the Matrix is not invertible, the program terminates
        FTL_VERIFY( cyMatSolver_.info( ) == Eigen::Success );

        return;
    }

private:
    //! covariance kernel
    Kernel                      kernel_;

    //! local hyperparameters 
    std::string                 hypIdPre_;
    std::vector<std::string>    hypIds_;

    //! the parameter map
    HpsValueMap                 hParams_;

    //! inputs and observations
    Eigen::MatrixXd             xMat_;
    Eigen::VectorXd             yVec_;

    //! auxiliary containers
    Eigen::LLT<Eigen::MatrixXd> cyMatSolver_;
};
//------------------------------------------------------------------------------

#endif
