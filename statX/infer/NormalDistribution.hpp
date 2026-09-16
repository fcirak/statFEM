// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//
//                 Computational Structural Mechanics Lab
//                         University of Cambridge
//
// This software is copyrighted and all rights are retained by the University of
// Cambridge. All use, disclosure, and/or reproduction of any part not expressly
// authorised by the University of Cambridge is prohibited. (C) 2010 --
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

//! @file NormalDistribution.hpp

#ifndef statX_infer_normaldistribution_h
#define statX_infer_normaldistribution_h

#include <random>
#include <Eigen/Dense>

#include <corlib/verify.hpp>

#include <statX/infer/NormalSampler.hpp>

//------------------------------------------------------------------------------
namespace statX {
    namespace infer {
        template <const unsigned DIM>
        class NormalDistribution;

        namespace eigenX = corlib::eigenX;
    }
}

//! Multivariate normal distribution
template <const unsigned DIM>
class statX::infer::NormalDistribution
{
public:
    static constexpr unsigned dim = DIM;
    typedef eigenX::VectorSd<dim> Vec;
    typedef eigenX::MatrixSd<dim> Mat;

    //! Constructor with mean vector of size dim and covariance matrix of size dim x dim
    NormalDistribution ( const Vec & mean, const Mat & covariance )
            : mean_( mean ), cov_( covariance )
    { }

    //! Evaluate PDF value at position x
    double givePDF( const Vec & x )
    {
        const Vec diff = x - mean_;

        Eigen::FullPivLU < Eigen::MatrixXd > covLU( cov_ );
        FTL_VERIFY( covLU.isInvertible( ) );
        const double tmp = diff.transpose( ) * covLU.inverse( ) * diff;
        const double val = 1. / ( std::sqrt( std::pow( 2. * M_PI, dim ) * covLU.determinant( ) ) ) *
                            std::exp( -0.5 * tmp );
        return val;
    }

    //! Give random samples with size numReadings
    void giveSamples( Eigen::MatrixXd & sampleMat, const unsigned numReadings,
                      const unsigned seed = std::random_device( )( ) )
    {
        statX::infer::NormalSampler normalSampler( seed );

        sampleMat.resize( dim, numReadings );
        for ( unsigned i = 0; i < numReadings; ++ i ) {
            Eigen::VectorXd sampleVec;
            normalSampler.giveSample( mean_, cov_, sampleVec );
            sampleMat.col( i ) = sampleVec;
        }

        return;
    }


private:
    const Vec mean_;
    const Mat cov_;

};

//! Univariate normal distribution
template <>
class statX::infer::NormalDistribution<1>
{
public:
    static constexpr unsigned dim = 1;
    typedef eigenX::VectorSd<1> Vec;
    typedef eigenX::MatrixSd<1> Mat;

    //! Constructor with mean vector of size 1 and variance matrix of size 1 x 1
    NormalDistribution ( const Vec & mean, const Mat & variance )
            : mean_( mean ), cov_( variance )
    { }

    //! Evaluate PDF value at position x
    double givePDF ( const Vec & x )
    {
        const double y = x( 0 );
        const double mu = mean_( 0 );
        const double sigma = std::sqrt( cov_( 0, 0 ) );

        const double val = 1. / ( std::sqrt( 2. * M_PI ) * sigma )
                * std::exp( -( y - mu ) * ( y - mu ) / ( 2. * sigma * sigma ) );
        return val;
    }

    //! Give random samples with size numReadings
    void giveSamples( Eigen::MatrixXd & sampleMat, const unsigned numReadings = 1,
                      const unsigned seed = std::random_device( )( ) )
    {
        statX::infer::NormalSampler normalSampler( seed );

        sampleMat.resize( 1, numReadings );

        for ( unsigned i = 0; i < numReadings; ++ i ) {
            Eigen::VectorXd sampleVec;
            normalSampler.giveSample( mean_, cov_, sampleVec );
            sampleMat.col( i ) = sampleVec;
        }

        return;
    }

private:
    const Vec mean_;
    const Mat cov_;
};

#endif
