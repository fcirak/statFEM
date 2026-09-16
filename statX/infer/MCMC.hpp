// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//
//                 Computational Structural Mechanics Lab
//                         University of Cambridge
//
// This software is copyrighted and all rights are retained by the University of
// Cambridge. All use, disclosure, and/or reproduction of any part not expressly
// authorised by the University of Cambridge is prohibited. (C) 2010 --
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

//! @file MCMC.hpp

#ifndef statX_infer_mcmc_h
#define statX_infer_mcmc_h

#include <cmath>
#include <utility>
#include <vector>
#include <random>
#include <functional>
#include <algorithm>

#include <boost/random/mersenne_twister.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/uniform_real_distribution.hpp>
#include <boost/random/variate_generator.hpp>

#include <corlib/eigenX.hpp>
#include <corlib/verify.hpp>

//------------------------------------------------------------------------------
namespace statX {
    namespace infer {
        template <unsigned NUMHPS>
        class MCMC;
    }
}

//------------------------------------------------------------------------------
/** \brief Sample from a log (pdf) given as a function PdfFunc
 *  \details
 */
template <unsigned NUMHPS>
class statX::infer::MCMC
{
public:
    //! template values
    static constexpr unsigned numHps = NUMHPS;

    //! public convenience typedefs
    typedef corlib::eigenX::VectorSd<numHps>               VecHpar;
    typedef std::pair<VecHpar, double>                     VecHparPair;
    typedef std::function<double ( VecHpar & )>            PdfFunc;
    //! boost statistics tools
    typedef boost::random::mt19937 mt19937;
    typedef boost::random::normal_distribution<>           NormalDist;
    typedef boost::random::variate_generator< mt19937&, NormalDist >
                                                           NormalGenerator;
    typedef boost::random::uniform_real_distribution<>     UniformDist;
    typedef boost::random::variate_generator< mt19937&, UniformDist >
                                                           UniformGenerator;


public:
    //! constructor
    MCMC ( const PdfFunc & pdfFunc, const unsigned seed = std::random_device( )( ) )
            : pdfFunc_ { pdfFunc }, seed_ { seed }, acceptanceRatio_ { 0. }
    {}

    /** Method to set proposal variance for the distribution of x
     *  This will iteratively call sample() and adjust the variance vector in proportions to
     *  its initialisation values. i.e. if you want larger variance in one direction,
     *  then set initial variance larger in that direction
     *  https://m-clark.github.io/docs/ld_mcmc/index_onepage.html#am
     *  BÉDARD, M. (2008). Optimal acceptance rates for Metropolis algorithms:
     *  Moving beyond 0.234. Stochastic Process. Appl. 118 2198–2222. MR2474348
     **/
    VecHpar computePropVariance ( VecHpar & propVariance,
                                  const VecHpar & varianceInc,
                                  const VecHpar & xInitial = VecHpar::Zero( ),
                                  const unsigned maxStep = 10,
                                  unsigned nSamples = 1000 )
    {
        double accOpt, accTol;

        // some empirical parameters
        if ( numHps < 5 ) {
            accOpt = 0.44;
            accTol = 0.1;
        }
        else {
            accOpt = 0.23;
            accTol = 0.05;
        }

        bool optimized = false;
        unsigned stepNum = 0;

        while ( ( not optimized ) && ( stepNum < maxStep ) ) {

            // user defined optimisation sample size
            sample( nSamples, propVariance, xInitial );

            // test if above or under acceptance limits, if not optimised yet, apply
            // one-sided bisection method
            if ( acceptanceRatio_ > ( accTol + accOpt ) ) {
                propVariance += varianceInc;
            }
            else if ( acceptanceRatio_ < ( accOpt - accTol ) ) {
                propVariance += 0.5 * varianceInc;
            }
            else {
                //condition for having optimized variance
                optimized = true;
            }
            stepNum ++;
        }

        return propVariance;
    }

    //! function to create specified number of randomly sampled points
    void sample ( unsigned nSamples, const VecHpar & propVariance,
                  const VecHpar & xInitial = VecHpar::Zero( ) )
    {
        samples_.clear( );

        mt19937 engine( seed_ );

        VecHpar xCurr = xInitial;
        double yCurr = pdfFunc_( xCurr );

        int numAccepted = 0;

        for ( unsigned i = 0; i < nSamples; ++ i ) {

            VecHpar xProp = transitionPDF_( xCurr, propVariance, engine );

            const double yProp = pdfFunc_( xProp );
            const double ratio = std::min( 0.0, yProp - yCurr );

            const double randNum = 
                UniformGenerator ( engine, UniformDist(0, 1) )( );

            if ( std::log( randNum ) <= ratio ) {
                xCurr = xProp;
                yCurr = yProp;

                ++ numAccepted;
            }

            samples_.push_back( std::make_pair( xCurr, yCurr ) );
        }

        acceptanceRatio_ = numAccepted / static_cast<double>( samples_.size( ) );

        return;
    }

    //! give sample with maximum value
    VecHparPair giveMax () const
    {
        FTL_VERIFY( not samples_.empty( ) );

        auto maxPair = std::max_element( samples_.begin( ), samples_.end( ),
                [] ( const VecHparPair & v1,
                     const VecHparPair & v2 ) { 
                     return v1.second < v2.second; } );
        return *maxPair;
    }

    //! give the acceptance ratio
    double giveAcceptanceRatio () const
    {
        return acceptanceRatio_;
    }

    //! give the samples
    std::vector<VecHparPair> giveSamples () const
    {
        return samples_;
    }

private:
    //! private method to evaluate the transition kernel
    VecHpar transitionPDF_ ( const VecHpar & x, const VecHpar & propVariance,
                            mt19937 & randomEngine )
    {
        VecHpar xProp;

        for ( unsigned int i = 0; i < numHps; ++ i ) {
            xProp[ i ] = 
                NormalGenerator ( randomEngine, NormalDist(x[ i ], propVariance[ i ]) )( );
        }

        return xProp;
    }


private:
    //! distribution to be sampled
    const PdfFunc     pdfFunc_;

    //! random seed
    const unsigned    seed_;

    //! acceptance ratio of the MCMC
    double      acceptanceRatio_;

    //! storage for the samples
    std::vector<VecHparPair>     samples_;
};


#endif
