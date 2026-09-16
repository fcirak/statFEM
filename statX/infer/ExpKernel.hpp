// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//
//                 Computational Structural Mechanics Lab
//                         University of Cambridge
//
// This software is copyrighted and all rights are retained by the University of
// Cambridge. All use, disclosure, and/or reproduction of any part not expressly
// authorised by the University of Cambridge is prohibited. (C) 2010 --
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

//! @fil ExpKernel.hpp

#ifndef statX_infer_expkernel_h
#define statX_infer_expkernel_h

#include <cmath>
#include <map>
#include <string>
#include <array>

#include <Eigen/Core>

#include <corlib/eigenX.hpp>
#include <corlib/verify.hpp>

//------------------------------------------------------------------------------
namespace statX {
    namespace infer {
        template <typename VECDIM>
        class ExpKernel;
    }
}

//------------------------------------------------------------------------------
/** \brief Wrapper class for the exponential covariance kernel
 *  \details
 *
 *  Provides kernel values and derivatives for given hyperParameter values (sigma
 *  and length)
 *
 *  \tparam VECDIM     coordinate container
 *
 */
template <typename VECDIM>
class statX::infer::ExpKernel
{
public:
    //! template parameters
    typedef VECDIM VecDim;

    //! number of hyperparameters
    static constexpr int numHps = 2;

    //! public convenience typedefs
    typedef std::map<std::string, double>       HpsValueMap;
    typedef std::array<std::string, numHps>     HypIdsArray;

public:
    //! constructor
    ExpKernel ()
    {
        //! hyperparameter names to be used in HpsValueMap
        hypIds_[ 0 ] = "sigma";
        hypIds_[ 1 ] = "length";
    }

    //! set preamble for hyperparameter names
    void setHypIdPre ( const std::string & preamble )
    {
        hypIdPre_ = preamble;
    }

    //! give kernel value for given hyperparameters
    double giveValue ( const HpsValueMap & hps, const VecDim & nCoor1,
                       const VecDim & nCoor2 ) const
    {
        //! find sigma value
        auto it = hps.find( hypIdPre_ + hypIds_[ 0 ] );
        FTL_VERIFY( it != hps.end( ) );
        const double sigma = it->second;

        //! find length value
        it = hps.find( hypIdPre_ + hypIds_[ 1 ] );
        FTL_VERIFY( it != hps.end( ) );
        const double length = it->second;

        const VecDim dist = nCoor1 - nCoor2;
        const double val = sigma * sigma
                * std::exp( -dist.squaredNorm( ) / ( 2. * length * length ) );

        return val;
    }

    //! give kernel derivatives for given hyperparameter map
    void giveDerivHParam ( const HpsValueMap & hps, const VecDim & nCoor1,
                           const VecDim & nCoor2, HpsValueMap & deriv ) const
    {
        //! find sigma value
        auto it = hps.find( hypIdPre_ + hypIds_[ 0 ] );
        FTL_VERIFY( it != hps.end( ) );
        const double sigma = it->second;

        //! find length value
        it = hps.find( hypIdPre_ + hypIds_[ 1 ] );
        FTL_VERIFY( it != hps.end( ) );
        const double length = it->second;

        const VecDim dist = nCoor1 - nCoor2;

        const double value = giveValue( hps, nCoor2, nCoor2 );

        const double dKerDk = 2. * value / sigma;
        const double dKerDl = -value * dist.squaredNorm( )
                / ( length * length * length );

        deriv[ hypIdPre_ + hypIds_[ 0 ] ] = dKerDk;
        deriv[ hypIdPre_ + hypIds_[ 1 ] ] = dKerDl;

        return;
    }

private:
    //! local hyperparameter names
    HypIdsArray     hypIds_;
    std::string     hypIdPre_;

};


#endif
