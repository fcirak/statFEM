// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//
//                 Computational Structural Mechanics Lab
//                         University of Cambridge
//
// This software is copyrighted and all rights are retained by the University of
// Cambridge. All use, disclosure, and/or reproduction of any part not expressly
// authorised by the University of Cambridge is prohibited. (C) 2010 --
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

//! @file MismatchCovariance.hpp

#ifndef statX_infer_mismatchcovariance_h
#define statX_infer_mismatchcovariance_h

#include <vector>
#include <map>
#include <iterator>

#include <Eigen/Core>

#include <corlib/eigenX.hpp>
#include <statX/infer/ExpKernel.hpp>

//------------------------------------------------------------------------------
namespace statX {
    namespace infer {
        template <typename KERNEL>
        class MismatchCovariance;

        namespace eigenX = corlib::eigenX;
    }
}

//------------------------------------------------------------------------------
/** \brief Square exponential mismatch and sensor covariance
 *  \details
 *
 *  Wrapper for computing the mismatch and sensor covariance \f$ C_d + C_e \f$.
 *
 *  \tparam KERNEL     type of covariance kernel to be used
 *
 */
template <typename KERNEL>
class statX::infer::MismatchCovariance
{
public:
    //! template parameters
    typedef KERNEL                      Kernel;

    //! public convenience typedefs and values
    static constexpr unsigned numHps = 1;
    typedef typename Kernel::VecDim                VecDim;
    typedef std::map<std::string, double>          HpsValueMap;
    typedef std::map<std::string, Eigen::MatrixXd> HpsMatMap;

public:
    /** constructor with references to
      *  Solvefunc: function to compute u for a given w
      *  PMat     : projection matrix
      */
    MismatchCovariance ( const std::vector<VecDim> & sensorCoords )
            : sensorCoords_( sensorCoords ), hypId_ ( "sensorEps" )
    {
        //! define the preamble for mismatch hyperparameters
        kernel_.setHypIdPre( "mismatch." );
    }

    //! set preamble for sensor parameter
    void setHypIdPre ( const std::string & preamble )
    {
        hypIdPre_ = preamble;
    }

    //! give mismatch covariance matrix for given hyperparameter map
    void giveCovariance ( const HpsValueMap & hyperParams,
                          Eigen::MatrixXd & kMat ) const
    {
        const unsigned numSensors = sensorCoords_.size( );
        kMat.resize( numSensors, numSensors );

        auto it = hyperParams.find( hypIdPre_ + hypId_ );
        FTL_VERIFY( it != hyperParams.end( ) );
        const double eps = it -> second;

        for ( unsigned k = 0; k < numSensors; ++ k ) {
            for ( unsigned l = 0; l < numSensors; ++ l ) {

                const double kval = kernel_.giveValue( hyperParams,
                                                       sensorCoords_[ k ],
                                                       sensorCoords_[ l ] );
                kMat( k, l ) = kval;
            }
        }

        kMat += eps * eps * Eigen::MatrixXd::Identity( numSensors, numSensors );

        return;
    }
	
	//! store the vector containing the derivatives of mismatch covariance
	//! with respect to the defined hyperparameters
    void giveCovarianceDerv ( const HpsValueMap & hyperParams,
                              HpsMatMap & mapHpsDMat ) const
    {
        const unsigned numSensors = sensorCoords_.size( );

        for ( unsigned k = 0; k < numSensors; ++ k ) {
            for ( unsigned l = 0; l < numSensors; ++ l ) {

                HpsValueMap deriv;
                kernel_.giveDerivHParam( hyperParams, sensorCoords_[ k ],
                                         sensorCoords_[ l ], deriv );

                for ( auto it = deriv.begin( ); it != deriv.end( ); ++ it ) {
                    mapHpsDMat[ it->first ]( k, l ) = it->second;
                }
            }
        }
        return;
    }


private:
    //! sensor coordinates
	std::vector<VecDim>     sensorCoords_;

	//! local  parameter names
	const std::string       hypId_;
	std::string             hypIdPre_;

	//! the kernel
	Kernel                  kernel_;

};

#endif 
