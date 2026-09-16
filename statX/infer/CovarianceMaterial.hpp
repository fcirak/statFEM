// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//
//                 Computational Structural Mechanics Lab
//                         University of Cambridge
//
// This software is copyrighted and all rights are retained by the University of
// Cambridge. All use, disclosure, and/or reproduction of any part not expressly
// authorised by the University of Cambridge is prohibited. (C) 2010 --
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

//! @file CovarianceMaterial.hpp

#ifndef statX_infer_covariancematerial_h
#define statX_infer_covariancematerial_h

#include <vector>
#include <string>
#include <map>

#include <Eigen/Core>

#include <corlib/verify.hpp>

//------------------------------------------------------------------------------
namespace statX {
    namespace infer {
        template <unsigned DOF, typename KERNEL>
        class CovarianceMaterial;
    }
}

//------------------------------------------------------------------------------
/** \brief Computes the covariance matrix \f$ C_{\kappa} \f$ at element barycentres
 *  \details
 *
 *  Computes the covariance matrix \f$ C_{\kappa} \f$ and
 *  \f[
 *    (C_{\kappa})_{ij} = c_{\kappa}(\vec x_i , \, \vec x_j)
 *  \f]
 *  where x_i and x_j are the barycentres of elements i and j
 *
 *  \tparam DOF    degress of freedom per node
 *  \tparam KERNEL a kernel like the exponential kernel ExpKernel
 */
template <unsigned DOF, typename KERNEL>
class statX::infer::CovarianceMaterial
{
public:
    //! template parameters and values
    typedef KERNEL Kernel;
    static const unsigned dof = DOF;

    //! public convenience typedefs
    typedef std::map<std::string, double>      HpsValueMap;
    typedef typename Kernel::VecDim            VecDim;

private:
    //! private convenience typedefs
    typedef std::vector<VecDim> CoordVec_;
    static const unsigned dim_ = VecDim::RowsAtCompileTime;

public:
    //! constructor with number of elements
    CovarianceMaterial ( )
    {
        kernel_.setHypIdPre( "material." );
    }

    //! set preamble for hyperparameter names
    void setHypIdPre ( const std::string & preamble )
    {
        hypIdPre_ = preamble;
    }

    //! insert element centroid (eCoords) one by one into the global vector coords_
    //! the numbering is according to element numbering
    void insertToCoordinates ( const VecDim &elemCoords )
    {
        for ( unsigned i = 0; i < dof; i ++ ) {
            coords_.push_back( elemCoords );
        }
        return;
    }

    //! give the material covariance \f$ \vec C_{\kappa} \f$
    void giveMaterialCovariance ( const HpsValueMap & params,
                                  Eigen::MatrixXd & materialCovar ) const
    {
        FTL_VERIFY( coords_.size( ) );
        const unsigned numElements = coords_.size( );
        materialCovar.resize( numElements, numElements );
        for ( unsigned k = 0; k < numElements; ++ k ) {
            for ( unsigned l = 0; l < numElements; ++ l ) {

                const double kval = kernel_.giveValue( params,
                                                       coords_[ k ],
                                                       coords_[ l ] );
                materialCovar( k, l ) = kval;
            }
        }

        return;
    }

private:
    //! preamble for local parameter names
    std::string      hypIdPre_;

    //! covariance kernel
    Kernel           kernel_;

    //! element centroids
    CoordVec_        coords_;
};

#endif
