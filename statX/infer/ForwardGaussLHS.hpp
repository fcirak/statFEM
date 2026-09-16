// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//
//                 Computational Structural Mechanics Lab
//                         University of Cambridge
//
// This software is copyrighted and all rights are retained by the University of
// Cambridge. All use, disclosure, and/or reproduction of any part not expressly
// authorised by the University of Cambridge is prohibited. (C) 2010 --
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

//! @file ForwardGaussLHS.hpp

#ifndef statX_infer_forwardgausslhs_h
#define statX_infer_forwardgausslhs_h

#include <string>
#include <functional>

#include <Eigen/Core>

#include <corlib/SystemSolveEigenSparse.hpp>
#include <corlib/eigenX.hpp>
#include <corlib/verify.hpp>

//------------------------------------------------------------------------------
namespace statX {
    namespace infer {
        class ForwardGaussLHS;
    }
}

//------------------------------------------------------------------------------
/** \brief Wrapper for the forward problem with only LHS uncertainty.
 *  \details
 *
 *  Provides so far only the mean \f$ u = A(w)^-1 f \f$ and will in future
 *  provide the density p(u)
 *
 */
class statX::infer::ForwardGaussLHS
{
public:
    //! public convenience typedefs
    typedef std::map<std::string, double> HpsValueMap;
    typedef std::function<void ( HpsValueMap &, Eigen::VectorXd & )> SolveFunc;

public:
    /** constructor with references to
     *  solveFunc: function to compute u for a given w
     *  pMat     : projection matrix
     */
    ForwardGaussLHS ( const SolveFunc & solveFunc, const Eigen::MatrixXd & pMat )
            : solveFunc_( solveFunc ), pMat_( pMat )
    { }

    //! set preamble for hyperparameter names
    void setHypIdPre ( const std::string & preamble )
    {
        hypIdPre_ = preamble;
    }

    //! give the projected mean \f$ P A^{-1} f \f$
    void giveProjectedMeanU ( Eigen::VectorXd & projMean ) const
    {
        FTL_VERIFY( feMean_.size( ) );
        projMean = feMean_;

        return;
    }

    //! this method is to be yet implemented
    void giveProjectedCovarU ( Eigen::MatrixXd & projCovmat ) const
    {
        const unsigned numSensors = pMat_.rows( );
        projCovmat = Eigen::MatrixXd::Zero( numSensors, numSensors );

        return;
    }

    //! compute mean \f$ P A(w)^{-1} f \f$ for given hyperparameters w
    void updateSolution ( HpsValueMap & hParam )
    {
        //! compute the finite element solution
        Eigen::VectorXd feSoln;
        solveFunc_( hParam, feSoln );

        //! project the finite element solution
        feMean_ = pMat_ * feSoln;

        return;
    }

private:

    //! function to compute u for a given w
    const SolveFunc           solveFunc_;

    //! preamble for local parameter names
    std::string               hypIdPre_;

    //! projected mean
    Eigen::VectorXd           feMean_;

    //! projection matrix P
    const Eigen::MatrixXd     pMat_;
};

#endif
