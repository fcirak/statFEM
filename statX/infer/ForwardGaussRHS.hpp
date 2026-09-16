// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//
//                 Computational Structural Mechanics Lab
//                         University of Cambridge
//
// This software is copyrighted and all rights are retained by the University of
// Cambridge. All use, disclosure, and/or reproduction of any part not expressly
// authorised by the University of Cambridge is prohibited. (C) 2010 --
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

//! @file ForwardGaussRHS.hpp

#ifndef statX_infer_forwardgaussrhs_h
#define statX_infer_forwardgaussrhs_h

#include <string>
#include <map>

#include <Eigen/Core>
#include <Eigen/Sparse>

#include <corlib/SystemSolveEigenSparse.hpp>
#include <corlib/eigenX.hpp>

#include <statX/infer/ExpKernel.hpp>

//------------------------------------------------------------------------------
namespace statX {
    namespace infer {
        template <typename SOLVER, typename COVARFORCING>
        class ForwardGaussRHS;
    }
}

//------------------------------------------------------------------------------
/** \brief Forward problem with uncertainty only on right hand side
 *  \details
 *
 *  Wrapper for the forward problem \f$ p(u) = N( A^{-1} f, A^{-1} C_f A^{-T} ) \f$
 *  with the density \f$ p(f) = N( A^{-1} f, C_f ) \f$
 *
 */
template <typename SOLVER, typename COVARFORCING>
class statX::infer::ForwardGaussRHS
{
public:
    //! template parameters
    typedef SOLVER                            Solver;
    typedef COVARFORCING                      CovarForcing;

    //! public convenience typedef
    typedef std::map<std::string, double>     HpsValueMap;

public:
    /** constructor with pointers and references to
     * solver      : direct solver, usually corlib::SystemSolveEigenSparse
     * covarForcing: class which gives force covariance \f$ C_f \f$
     * pMat        : projection matrix
     * regEps      : regularisation factor as a fraction of \f$ \sigma_f \f$
     */
    ForwardGaussRHS ( Solver * solver, CovarForcing * covarForcing, const double regEps = 1.e-4 )
            : solver_( solver ), covarForcing_( covarForcing ), regEps_( regEps )
    { }

    //! set preamble for hyperparameter names
    void setHypIdPre ( const std::string & preamble )
    {
        hypIdPre_ = preamble;
    }

    void setProjectionMatrix ( Eigen::MatrixXd & pMat )
    {
        pMat_ = pMat;
    }

    //! give the projected mean \f$ P A^{-1} f \f$
    void giveProjectedMeanU ( Eigen::VectorXd & projMean ) const
    {
        FTL_VERIFY( feMean_.size( ) );
        projMean = feMean_;

        return;
    }

    //! give the projected covariance \f$ P A^{-1} C_f A^{-T} P^T \f$
    void giveProjectedCovarU ( Eigen::MatrixXd & projCovmat ) const
    {
        FTL_VERIFY( feCovMat_.size( ) );
        projCovmat = feCovMat_;

        return;
    }

    //! compute the mean and projected covariance
    void updateSolution ( const HpsValueMap & params )
    {
        //! solve auxiliary problem \f$ A^T * R^T = P^T \f$
        FTL_VERIFY( pMat_.size( ) );
        const unsigned numObservations = pMat_.rows( );
        const unsigned numDofs = pMat_.cols( );

        Eigen::MatrixXd rMat( numObservations, numDofs );

        //! one solve per column of \$ P \f$
        bool reuseFactorisation = false;
        for ( unsigned i = 0; i < numObservations; ++ i ) {
            Eigen::VectorXd tmp;
            solver_->solveSystemTranspose( pMat_.row( i ), tmp,
                                           reuseFactorisation );
            reuseFactorisation = true;
            rMat.row( i ) = tmp;
        }

        //! get the force from solver for all dofs
        std::vector<unsigned> dofIndices( numDofs );
        for ( unsigned i = 0; i < numDofs; ++ i ) {
            dofIndices[ i ] = i;
        }

        Eigen::VectorXd force( numDofs );
        solver_->giveForce( dofIndices, force );

        //! compute projected mean of finite element solution \f$ P A^{-1} f \f$
        feMean_ = rMat * force;

        //! compute covariance of FE solution and store it in feCovMat_
        //! \f$ P A^{-1}  C_f A^{-T} P  \f$ for given hyperparameters
        covarForcing_ -> giveProjectedCovar( params, rMat, feCovMat_ );

        //! compute regularisation term
        auto it = params.find( "forcing.sigma" );
        FTL_VERIFY( it != params.end( ) );
        const double addDiag = ( it -> second ) * regEps_;

        //! add regularisation term to \f$ C_u \f$
        feCovMat_ += ( addDiag * addDiag * pMat_ * pMat_.transpose( ) );

        return;
    }

private:
    Solver             * const solver_;
    std::string          hypIdPre_;

    const CovarForcing * const covarForcing_;

    Eigen::VectorXd      feMean_;
    Eigen::MatrixXd      feCovMat_;

    Eigen::MatrixXd      pMat_;

    //! regularisation factor as fraction of \f$ \sigma_f \f$
    const double         regEps_;

};

#endif
