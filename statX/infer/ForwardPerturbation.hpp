// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//
//                 Computational Structural Mechanics Lab
//                         University of Cambridge
//
// This software is copyrighted and all rights are retained by the University of
// Cambridge. All use, disclosure, and/or reproduction of any part not expressly
// authorised by the University of Cambridge is prohibited. (C) 2010 --
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

//! @file ForwardPerturbation.hpp

#ifndef statX_infer_forwardperturbation_h
#define statX_infer_forwardperturbation_h

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
        template <typename SOLVER, typename COVARFORCING, typename COVARMATERIAL>
        class ForwardPerturbation;
    }
}

//------------------------------------------------------------------------------
/** \brief Determine the probability density of the forward solution due to
 *         LHS and RHS uncertainty using first-order perturbation analysis
 *  \details
 *
 *  Wrapper for the forward problem \f$ p(u) = N( \overline{u}, C_u ) \f$
 *  with mean \f$ \overline{u} = A^{-1}(\overline{\kappa}) \overline{f} \f$, and
 *  covariance \f$ C_u = A(\overline{\kappa})^{-1} C_r A(\overline{\kappa})^{-T} +
 *                   \sum_e \sum_d C_{\kappa}_{ed} A(\overline {\kappa})^{-1}
 *                   \frac{\partial A(\overline {\kappa})}{\partial \kappa_e} A(\overline {\kappa})^{-1}
 *                   (C_f + \overline{\vec f} \otimes \overline{f})
 *                   \frac{\partial A(\overline {\kappa})^T}{\partial \kappa_d} A(\overline {\kappa})^{-T}
 *              \f$
 *
 */

template <typename SOLVER, typename COVARFORCING, typename COVARMATERIAL>
class statX::infer::ForwardPerturbation
{
public:
    //! template parameters
    typedef SOLVER                            Solver;
    typedef COVARFORCING                      CovarForcing;
    typedef COVARMATERIAL                     CovarMaterial;

    //! public convenience typedefs
    typedef std::map<std::string, double>              HpsValueMap;
    typedef std::vector < Eigen::Triplet<double> >     TripletList;

public:
    //! constructor
    ForwardPerturbation ( Solver * const solver, const CovarForcing * const covarForcing,
                          CovarMaterial * covarMaterial, const double regFactor = 1.e-4 )
            : solver_( solver ), covarForcing_( covarForcing ),
              covarMaterial_( covarMaterial ), regFactor_( regFactor )
    { }

    //! set preamble for hyperparameter names
    void setHypIdPre ( const std::string & preamble )
    {
        hypIdPre_ = preamble;
    }

    void setProjectionMatrix ( const Eigen::MatrixXd & pMat )
    {
        pMat_ = pMat;
    }

    //! set the vector of TripletList containing the derivatives of the stiffness matrix
    //! with respect to the material parameters
    void setStiffnessDerivatives( const std::vector< TripletList > & stiffnessDerivs )
    {
        stiffnessDerivs_ = stiffnessDerivs;
    }

    //! give the mean \f$ \Phi A^{-1}(\overline w) \overline f \f$
    void giveProjectedMeanU ( Eigen::VectorXd & projMean ) const
    {
        FTL_VERIFY( feMean_.size( ) );
        projMean = pMat_ * feMean_;

        return;
    }

    //! give projected covariance \f$ \vec P \vec C_u \vec P^T \f$
    void giveProjectedCovarU ( Eigen::MatrixXd & projCovmat ) const
    {
        FTL_VERIFY( feCovMat_.size( ) );
        projCovmat = feCovMat_;

        return;
    }

    //! set update solution or use the previously computed solution
    void updateSolution ( const HpsValueMap &params )
    {
        //! get numDofs and populate dofIndices
        const unsigned numDofs = pMat_.cols( );
        std::vector<unsigned> dofIndices( numDofs );
        for ( unsigned i = 0; i < numDofs; ++ i ) {
            dofIndices[ i ] = i;
        }

        //! get solution at all nodes
        feMean_.resize( numDofs );
        solver_->giveSolution( dofIndices, feMean_ );

        //! solve auxiliary problem \f$ A^T * R^T = P^T \f$
        FTL_VERIFY( pMat_.size( ) );
        const unsigned numObservations = pMat_.rows( );

        Eigen::MatrixXd pMatId = Eigen::MatrixXd::Identity( numDofs, numDofs );
        Eigen::MatrixXd aInv( numDofs, numDofs );
        Eigen::MatrixXd rMat( numObservations, numDofs );

        //! one solve per column of \$ P \f$
        bool reuseFactorisation = false;
        for ( unsigned i = 0; i < numObservations; ++ i ) {
            Eigen::VectorXd tmp, tmpId;
            solver_->solveSystemTranspose( pMat_.row( i ), tmp,
                                           reuseFactorisation );
            reuseFactorisation = true;
            solver_->solveSystemTranspose( pMatId.row( i ), tmpId,
                                           reuseFactorisation );
            rMat.row( i )   = tmp;
            aInv.row( i ) = tmpId;
        }

        //! compute covariance of FE solution and store in feCovMat_
        //! first the second term of equation (27) in the statFem paper is computed
        //! \f$ A^{-1} C_f A^{-1} + ( u \otimes u ) \f$
        Eigen::MatrixXd forceCovMat;
        covarForcing_ -> giveProjectedCovar( params, aInv, forceCovMat );
        forceCovMat += ( feMean_ * feMean_.transpose( ) );

        //! get material covariance
        Eigen::MatrixXd materialCovMat;
        covarMaterial_ -> giveMaterialCovariance( params, materialCovMat );

        FTL_VERIFY( stiffnessDerivs_.size( ) );
        unsigned nDeriv = stiffnessDerivs_.size( );

        Eigen::MatrixXd covarUmaterial;
        covarUmaterial.resize( numObservations, numObservations );
        covarUmaterial.setZero( );
        for ( unsigned e = 0; e < nDeriv; ++e ) {
            //! create sparse matrix from triplets
            Eigen::SparseMatrix<double> derivE ( numDofs, numDofs );
            TripletList stiffnessDerivWRTe = stiffnessDerivs_[ e ];
            derivE.setFromTriplets( stiffnessDerivWRTe.begin( ), stiffnessDerivWRTe.end( ) );

            for ( unsigned d = 0; d < nDeriv; ++d ) {
                //! create sparse matrix from triplets
                Eigen::SparseMatrix<double> derivD ( numDofs, numDofs );
                TripletList stiffnessDerivWRTd = stiffnessDerivs_[ d ];
                derivD.setFromTriplets( stiffnessDerivWRTd.begin( ), stiffnessDerivWRTd.end( ) );

                covarUmaterial += materialCovMat( e, d ) * rMat * derivE * forceCovMat *
                        derivD.transpose( ) * rMat.transpose( );

            }
        }

        //! the projected covariance using perturbation method
        feCovMat_.resize( numObservations, numObservations );
        feCovMat_.setZero( );

        //! add the first term in equation (27) in the statFEM paper
        covarForcing_ -> giveProjectedCovar( params, rMat, feCovMat_ );
        feCovMat_ += covarUmaterial;

        //! add regularisation term to improve conditioning of the covariance matrix
        auto it = params.find( "material.sigma" );
        FTL_VERIFY( it != params.end( ) );
        const double sigma = it -> second;
        double addDiag = regFactor_ * sigma;

        feCovMat_ += ( addDiag * addDiag * pMat_ * pMat_.transpose( ) );

        return ;
    }

private:
    Solver                  * const solver_;
    std::string               hypIdPre_;

    const CovarForcing      * const covarForcing_;
    const CovarMaterial     * const covarMaterial_;
    Eigen::VectorXd           feMean_;

    Eigen::MatrixXd           pMat_;
    Eigen::MatrixXd           feCovMat_;

    std::vector<TripletList > stiffnessDerivs_;

    //! default regularisation factor as fraction of \f$ \sigma_{\kappa} \f$
    const double              regFactor_;

};

#endif
