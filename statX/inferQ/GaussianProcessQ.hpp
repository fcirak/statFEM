// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//
//                 Computational Structural Mechanics Lab
//                         University of Cambridge
//
// This software is copyrighted and all rights are retained by the University of
// Cambridge. All use, disclosure, and/or reproduction of any part not expressly
// authorised by the University of Cambridge is prohibited. (C) 2010 --
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

//! @file GaussianProcessQ.hpp

#ifndef statX_inferQ_gaussianprocess_h
#define statX_inferQ_gaussianprocess_h

//------------------------------------------------------------------------------
// system includes
#include <vector>
#include <string>
#include <map>
#include <cmath>
// Eigen includes
#include <Eigen/Core>
#include <Eigen/Sparse>
// corlib includes
#include <corlib/verify.hpp>

//------------------------------------------------------------------------------
namespace statX {
    namespace inferQ {
        class GaussianProcessQ;
    }
}

//------------------------------------------------------------------------------
/** \brief Gaussian process regression using sparse precision matrix
 *  \details
 *
 *  - By establishing that \f$ u = P_r v \f$, the statistical model reads
 *    \f[
 *        y = P u + e = ( P P_r ) v + e = P_v v    \quad \text{with} \quad
 *        v \sim \mathcal{N}( 0, Q_v^{-1}     )     \,, \quad
 *        e \sim \mathcal{N}( 0, \sigma_e^2 I )
 *    \f]
 *    where \f$ Q_v = P_l^\trans Q P_l \f$ is the sparse precision matrix.
 *    The strategy is to first infer \f$ v \f$. The actual variable is inferred
 *    as \f$ u = P_r v \f$.
 *
 *  - For the integer exponent case \f$\beta = 1, 2, \dotsc \f$, the right and
 *    left fractional matrices should be set as
 *    \f[
 *        P_r := I \quad \text{and} \quad P_l := I,
 *    \f]
 *    so that the statistical model simplifies into
 *    \f[
 *        y = P u + e                              \quad \text{with} \quad
 *        u \sim \mathcal{N}( 0, Q^{-1}       )     \,, \quad
 *        e \sim \mathcal{N}( 0, \sigma_e^2 I )
 *    \f]
 *
 *  - The observation matrix \f$ Y = ( y_1, \dotsc, y_{n_o} ) \in
 *    \mathbb{R}^{n_y \times n_o} \f$ represents \f$ n_o \f$ readings taken
 *    at each of the \f$ n_y \f$ data points.
 *
 *  - Provides evaluation of log marginal likelihood \f$ \log p(Y) \f$ and
 *    posterior \f$ p(u \vert Y) \f$.
 *
 *  - By default, Eigen sparse Cholesky solver applies permutation onto
 *    the matrix to reduce fill-in for \f$ LL^T \f$ decomposition.
 *    Mathematically, the permutation does not affect the computation of log
 *    determinant of the matrix because the determinant of a permutation
 *    matrix is either 1 or -1.
 */
class statX::inferQ::GaussianProcessQ
{
public:

    //! @name Convenience typedefs, publicly available
    ///@{
    typedef std::map<std::string, double>    HpsValueMap;
    ///@}

private:

    //! @name Convenience typedefs
    ///@{
    typedef Eigen::SimplicialLLT
            <Eigen::SparseMatrix<double>>    SparseLLT_;
    ///@}

public:

    //! @name Constructor and destructor
    ///@{

    //! Construct with observation values and projection matrix
    GaussianProcessQ ( const Eigen::MatrixXd & obsValues,
                       const Eigen::SparseMatrix<double> & pMat ) :
        obsValues_( obsValues ),
        pMat_( pMat )
    {
        // size check
        FTL_VERIFY_DESCRIPTIVE( pMat_.rows( ) == obsValues_.rows( ),
                                "Number of data points not matched.\n" );

        // observation error \sigma_e
        hypIds_.push_back( "eps" );
    }

    //! Set preamble for hyperparameter names
    void setHypIdPre ( const std::string & preamble );

    ///@}

public:

    //! Give the log marginal likelihood \f$ \log p( Y ) \f$
    double giveLogLikelihood ( const HpsValueMap & hParams,
                               const Eigen::SparseMatrix<double> & qMat,
                               const Eigen::SparseMatrix<double> & prMat );

    //! Give the mean vector of the posterior \f$  p( u \vert Y ) \f$
    void predictMean ( const HpsValueMap & hParams,
                       const Eigen::SparseMatrix<double> & qMat,
                       const Eigen::SparseMatrix<double> & prMat,
                       Eigen::VectorXd & uVec );

    //! Give the covariance matrix of the posterior \f$ p( u \vert Y ) \f$
    void predictCovariance ( const HpsValueMap & hParams,
                             const Eigen::SparseMatrix<double> & qMat,
                             const Eigen::SparseMatrix<double> & prMat,
                             Eigen::MatrixXd & uCovMat );

    //! Give the variance vector of the posterior \f$ p( u \vert Y ) \f$
    void predictVariance ( const HpsValueMap & hParams,
                           const Eigen::SparseMatrix<double> & qMat,
                           const Eigen::SparseMatrix<double> & prMat,
                           Eigen::VectorXd & uVarVec );

private:

    //! Helper function to give observation error \f$ \sigma_e \f$
    double extractEps_ ( const HpsValueMap & hParams ) const;

private:
    std::string                  hypIdPre_;  //!< hyperparameters preamble
    std::vector<std::string>     hypIds_;    //!< hyperparameters name

    Eigen::MatrixXd              obsValues_; //!< observation values
    Eigen::SparseMatrix<double>  pMat_;      //!< projection matrix
    SparseLLT_                   solver_;    //!< direct sparse Cholesky solver

};

//------------------------------------------------------------------------------
#include "GaussianProcessQ.ipp"
//------------------------------------------------------------------------------

#endif
