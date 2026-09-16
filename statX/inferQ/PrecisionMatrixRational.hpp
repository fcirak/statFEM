// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//
//                 Computational Structural Mechanics Lab
//                         University of Cambridge
//
// This software is copyrighted and all rights are retained by the University of
// Cambridge. All use, disclosure, and/or reproduction of any part not expressly
// authorised by the University of Cambridge is prohibited. (C) 2010 --
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

//! @file PrecisionMatrixRational.hpp

#ifndef statX_inferQ_precisionmatrixrational_h
#define statX_inferQ_precisionmatrixrational_h

//------------------------------------------------------------------------------
// system includes
#include <array>
#include <cmath>
#include <map>
#include <vector>
#include <string>
#include <sstream>
#include <fstream>
#include <cstdlib>
#include <iomanip>
#include <limits>
// Eigen includes
#include <Eigen/Core>
#include <Eigen/Sparse>

//------------------------------------------------------------------------------
namespace statX {
    namespace inferQ {
        class PrecisionMatrixRational;
    }
}


//------------------------------------------------------------------------------
/** \brief Precision matrix for real-valued exponent of linear operator
 *  \details
 *
 *  - Applicable to all smoothness \f$ \nu > 0 \f$, i.e.
 *    \f$ \beta > \frac{d}{4} \f$ where \f$ d \f$ is the spatial dimension.
 *
 *  - Provides query of the precision matrix given the Matern parameters
 *    \f$ \nu \f$, \f$ \ell \f$ (or \f$ \kappa \f$) and \f$ \sigma \f$.
 *
 *  - Computes the final precision matrix in three layers:
 *    -# \f$ Q \f$: sparse and is equal to Lindgren precision matrix
 *                  for \f$ \beta \in \mathbb{N} \f$
 *
 *    -# \f$ P_l^\trans Q P \f$: sparse
 *
 *    -# \f$ P_r^{-\trans} P_l^\trans Q P_l P_r^{-1} \f$: dense
 *
 *  - Follows the paper by Bolin and Kirchner (2020) with some exceptions:
 *
 *    -# We use the same polynomial degree for both numerator and denominator
 *       of the rational interpolant.
 *
 *    -# We generalise Lindgren et al. (2011) by computing the final precision
 *       matrix in three layers stated above.
 *
 *    -# We use the rational interpolation algorithm by Hofreither (2020). See
 *       tools/baryrat for installation guide.
 *
 *  - Bolin, David, and Kristin Kirchner. "The rational SPDE approach for
 *    Gaussian random fields with general smoothness." Journal of Computational
 *    and Graphical Statistics 29.2 (2020): 274-285.
 *
 *  - Hofreither, Clemens. "An algorithm for best rational approximation based
 *    on barycentric rational interpolation." Numerical Algorithms (2021): 1-24.
 *
 */
class statX::inferQ::PrecisionMatrixRational
{
public:

    //! @name Basic attributes
    ///@{
    static constexpr unsigned  numHps = 3;    //!< number of hyperparameters
    ///@}

    //! @name Convenience typedefs, publicly available
    ///@{
    typedef std::map<std::string, double>                    HpsValueMap;
    typedef std::array<std::string, numHps>                  HypIdsArray;
    typedef std::array<double, numHps>                       HypArray;
    ///@}

public:

    //! @name Constructor and destructor
    ///@{

    //! Construct with dimension and rational approximation polynomial degree
    PrecisionMatrixRational ( unsigned dim,
                              unsigned degree = 4 ) :
        dim_( dim ),
        degree_( degree )
    {
        // hyperparameter names
        hypIds_[ 0 ] = "sigma";
        hypIds_[ 1 ] = "kappa"; // the alternative for "kappa" is "length"
        hypIds_[ 2 ] = "nu";

        // check if baryrat is installed in user's python environment
        const int status = std::system( "python3 -c 'import baryrat'" );
        FTL_VERIFY_DESCRIPTIVE( status == 0,
                                "Cannot call python3 -c 'import baryrat'." );
    }
    ///@}

public:

    //! Set preamble for hyperparameter names
    void setHypIdPre ( const std::string & preamble );

    //! Give precision matrix \f$ P_l^\trans Q P_l \f$
    void givePrecisionMatrix ( const HpsValueMap & hps,
                               const Eigen::SparseMatrix<double> & maternMat,
                               const Eigen::VectorXd & massVec,
                               Eigen::SparseMatrix<double> & qMat );

    //! Give the non-symmetric matrix \f$ P_l \f$
    void giveLeftFracMatrix  ( const HpsValueMap & hps,
                               const Eigen::SparseMatrix<double> & maternMat,
                               const Eigen::VectorXd & massVec,
                               Eigen::SparseMatrix<double> & plMat );

    //! Give the non-symmetric matrix \f$ P_r \f$
    void giveRightFracMatrix ( const HpsValueMap & hps,
                               const Eigen::SparseMatrix<double> & maternMat,
                               const Eigen::VectorXd & massVec,
                               Eigen::SparseMatrix<double> & prMat );

private:

    //! Give precision matrix \f$ Q \f$
    void givePrecisionMatrixInteger_ (
            const HpsValueMap & hps,
            const Eigen::SparseMatrix<double> & maternMat,
            const Eigen::VectorXd & massVec,
            Eigen::SparseMatrix<double> & qIntMat );

    //! Give the array of Matern hyperparameter values
    void giveHyperparameters_ ( const HpsValueMap & hps,
                                HypArray & hpsArray ) const;

    //! Return the scaling constant \f$ \tau \f$
    double giveTau_ ( const HpsValueMap & hps ) const;

    //! Compute the precision matrix \f$ Q \f$ recursively
    void applyRecursion_ ( unsigned j,
                           const Eigen::SparseMatrix<double> & maternMat,
                           const Eigen::DiagonalMatrix
                               <double, Eigen::Dynamic>      & mMatInv,
                           Eigen::SparseMatrix<double>       & qIntMat ) const;

    //! Perform rational approximation using the brasil algorithm
    void brasil_ ( const double & exponent,
                   std::vector<double> & rootsPR,
                   std::vector<double> & rootsPL );

    //! Helper function to read a vector of roots from a file stream
    void readRoots_ ( const std::string & fileName, std::vector<double> & out );

private:

    HypIdsArray         hypIds_;        //!< array of hyperparameter names
    std::string         hypIdPre_;      //!< preamble of hyperparameter names

    const unsigned      dim_;           //!< spatial dimension
    const unsigned      degree_;        //!< polynomial degree

};

//------------------------------------------------------------------------------
#include "PrecisionMatrixRational.ipp"
//------------------------------------------------------------------------------

#endif
