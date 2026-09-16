// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//
//                 Computational Structural Mechanics Lab
//                         University of Cambridge
//
// This software is copyrighted and all rights are retained by the University of
// Cambridge. All use, disclosure, and/or reproduction of any part not expressly
// authorised by the University of Cambridge is prohibited. (C) 2010 --
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

//! @file PrecisionMatrixRational.ipp

#include <corlib/verify.hpp>

//------------------------------------------------------------------------------
void statX::inferQ::PrecisionMatrixRational::
setHypIdPre( const std::string & preamble )
{
    hypIdPre_ = preamble;
    return;
}

//------------------------------------------------------------------------------
void statX::inferQ::PrecisionMatrixRational::
givePrecisionMatrix( const HpsValueMap & hps,
                     const Eigen::SparseMatrix<double> & maternMat,
                     const Eigen::VectorXd & massVec,
                     Eigen::SparseMatrix<double> & qMat )
{
    Eigen::SparseMatrix<double> qIntMat;
    this->givePrecisionMatrixInteger_( hps, maternMat, massVec, qIntMat );

    Eigen::SparseMatrix<double> plMat;
    this->giveLeftFracMatrix( hps, maternMat, massVec, plMat );

    qMat = plMat.transpose( ) * qIntMat * plMat;

    return;
}

//------------------------------------------------------------------------------
/**
 * The non-symmetric matrix \f$ P_l \f$ is defined as
 * \f[
 *       P_l = b_m \prod_{ j = 0 }^{ m - 1 } \left( I - d_j M^{-1} A \right)
 * \f]
 */
void statX::inferQ::PrecisionMatrixRational::
giveLeftFracMatrix( const HpsValueMap & hps,
                    const Eigen::SparseMatrix<double> & maternMat,
                    const Eigen::VectorXd & massVec,
                    Eigen::SparseMatrix<double> & plMat )
{
    // get hyperparameter values
    HypArray hpsArray;
    this->giveHyperparameters_( hps, hpsArray );
    const double sigma = hpsArray[ 0 ];
    const double kappa = hpsArray[ 1 ];
    const double nu    = hpsArray[ 2 ];

    // check validity of hyperparameters
    FTL_VERIFY( ( sigma > 0. ) and ( kappa > 0. ) and ( nu > 0. ) );

    // get the fractional part of \beta
    const double beta = ( nu + ( double ) dim_ / 2. ) / 2.;
    const double alpha = std::max( 1., std::floor( beta ) );
    const double exponent = beta - alpha;

    // set to identity matrix if the exponent is zero
    if ( std::fabs( exponent ) < 1.e-12 )
    {
        plMat.resize( massVec.size( ), massVec.size( ) );
        plMat.setIdentity( );

        // exit this function
        return;
    }

    // rational approximation of x^exponent
    std::vector<double> rootsPR; // not used
    std::vector<double> rootsPL;
    this->brasil_( exponent, rootsPR, rootsPL );

    // inverse of lumped mass matrix
    Eigen::DiagonalMatrix<double, Eigen::Dynamic> mMatInv( massVec.size( ) );
    mMatInv.diagonal( ) = massVec.cwiseInverse( );

    // identity matrix
    Eigen::SparseMatrix<double> idenMat( massVec.size( ), massVec.size( ) );
    idenMat.setIdentity( );

    // compute the matrix P_l
    plMat.resize( massVec.size( ), massVec.size( ) );
    plMat.setIdentity( );
    for ( unsigned j = 0; j < degree_; ++ j )
        plMat = plMat * ( idenMat - rootsPL[ j + 1 ] * mMatInv * maternMat );
    plMat = rootsPL[ 0 ] * plMat;

    return;
}

//------------------------------------------------------------------------------
/**
 * The non-symmetric matrix \f$ P_r \f$ is defined as
 * \f[
 *       P_r = a_m \prod_{ i = 0 }^{ m - 1 } \left( I - c_i M^{-1} A \right)
 * \f]
 */
void statX::inferQ::PrecisionMatrixRational::
giveRightFracMatrix( const HpsValueMap & hps,
                     const Eigen::SparseMatrix<double> & maternMat,
                     const Eigen::VectorXd & massVec,
                     Eigen::SparseMatrix<double> & prMat )
{
    // get hyperparameter values
    HypArray hpsArray;
    this->giveHyperparameters_( hps, hpsArray );
    const double sigma = hpsArray[ 0 ];
    const double kappa = hpsArray[ 1 ];
    const double nu    = hpsArray[ 2 ];

    // check validity of hyperparameters
    FTL_VERIFY( ( sigma > 0. ) and ( kappa > 0. ) and ( nu > 0. ) );

    // get the fractional part of \beta
    const double beta = ( nu + ( double ) dim_ / 2. ) / 2.;
    const double alpha = std::max( 1., std::floor( beta ) );
    const double exponent = beta - alpha;

    // set to identity matrix if the exponent is zero
    if ( std::fabs( exponent ) < 1.e-12 )
    {
        prMat.resize( massVec.size( ), massVec.size( ) );
        prMat.setIdentity( );

        // exit this function
        return;
    }

    // rational approximation of x^exponent
    std::vector<double> rootsPR;
    std::vector<double> rootsPL; // not used
    this->brasil_( exponent, rootsPR, rootsPL );

    // inverse of lumped mass matrix
    Eigen::DiagonalMatrix<double, Eigen::Dynamic> mMatInv( massVec.size( ) );
    mMatInv.diagonal( ) = massVec.cwiseInverse( );

    // identity matrix
    Eigen::SparseMatrix<double> idenMat( massVec.size( ), massVec.size( ) );
    idenMat.setIdentity( );

    // compute the matrix P_r
    prMat.resize( massVec.size( ), massVec.size( ) );
    prMat.setIdentity( );
    for ( unsigned i = 0; i < degree_; ++ i )
        prMat = prMat * ( idenMat - rootsPR[ i + 1 ] * mMatInv * maternMat );
    prMat = rootsPR[ 0 ] * prMat;

    return;
}

//------------------------------------------------------------------------------
/**
 * Given \f$ \mathcal{L}^\alpha \f$ with \f$ \alpha = \max( 1,\lfloor \beta
 * \rfloor ) \in \mathbb{N} \f$, the precision matrix \f$ Q \equiv Q_j \f$ is
 * computed recursively as
 * \f{align}{
 *         Q_1 = \kappa^{4 \beta} \tau^2 A M^{-1} A  \quad &\text{for } j = 1
 *         \\
 *         Q_j = A M^{-1} Q_{j-1} M^{-1} A
 *         \quad &\text{for } j = 2,3,\dotsc, \alpha \in \mathbb{N}
 * \f}
 * where \f$ A =  M + \kappa^{-2} G \f$ is the scaled Matern stiffness matrix,
 * \f$ G \f$ is the Laplacian matrix and \f$ M \f$ is the lumped mass matrix.
 */
void statX::inferQ::PrecisionMatrixRational::
givePrecisionMatrixInteger_( const HpsValueMap & hps,
                             const Eigen::SparseMatrix<double> & maternMat,
                             const Eigen::VectorXd & massVec,
                             Eigen::SparseMatrix<double> & qIntMat )
{
    // get hyperparameter values
    HypArray hpsArray;
    this->giveHyperparameters_( hps, hpsArray );
    const double sigma = hpsArray[ 0 ];
    const double kappa = hpsArray[ 1 ];
    const double nu    = hpsArray[ 2 ];

    // check validity of hyperparameters
    FTL_VERIFY( ( sigma > 0. ) and ( kappa > 0. ) and ( nu > 0. ) );

    // get the integer part of \beta
    const double beta = ( nu + ( double ) dim_ / 2. ) / 2.;
    const double alpha = std::max( 1., std::floor( beta ) );

    // inverse of lumped mass matrix
    Eigen::DiagonalMatrix<double, Eigen::Dynamic> mMatInv( massVec.size( ) );
    mMatInv.diagonal( ) = massVec.cwiseInverse( );

    // compute the precision matrix recursively
    this->applyRecursion_( ( unsigned int ) alpha, maternMat, mMatInv, qIntMat );

    // the scaling \kappa ^ { 4 \beta } is applied now instead of next layers
    const double tau = this->giveTau_( hps );
    qIntMat = tau * tau * std::pow( kappa, 4. * beta ) * qIntMat;

    return;
}

//------------------------------------------------------------------------------
void statX::inferQ::PrecisionMatrixRational::
giveHyperparameters_( const HpsValueMap & hps,
                      HypArray & hpsArray ) const
{
    // find \sigma value
    auto it = hps.find( hypIdPre_ + hypIds_[ 0 ] );
    FTL_VERIFY( it != hps.end( ) );
    const double sigma = it->second;

    // find \nu value
    it = hps.find( hypIdPre_ + hypIds_[ 2 ] );
    FTL_VERIFY( it != hps.end( ) );
    const double nu = it->second;

    // find \kappa or \ell value
    it = hps.find( hypIdPre_ + hypIds_[ 1 ] );
    double kappa;
    if ( it != hps.end( ) )
        kappa = it->second;
    else {
        auto itr = hps.find( hypIdPre_ + "length" );
        FTL_VERIFY( itr != hps.end( ) );
        const double length = itr->second;
        kappa = std::sqrt( 2. * nu ) / length;
    }

    // update hyperparameter values
    hpsArray[ 0 ] = sigma;
    hpsArray[ 1 ] = kappa;
    hpsArray[ 2 ] = nu;

    return;
}

//------------------------------------------------------------------------------
/**
 * - Given \f$ \sigma \f$, \f$ \nu \f$ and \f$ \kappa \f$, the scaling
 *   constant \f$ \tau \f$ is defined as
 *   \f[
 *          \tau^2 = \frac{ \Gamma( \nu ) }{ \sigma^2
 *                   \Gamma( \nu + d / 2 ) ( 4 \pi )^{d/2} \kappa ^ {2 \nu} }
 *   \f]
 *   where \f$ d \f$ is the spatial dimension.
 * - Calling the gamma function \f$ \Gamma \f$ requires c++17.
 *
 */
double statX::inferQ::PrecisionMatrixRational::
giveTau_( const HpsValueMap & hps ) const
{
    // get hyperparameter values
    HypArray hpsArray;
    this->giveHyperparameters_( hps, hpsArray );
    const double sigma = hpsArray[ 0 ];
    const double kappa = hpsArray[ 1 ];
    const double nu    = hpsArray[ 2 ];

    // return \tau
    double tauSquared = std::tgamma( nu )
            / ( ( sigma * sigma ) * std::tgamma( nu + ( double ) dim_ / 2. )
                    * std::pow( 4. * M_PI, ( double ) dim_ / 2. )
                    * std::pow( kappa, 2. * nu ) );

    return std::sqrt( tauSquared );
}

//------------------------------------------------------------------------------
void statX::inferQ::PrecisionMatrixRational::
applyRecursion_ ( unsigned j,
                  const Eigen::SparseMatrix<double> & maternMat,
                  const Eigen::DiagonalMatrix<double, Eigen::Dynamic> & mMatInv,
                  Eigen::SparseMatrix<double> & qIntMat ) const
{
    // \floor{ \beta } = 1
    if ( j == 1 )
    {
        qIntMat = maternMat * mMatInv * maternMat;
        return;
    }

    // recursive cases
    applyRecursion_( j - 1, maternMat, mMatInv, qIntMat );

    // update
    qIntMat = maternMat * mMatInv * qIntMat * mMatInv * maternMat;

    return;
}

//------------------------------------------------------------------------------
/**
 * - Performs rational approximation of the power function \f$ g(x) =
 *  x^{\beta - \alpha} \f$ with the exponent \f$ \beta -\alpha \in ( -1, 0 )
 *  \cup ( 0, +1 ) \subset \mathbb R \f$ using the brasil algorithm.
 *
 * - The rational approximation reads
 *  \f[
 *     g(x) \approx \frac{ \sum_{ i = 0 }^{ m } a_i x^i }
 *                       { \sum_{j=0}^m b_j  x^j } =
 *                  \frac{ a_m \prod_{ i = 0 }^{ m - 1 } ( x- c_i ) }
 *                       { b_m \prod_{ j = 0 }^{ m - 1 } ( x- d_j ) }
 *  \f]
 *
 * - The \f$ m \f$ numerator factors \f$ a_m \f$ and \f$ c_i \f$ lead to the
 *   matrix \f$ P_r \f$
 *
 * - The \f$ m \f$ denominator factors \f$ b_m \f$ and \f$ d_j \f$ lead to the
 *   the matrix \f$ P_l \f$.
 *
 * - The storage sequence is \f$ a_m \f$ (resp. \f$ b_m \f$) followed by the
 *   \f$ m - 1\f$ roots \f$ c_i \f$ (resp. \f$ d_j \f$).
 */
void statX::inferQ::PrecisionMatrixRational::
brasil_( const double & exponent,
         std::vector<double> & rootsPR,
         std::vector<double> & rootsPL  )
{
    // exponent of power function must not be zero
    FTL_VERIFY( std::fabs( exponent ) > 0. );

    // convert double to std::string with maximum precision
    std::stringstream stream;
    stream << std::fixed
           << std::setprecision( std::numeric_limits<double>::digits10 + 1 )
           << exponent;

    // call python script from command line
    const std::string brasil =
            "python3 $STATFEMROOT/tools/baryrat/baryrat-1.4.0/roots.py";
    const std::string argv1 = std::to_string( degree_ );
    const std::string argv2 = stream.str( );
    const std::string line = brasil + " " + argv1 + " " + argv2;
    const int status = std::system( line.c_str( ) );
    FTL_VERIFY_DESCRIPTIVE(
            status == 0,
            "Cannot call python3 $STATFEMROOT/tools/baryrat/baryrat-1.4.0/roots.py." );

    // read approximation result
    this->readRoots_( "top.out", rootsPR );
    this->readRoots_( "bot.out", rootsPL );

    return;
}

//------------------------------------------------------------------------------
void statX::inferQ::PrecisionMatrixRational::
readRoots_ ( const std::string & fileName, std::vector<double> & roots )
{
    std::ifstream fileStream( fileName );
    FTL_VERIFY( fileStream.is_open( ) );
    roots.reserve( degree_ + 1 );
    for( unsigned i = 0; i < degree_ + 1; ++ i )
    {
        double root;
        fileStream >> root;
        roots.push_back( root );
    }

    fileStream.close( );

    return;
}
