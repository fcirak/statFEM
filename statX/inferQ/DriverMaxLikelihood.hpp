// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//
//                 Computational Structural Mechanics Lab
//                         University of Cambridge
//
// This software is copyrighted and all rights are retained by the University of
// Cambridge. All use, disclosure, and/or reproduction of any part not expressly
// authorised by the University of Cambridge is prohibited. (C) 2010 --
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

//! @file DriverMaxLikelihood.hpp

#ifndef statX_inferQ_drivermaxlikelihood_h
#define statX_inferQ_drivermaxlikelihood_h

// system headers
#include <cmath>
#include <vector>
#include <string>

// boost headers
#include <boost/function.hpp>

// corlib headers
#include <corlib/Integrator.hpp>
#include <corlib/ComputeElementMatrix.hpp>
#include <corlib/MatrixAssembler.hpp>
#include <corlib/VectorAssembler.hpp>

// eigen headers
#include <Eigen/Core>
#include <Eigen/Sparse>

//------------------------------------------------------------------------------
namespace statX {
    namespace inferQ {
        template <typename MESH, typename QUAD, typename MATERNASSEMBLER,
                  typename PRECISION, typename GP>
        class DriverMaxLikelihood;
    }
}

//------------------------------------------------------------------------------
/**
 * \brief Maximum likelihood estimate driver for gaussian process regression
 *        using sparse precision matrix.
 * \details Provides convenience functions for interfacing gaussian process,
 *          precision matrix and optimisation.
 */
template <typename MESH, typename QUAD, typename MATERNASSEMBLER,
          typename PRECISION, typename GP>
class statX::inferQ::DriverMaxLikelihood
{
public:
    typedef MESH                               Mesh;
    typedef QUAD                               Quad;
    typedef MATERNASSEMBLER                    MaternAssembler;
    typedef PRECISION                          Precision;
    typedef GP                                 GaussianProcessQ;

public:
    typedef typename Mesh::Element             Element;
    typedef typename Element::Node             Node;
    typedef typename Element::MatDimDim        MatDimDim;
    typedef typename Node::VecDim              VecDim;
    typedef typename Precision::HpsValueMap    HpsValueMap;

public:

    //! Construct with interfacing objects and pointers
    DriverMaxLikelihood( Mesh *             mesh,
                         Quad               quadrature,
                         MaternAssembler *  assembler,
                         Precision *        precision,
                         GaussianProcessQ * gp ) :
        mesh_      ( mesh       ),
        quadrature_( quadrature ),
        assembler_ ( assembler  ),
        precision_ ( precision  ),
        gp_        ( gp         )
    {
        // set preamble
        precision->setHypIdPre( "gaussianProcess." );
        gp_->setHypIdPre( "gaussianProcess." );

        // lump mass needs only to compute once
        this->computeLumpedMassVector_( );
    }

public:

    //! Set all hyperparameters
    void setAllHps( const HpsValueMap & hps ) { hps_ = hps; }

    //! Give all hyperparameters
    void giveAllHps( HpsValueMap & hps ) const { hps = hps_; }

    //! Add the name of a hyperparameter involved in optimisation
    void addOptimHps( const std::string & name ) { names_.push_back( name ); }

    //! Update hyperparameters involved in optimisation
    void updateHps( const std::vector<double> & params );

    //! Convenience function to compute log marginal likelihood
    double computeLogLikelihood( );

    //! Convenience function to compute mean values at mesh points
    void predictMean( Eigen::VectorXd & mean );

    //! Convenience function to compute covariance of all mesh points
    void predictCovariance( Eigen::MatrixXd & covariance );

private:

    //! Convenience function to compute global Matern stiffness matrix
    void computeMaternMatrix_ ( const double & kappa,
                                const MatDimDim & diffusivity =
                                       MatDimDim::Identity( ),
                                const bool & isScaled = true );

    //! Convenience function to compute global lumped mass vector
    void computeLumpedMassVector_( );

    //! Convenience function to compute precision matrix
    void computePrecisionMatrix_( Eigen::SparseMatrix<double> & qMat,
                                  Eigen::SparseMatrix<double> & prMat );

    //! Convenience function to extract kappa
    double extractKappa_ ( ) const;

private:

    Mesh *                mesh_;
    MaternAssembler *     assembler_;
    Quad                  quadrature_;
    Precision *           precision_;
    GaussianProcessQ *    gp_;

private:

    std::vector<std::string>    names_; //< vector of optimisation parameter names
    HpsValueMap                 hps_;   //< map of hyperparameters
};

//------------------------------------------------------------------------------
template <typename MESH, typename QUAD, typename MATERNASSEMBLER,
          typename PRECISION, typename GP>
void statX::inferQ::DriverMaxLikelihood<MESH, QUAD, MATERNASSEMBLER, PRECISION, GP>::
updateHps( const std::vector<double> & params )
{
    for( unsigned k = 0; k < params.size( ); ++k )
        hps_[ names_[ k ] ] = params[ k ];

    return;
}

//------------------------------------------------------------------------------
template <typename MESH, typename QUAD, typename MATERNASSEMBLER,
          typename PRECISION, typename GP>
double statX::inferQ::DriverMaxLikelihood<MESH, QUAD, MATERNASSEMBLER, PRECISION, GP>::
computeLogLikelihood( )
{
    Eigen::SparseMatrix<double> qMat;
    Eigen::SparseMatrix<double> prMat;
    this->computePrecisionMatrix_( qMat, prMat );

    return gp_->giveLogLikelihood( hps_, qMat, prMat );
}

//------------------------------------------------------------------------------
template <typename MESH, typename QUAD, typename MATERNASSEMBLER,
          typename PRECISION, typename GP>
void statX::inferQ::DriverMaxLikelihood<MESH, QUAD, MATERNASSEMBLER, PRECISION, GP>::
predictMean( Eigen::VectorXd & mean )
{
    Eigen::SparseMatrix<double> qMat;
    Eigen::SparseMatrix<double> prMat;
    this->computePrecisionMatrix_( qMat, prMat );
    gp_->predictMean( hps_, qMat, prMat, mean );

    return;
}

//------------------------------------------------------------------------------
template <typename MESH, typename QUAD, typename MATERNASSEMBLER,
          typename PRECISION, typename GP>
void statX::inferQ::DriverMaxLikelihood<MESH, QUAD, MATERNASSEMBLER, PRECISION, GP>::
predictCovariance( Eigen::MatrixXd & covariance )
{
    Eigen::SparseMatrix<double> qMat;
    Eigen::SparseMatrix<double> prMat;
    this->computePrecisionMatrix_( qMat, prMat );
    gp_->predictCovariance( hps_, qMat, prMat, covariance );

    return;
}

//------------------------------------------------------------------------------
template <typename MESH, typename QUAD, typename MATERNASSEMBLER,
          typename PRECISION, typename GP>
void statX::inferQ::DriverMaxLikelihood<MESH, QUAD, MATERNASSEMBLER, PRECISION, GP>::
computeMaternMatrix_( const double & kappa,
                      const MatDimDim & diffusivity,
                      const bool & isScaled )
{
    // clear global triplet storage
    assembler_->clearMaternMatrix( );

    // set kappa and diffusivity
    mesh_->iterateOverElements( boost::bind( &Element::setKappa, _1, kappa ) );
    mesh_->iterateOverElements(
            boost::bind( &Element::setDiffusivity, _1, diffusivity ) );

    // scaling factor
    double factor = 1.;
    if( isScaled )
        factor = 1. / ( kappa * kappa );

    // assemble global Matern stiffness matrix
    corlib::MatrixComputeAndAssembleFun<const Element, Quad, MaternAssembler>
        maternComputer( &Element::maternIntegrand,
                        &Element::getDofIndices,
                        &Element::getDofIndices,
                        quadrature_, assembler_, factor );
    mesh_->iterateOverElements( maternComputer );

    return;
}

//------------------------------------------------------------------------------
template <typename MESH, typename QUAD, typename MATERNASSEMBLER,
          typename PRECISION, typename GP>
void statX::inferQ::DriverMaxLikelihood<MESH, QUAD, MATERNASSEMBLER, PRECISION, GP>::
computeLumpedMassVector_( )
{
    // clear nodal quantities
    mesh_->iterateOverNodes( boost::bind( &Node::clearForce, _1 ) );

    // unit lumped mass integrand
    typedef boost::function<void( Element *,
                                  const typename Element::VecLDim &,
                                  const double & )>              Integrand;
    Integrand lumpedMassIntegrand(
            boost::bind( &Element::lumpedMassIntegrand, _1, _2, _3 ) );
    corlib::Integrator<Quad, Integrand>
        lumpedMassIntegrator( quadrature_, lumpedMassIntegrand );
    mesh_->iterateOverElements( lumpedMassIntegrator );

    // assemble lumped mass values into a global vector
    mesh_->iterateOverNodes( corlib::vectorAssemblerFun( &Node::getForce,
                                                         &Node::copyDofArray,
                                                         assembler_ ) );

    // clear nodal quantities after storing in global vector
    mesh_->iterateOverNodes( boost::bind( &Node::clearForce, _1 ) );

    return;
}

//------------------------------------------------------------------------------
template <typename MESH, typename QUAD, typename MATERNASSEMBLER,
          typename PRECISION, typename GP>
void statX::inferQ::DriverMaxLikelihood<MESH, QUAD, MATERNASSEMBLER, PRECISION, GP>::
computePrecisionMatrix_( Eigen::SparseMatrix<double> & qMat,
                         Eigen::SparseMatrix<double> & prMat )
{
    // compute Matern stiffness matrix
    const double kappa = this->extractKappa_( );
    this->computeMaternMatrix_( kappa );

    // get Matern stiffness matrix and lumped mass vector
    Eigen::SparseMatrix<double> maternMat;
    Eigen::VectorXd lumpMassVec;
    assembler_->giveMaternMatrix( maternMat );
    assembler_->giveLumpedMassVector( lumpMassVec );

    // get sparse precision matrix and right fractional matrix
    precision_->givePrecisionMatrix( hps_, maternMat, lumpMassVec, qMat  );
    precision_->giveRightFracMatrix( hps_, maternMat, lumpMassVec, prMat );

    return;
}

//------------------------------------------------------------------------------
template <typename MESH, typename QUAD, typename MATERNASSEMBLER,
          typename PRECISION, typename GP>
double statX::inferQ::DriverMaxLikelihood<MESH, QUAD, MATERNASSEMBLER, PRECISION, GP>::
extractKappa_( ) const
{
    // find \nu value
    auto it = hps_.find( "gaussianProcess.nu" );
    FTL_VERIFY( it != hps_.end( ) );
    const double nu = it->second;

    // find \kappa or \ell value
    it = hps_.find( "gaussianProcess.kappa" );
    double kappa;
    if ( it != hps_.end( ) )
    {
        kappa = it->second;
    }
    else
    {
        auto itr = hps_.find( "gaussianProcess.length" );
        FTL_VERIFY( itr != hps_.end( ) );
        const double length = itr->second;
        kappa = std::sqrt( 2. * nu ) / length;
    }

    return kappa;
}

#endif
