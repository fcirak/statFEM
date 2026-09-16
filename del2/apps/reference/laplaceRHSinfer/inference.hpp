// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//
//                 Computational Structural Mechanics Lab
//                         University of Cambridge
//
// This software is copyrighted and all rights are retained by the University of
// Cambridge. All use, disclosure, and/or reproduction of any part not expressly
// authorised by the University of Cambridge is prohibited. (C) 2010 --
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

//! @file inference.hpp

#ifndef del2_apps_reference_laplaceRHSinfer_inference_h
#define del2_apps_reference_laplaceRHSinfer_inference_h

#include <vector>
#include <map>
#include <string>
#include <fstream>
#include <cmath>
#include <functional>
#include <utility>

#include <corlib/VectorAssembler.hpp>
#include <corlib/PropertiesParser.hpp>

#include <statX/infer/ExpKernel.hpp>
#include <statX/infer/CovarianceForcing.hpp>
#include <statX/infer/ForwardGaussRHS.hpp>
#include <statX/infer/MismatchCovariance.hpp>
#include <statX/infer/LikelihoodKOH.hpp>
#include <statX/infer/MCMC.hpp>

#include "ioHelpers.hpp"
#include "hyperParamLearning.hpp"

namespace infer  = statX::infer;
namespace eigenX = corlib::eigenX;


//------------------------------------------------------------------------------
/**
 * Bayesian inference of the true generating process given the observations using the
 * hierarchical approach. First the posterior of hyperparameters p(w|y) are sampled
 * through MCMC sampling. Point estimates of the hyperparameters are then used to
 * evaluate the finite element posterior density p(u|y) and true process density p(z|y).
 * Refer to: Girolami, et.al. (2020). The statistical finite element method (statFEM)
 * for coherent synthesis of observation data and model predictions.
 **/
template<typename Mesh, typename Quad, typename SysSolve, typename HpsValueMap>
void inference( int numDofs, Mesh & mesh, Quad & quadrature,
                SysSolve * sysmat, const std::string & inferenceInputProps,
                HpsValueMap & hParamMap) {

    typedef typename Mesh::Node       Node;
    typedef typename Mesh::Element    Element;
    typedef typename Node::VecDim     VecDim;

    constexpr auto dim = Node::dim;
    constexpr auto dof = Node::dof;

    //! reading inference parameters
    corlib::PropertiesParser *inferenceProps = new corlib::PropertiesParser;

    //! number of readings
    unsigned numReadings;
    inferenceProps -> registerPropertiesVar( "numReadings",      numReadings     );
    //! input file names
    std::string phiMatFile, sensorValFile, sensorPosFile, conductCoordFile;
    inferenceProps -> registerPropertiesVar( "phiMatFile",       phiMatFile      );
    inferenceProps -> registerPropertiesVar( "sensorValFile",    sensorValFile   );
    inferenceProps -> registerPropertiesVar( "sensorPosFile",    sensorPosFile   );
    inferenceProps -> registerPropertiesVar( "conductCoordFile", conductCoordFile);
    //! output file names
    std::string feCondDataFile, meanCovarRHSFile, meanCovarGenFile;
    inferenceProps -> registerPropertiesVar( "feCondDataFile",   feCondDataFile  );
    inferenceProps -> registerPropertiesVar( "meanCovarRHSFile", meanCovarRHSFile);
    inferenceProps -> registerPropertiesVar( "meanCovarGenFile", meanCovarGenFile);

    //! read variables from the inferenceInput.dat file
    std::ifstream inferenceInputFile( inferenceInputProps );
    FTL_VERIFY( inferenceInputFile.is_open( ) );
    inferenceProps -> readValues( inferenceInputFile );

    //! delete properties parser and close the input file
    delete inferenceProps;
    inferenceInputFile.close( );

    //! typedef for CovarForcing
    typedef infer::ExpKernel<VecDim>                           ExpKernel;
    typedef infer::CovarianceForcing<Node::dof, ExpKernel>     CovarForcing;

    //! create an instance of CovarForcing
    CovarForcing covarForcing( mesh.numNodes() );

    //! set the nodal coordinates to covarForcing
    mesh.iterateOverNodes( [ & covarForcing ]( const auto * const node ) {
        const auto coord = node->giveCoordinates( );
        covarForcing.insertToCoordinates( coord );
    } );

    //! compute the integral of the shape function
    typedef eigenX::VectorSd<1>                        Vec1;
    typedef std::function<Vec1( const VecDim )>        VecDim2Vec1;
    typedef corlib::BodyForce<Element,VecDim2Vec1>     UnitForce;

    //! the shape function integral is firstly computed at the nodes
    mesh.iterateOverNodes( [ ]( auto * node ){ node -> clearForce( ); } );
    UnitForce shapeFunIntegrand( [  ]( const VecDim & x ) { eigenX::VectorSd<1> f;
                                    f( 0 ) = 1.; return f; }
    );

    corlib::Integrator<Quad, UnitForce> shapeFunIntegrator( quadrature, shapeFunIntegrand );
    mesh.iterateOverElements( shapeFunIntegrator );

    //! assemble the shape function integral to a vector saved in covarForcing
    mesh.iterateOverNodes( corlib::vectorAssemblerFun( & Node::getForce,
                                                       & Node::copyDofArray,
                                                       & covarForcing ) );

    //! Apply dirichlet boundary conditions in the covarForcing
    corlib::ConstraintFunctor< Node, CovarForcing > constraint( & covarForcing );
    constraint = mesh.iterateOverNodes( constraint );
    constraint.applyConstraints( );

    //! read the projection matrix from a file
    Eigen::MatrixXd projMat;
    apps::readMatrixPhi( phiMatFile, projMat, numDofs );
    const unsigned numSensor = projMat.cols( );

    //! create an instance of forward solver with only rhs uncertainty
    typedef infer::ForwardGaussRHS<SysSolve, CovarForcing>     ForwardGaussRHS;

    //! add regularisation to covariance matrix \f$ C_u \f$
    //! defined as \f$ ( reg * \sigma_f )^2 I \f$
    const double reg = 0.001;
    ForwardGaussRHS forwardGaussRHS( sysmat, & covarForcing, reg );

    //! set projection matrix
    forwardGaussRHS.setProjectionMatrix( projMat );

    //! update the solution in the first instance
    forwardGaussRHS.updateSolution( hParamMap );

    //! read sensor coordinates from a file
    std::vector<VecDim> sensorCoord;
    apps::readCoordData<VecDim>(sensorPosFile, sensorCoord);

    //! create an instance of class MismatchCov
    typedef infer::MismatchCovariance<ExpKernel>      MismatchCov;
    MismatchCov misMatchCov( sensorCoord );

    //! the class for Kennedy O'Hagan inference and its instance
    typedef infer::LikelihoodKOH<ForwardGaussRHS, MismatchCov>  LikelihoodKOH;
    LikelihoodKOH likelihoodKOH( & forwardGaussRHS, & misMatchCov );

    //! read sample from file
    Eigen::MatrixXd sensorValues;
    apps::readMatrixData( sensorValFile, sensorValues );

    //! set the sampled sensor values to the inference module
    likelihoodKOH.setSensorValues( sensorValues );

    //! hyperparameter learning through MCMC sampling of \f$ p( w | y ) \f$
    hyperParamLearning<LikelihoodKOH, HpsValueMap>( likelihoodKOH, inferenceInputProps, hParamMap );

    //! compute finite element mean and covariance given data
    Eigen::VectorXd meanFECondData;
    Eigen::MatrixXd covarFECondData;
    likelihoodKOH.setProjectionMatrix( projMat );
    likelihoodKOH.giveFEconditionedDataWoodbury( hParamMap, meanFECondData, covarFECondData );

    //! print \f$ p(u|y) \f$
    std::ofstream feCondDataStream( feCondDataFile.c_str( ) );
    FTL_VERIFY( feCondDataStream.is_open( ) );
    feCondDataStream << meanFECondData;
    feCondDataStream.close();

}

#endif
