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

#ifndef del2_apps_reference_laplaceInverse_inference_h
#define del2_apps_reference_laplaceInverse_inference_h

#include <vector>
#include <map>
#include <string>
#include <fstream>
#include <cmath>
#include <functional>
#include <utility>
#include <random>
#include <algorithm>

#include <corlib/VectorAssembler.hpp>
#include <corlib/PropertiesParser.hpp>

#include <statX/infer/ExpKernel.hpp>
#include <statX/infer/GaussianProcess.hpp>
#include <statX/infer/ForwardGaussLHS.hpp>
#include <statX/infer/MismatchCovariance.hpp>
#include <statX/infer/LikelihoodKOH.hpp>
#include <statX/infer/MCMC.hpp>

#include "ioHelpers.hpp"
#include "hyperParamLearning.hpp"

namespace infer  = statX::infer;
namespace eigenX = corlib::eigenX;


//------------------------------------------------------------------------------
/**
 * Bayesian inference of the material coefficients and sensor uncertainty given
 * the observations of potentials (inverse problem). The posterior of hyperparameters
 * p(w|y) are sampled through MCMC sampling. For each iteration, a continuous material
 * property is assumed as the mean of a Gaussian process with sampled coefficients as
 * anchor points. Refer to Girolami, et.al. (2020). The statistical finite element method
 * (statFEM) for coherent synthesis of observation data and model predictions.
 **/
template<typename Mesh, typename Quad, typename SysSolve, typename HpsValueMap, typename GaussProc>
void inference( int numDofs, Mesh & mesh, Quad & quadrature,
                SysSolve * sysmat, const std::string & inferenceInputProps,
                HpsValueMap & hParamMap, GaussProc & gaussProc) {

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

    //! Lambda function to interpolate material property at the element's centroid
    auto setConductivityVar = [ & gaussProc ] ( auto * element ) {
        //! compute the centroid
        typename Element::NodeConstIterator first = element->nodesBegin( );
        typename Element::NodeConstIterator last  = element->nodesEnd( );

        unsigned numNodes( 0 );
        VecDim centroid{ 0. };
        for ( ; first != last; ++first ) {
            centroid += (*first) -> giveCoordinates( );
            numNodes ++;
        }
        centroid /= ( double ) numNodes;

        const double interpConductivity = gaussProc.predictMean( centroid );
        element -> setConductivity( interpConductivity );

    };

    std::vector< unsigned > dofIndices( numDofs );
    std::iota( dofIndices.begin( ), dofIndices.end(), 0 );

    auto giveFEsolution = [ & mesh, & quadrature, & sysmat, & dofIndices,
                            & setConductivityVar, & gaussProc ]
              ( HpsValueMap & hParam, Eigen::VectorXd & feSolution ){
        //! clean the matrix system
        sysmat -> clearMatrix( );
        sysmat -> clearRhs( );

        gaussProc.setObservationValues( hParam );
        mesh.iterateOverElements( std::bind( setConductivityVar, std::placeholders::_1 ) );
        fem( mesh, quadrature, sysmat );
        sysmat -> solveSystem( );

        feSolution.resize( dofIndices.size( ) );
        sysmat -> giveSolution( dofIndices, feSolution );
    };

    //! read the projection matrix from a file
    Eigen::MatrixXd projMat;
    apps::readMatrixPhi( phiMatFile, projMat, numDofs );
    const unsigned numSensor = projMat.cols( );

    //! create an instance of class ForwardGaussLHS
    typedef infer::ForwardGaussLHS ForwardGaussLHS;
    ForwardGaussLHS forwardGaussLHS( giveFEsolution, projMat );

    //! compute finite element solution in the first instance
    forwardGaussLHS.updateSolution( hParamMap );

    //! read sensor coordinates from a file
    std::vector<VecDim> sensorCoord;
    apps::readCoordData<VecDim>( sensorPosFile, sensorCoord );

    //! create an instance of class MismatchCov
    typedef infer::ExpKernel<VecDim>                  ExpKernel;
    typedef infer::MismatchCovariance<ExpKernel>      MismatchCov;
    MismatchCov mismatchCov( sensorCoord );

    //! the class for Kennedy O'Hagan inference and its instance
    typedef infer::LikelihoodKOH<ForwardGaussLHS, MismatchCov>  LikelihoodKOH;
    LikelihoodKOH likelihoodKOH( & forwardGaussLHS, & mismatchCov );

    //! read sample from file
    Eigen::MatrixXd sensorValues;
    apps::readMatrixData( sensorValFile, sensorValues );

    //! set the sampled sensor values to the inference module
    likelihoodKOH.setSensorValues( sensorValues );

    //! hyperparameter learning through MCMC sampling of \f$ p( w | y ) \f$
    hyperParamLearning<LikelihoodKOH, HpsValueMap>( likelihoodKOH, inferenceInputProps, hParamMap );

}

#endif
