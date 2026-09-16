// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//
//                 Computational Structural Mechanics Lab
//                         University of Cambridge
//
// This software is copyrighted and all rights are retained by the University of
// Cambridge. All use, disclosure, and/or reproduction of any part not expressly
// authorised by the University of Cambridge is prohibited. (C) 2010
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

//! @file gaussianProcess.cpp
#include <cmath>
#include <vector>
#include <iostream>
#include <fstream>
#include <string>
#include <Eigen/Dense>

#include <corlib/verify.hpp>
#include <corlib/eigenX.hpp>
#include <corlib/PropertiesParser.hpp>

#include <statX/infer/GaussianProcess.hpp>
#include <statX/infer/ExpKernel.hpp>

namespace eigenX = corlib::eigenX;

//------------------------------------------------------------------------------
//! read data from a given fileName.dat file and store it as Eigen matrix
void readMatrixData ( const std::string & fileName, Eigen::MatrixXd & data )
{
    std::ifstream fileStream( fileName );
    FTL_VERIFY( fileStream.is_open( ) );

    std::vector<double> arrayData;
    unsigned numRows, numCols;
    double x;

    fileStream >> numRows >> numCols;
    while ( fileStream >> x ) {
        arrayData.push_back( x );
    }
    fileStream.close( );
    
    // check the data size fit the matrix
    FTL_VERIFY( arrayData.size() == numRows * numCols );
    
    data = Eigen::Map< Eigen::Matrix< double, Eigen::Dynamic, Eigen::Dynamic > >
        ( &(arrayData[0]), numRows, numCols );

    return;
}

//! read parameters and insert them in a std::map
void readParamToMap ( const std::string & fileName,
                        std::map<std::string, double> & hParam )
{
    std::ifstream fileStream( fileName );
    FTL_VERIFY( fileStream.is_open( ) );

    std::string key;
    double value;
    while ( fileStream >> key >> value ) {
        hParam[ key ] = value; // input them into the map
    }
    fileStream.close( );

    return;
}
//------------------------------------------------------------------------------


//------------------------------------------------------------------------------
int main( )
{ 
    //! spatial dimension of the problem
    constexpr unsigned dim = 1;
    //! input file names
    std::string gpHpsFileName, trainingCoordFileName;
    //! output file names
    std::string likelihoodFileName, predictionCoordFileName, 
                predictionMeanFileName, predictionVarFileName;
    //! instantiate and feed properties parser with the variables to read
    corlib::PropertiesParser *prop = new corlib::PropertiesParser;
    prop -> registerPropertiesVar( "gpHpsFile",           gpHpsFileName           );
    prop -> registerPropertiesVar( "trainingCoordFile",   trainingCoordFileName   );
    prop -> registerPropertiesVar( "likelihoodFile",      likelihoodFileName      );
    prop -> registerPropertiesVar( "predictionCoordFile", predictionCoordFileName );
    prop -> registerPropertiesVar( "predictionMeanFile",  predictionMeanFileName  );
    prop -> registerPropertiesVar( "predictionVarFile",   predictionVarFileName   );

    //! read variables from the input.dat file
    std::ifstream inputFile( "./input.dat" );
    FTL_VERIFY( inputFile.good() );
    prop -> readValues( inputFile );

    //! delete properties parser and close the input file
    delete prop;
    inputFile.close( );

    //! get the parameters
    typedef eigenX::VectorSd<dim>                               VecDim;
    typedef statX::infer::ExpKernel<VecDim>                      ExpKernel;
    typedef typename ExpKernel::HpsValueMap                     HpsValueMap;
    HpsValueMap hyperParams;
    readParamToMap( gpHpsFileName, hyperParams );

    //! read training locations
    Eigen::MatrixXd trainCoords;
    readMatrixData( trainingCoordFileName, trainCoords );

    //! construct Gaussian process (numPts = 6)
    constexpr unsigned numPts = 6;
    FTL_VERIFY( numPts == trainCoords.rows() );
    typedef statX::infer::GaussianProcess<ExpKernel, numPts>     GaussProcess;

    GaussProcess gaussProc;
    gaussProc.setHypIdPre( "gaussianProcess." );
    gaussProc.createGaussianProcess( hyperParams, trainCoords );
    gaussProc.setObservationValues( hyperParams );

    //! precision
    const unsigned outputPrecision = 10;

    //! output files for unit testing created with numPts=6
    //! likelihood of the GP model
    std::ofstream likelihoodFile( likelihoodFileName );
    FTL_VERIFY( likelihoodFile.good() );
    likelihoodFile << std::setprecision(outputPrecision)
                   << gaussProc.giveLogLikelihood() << '\n';
    likelihoodFile.close();

    //! mean and variance at prediction positions
    //! read prediction positions
    std::ifstream predictionCoordFile( predictionCoordFileName );
    FTL_VERIFY( predictionCoordFile.good() );
    std::ofstream predictionMeanFile( predictionMeanFileName );
    FTL_VERIFY( predictionMeanFile.good() );
    std::ofstream predictionVarFile( predictionVarFileName );
    FTL_VERIFY( predictionVarFile.good() );
    double predX;
    //! make predictions
    while ( predictionCoordFile >> predX ) {
        Eigen::VectorXd xStarVec(dim); 
        xStarVec << predX;
        predictionMeanFile << std::setprecision(outputPrecision) 
                           << gaussProc.predictMean( xStarVec )     << '\n';
        predictionVarFile  << std::setprecision(outputPrecision) 
                           << gaussProc.predictVariance( xStarVec ) << '\n';
    }
    predictionCoordFile.close();
    predictionMeanFile.close();
    predictionVarFile.close();

    return EXIT_SUCCESS;
}
