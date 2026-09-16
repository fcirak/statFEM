// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//
//                 Computational Structural Mechanics Lab
//                         University of Cambridge
//
// This software is copyrighted and all rights are retained by the University of
// Cambridge. All use, disclosure, and/or reproduction of any part not expressly
// authorised by the University of Cambridge is prohibited. (C) 2010
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

//! @file mcmcGaussian.cpp
#include <fstream>
#include <cmath>
#include <vector>

#include <corlib/PropertiesParser.hpp>

#include <statX/infer/MCMC.hpp>

//------------------------------------------------------------------------------
class NormalPDF {
public:
    NormalPDF ( double mu, double stdDev ) : mu_( mu ), stdDev_( stdDev ) { }

    double giveValue ( const double & x ) const {
        const double val = 1. / ( std::sqrt( 2. * M_PI ) * stdDev_ ) *
                std::exp( - ( x - mu_ ) * ( x - mu_ ) / ( 2. * stdDev_ * stdDev_ ) );
        return val;
    }

private:
    const double mu_;
    const double stdDev_;
};

//------------------------------------------------------------------------------
int main( )
{ 
    //! variables to be read from input file "./input.dat"
    //! predefined distribution parameters
    double mean, stdDev;
    //! MCMC setting parameters
    double initialProp, stdDevProp, burnInRatio;
    unsigned seed, numSamples;
    //! output files
    std::string acceptRatioFileName, sampleFileName, likelihoodFileName;

    //! instantiate and feed properties parser with the variables to read
    corlib::PropertiesParser *prop = new corlib::PropertiesParser;
    prop -> registerPropertiesVar( "mean",            mean                );
    prop -> registerPropertiesVar( "stdDev",          stdDev              );
    prop -> registerPropertiesVar( "initialProp",     initialProp         );
    prop -> registerPropertiesVar( "stdDevProp",      stdDevProp          );
    prop -> registerPropertiesVar( "seed",            seed                );
    prop -> registerPropertiesVar( "numSamples",      numSamples          );
    prop -> registerPropertiesVar( "burnInRatio",     burnInRatio         );
    prop -> registerPropertiesVar( "sampleFile",      sampleFileName      );
    prop -> registerPropertiesVar( "acceptRatioFile", acceptRatioFileName );
    prop -> registerPropertiesVar( "likelihoodFile",  likelihoodFileName  );

    //! read variables from the input.dat file
    std::ifstream inputFile( "./input.dat" );
    prop -> readValues( inputFile );

    //! delete properties parser and close the input file
    delete prop;
    inputFile.close( );

    //! predefined univariate normal pdf
    NormalPDF normal( mean, stdDev );
    
    //! instantiate MCMC
    typedef statX::infer::MCMC<1>     MCMC;
    typedef MCMC::VecHpar            VecHpar;
    
    //! lambda function for changing to log
    //! Note that our MCMC implementation expects the log of the pdf 
    auto pdf = [ & normal ]( VecHpar x ) {
        return std::log( normal.giveValue( x[0] ) );
    };
    
    //! construct MCMC
    MCMC mcmc( pdf, seed );
    VecHpar initialPropVec { initialProp };
    VecHpar stdDevPropVec { stdDevProp };
    mcmc.sample( numSamples, stdDevPropVec, initialPropVec );

    //! precision
    const unsigned precision = 10;

    //! write acceptance ratio
    std::ofstream acceptRatioFile( acceptRatioFileName );
    acceptRatioFile << std::setprecision( precision ) 
                    << mcmc.giveAcceptanceRatio( ) << '\n';
    acceptRatioFile.close();

    //! write samples
    auto samples = mcmc.giveSamples( );
    std::ofstream sampleFile( sampleFileName );
    FTL_VERIFY( sampleFile.good() );
    std::ofstream likelihoodFile( likelihoodFileName );
    FTL_VERIFY( likelihoodFile.good() );
    for ( unsigned i = (unsigned) (burnInRatio*numSamples); i < numSamples; ++ i ) {
        sampleFile     << std::setprecision( precision ) 
                       << samples[i].first[0] << '\n';
        likelihoodFile << std::setprecision( precision ) 
                       << samples[i].second   << '\n';
    }
    sampleFile.close();
    likelihoodFile.close();

    return EXIT_SUCCESS;
}
