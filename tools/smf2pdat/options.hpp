//
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//
//                 Computational Structural Mechanics Lab
//                         University of Cambridge
//
// This software is copyrighted and all rights are retained by the CSMLab. All
// use, disclosure, and/or reproduction of any part not expressly authorized by
// F Cirak is prohibited. (C) 2010.
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//

//! @authors Kosala Bandara, Jakub Sistek
//! @date   12/2012

/** \brief Convertor to generate a tag file for a mesh
 * \info Compares the coordinates of the given mesh region with that of the full 
 * mesh to identify the tagged node indices.
 */

#ifndef tools_helpers_kkmb2_smf2tg_options_h
#define tools_helpers_kkmb2_smf2tg_options_h

#include <iostream>
#include <algorithm>
#include <iterator>
#include <vector>
#include <string>

#include <boost/program_options.hpp>

//==============================================================================
// declarations
namespace tools{
    namespace helpers{
        namespace smf2tg{
    
            //------------------------------------------------------------------
            /** \brief Read command-line arguments to find input file name etc
             *
             * \param[in]    argc            Number of command-line args
             * \param[in]    argv            Command-line arguments
             * \return                       Status
             */
            int readOptions( int argc,
                             char ** argv,
                             std::string & meshFilename,
                             std::string & regionFilename,
                             std::string & outFilename,
                             double & tol,
                             bool & inverse, bool & writeTag, 
                             bool & writeETag, bool & writeFix,
                             bool beQuiet = false);

        }
    }
}

//==============================================================================
// definitions

//------------------------------------------------------------------------------
// Read command-line arguments to find input file name etc
int tools::helpers::smf2tg::readOptions( int argc,
                                         char ** argv,
                                         std::string & meshFilename,
                                         std::string & regionFilename,
                                         std::string & outFilename,
                                         double & tol,
                                         bool & inverse,
                                         bool & writeTag,
                                         bool & writeETag,  
                                         bool & writeFix,
                                         bool beQuiet)
{
    // Declare the supported options.
    namespace po = boost::program_options;
    po::options_description desc( "Allowed options" );
    desc.add_options( )
        ( "help,h", "Produce help message" )
        ( "inverse,i", "inverse region")
        ( "tag,t", "write tag file")
        ( "edgetag,et", "write edgetag file")
        ( "fix,f", "write fix file")
        ( "quiet,q", "Suppress messages" )
        ( "output-file,o", po::value< std::string >( ), "output file" )
        ( "delta,d", po::value< double >( ), "tolarence" );

    // Declare options which are not written 
    po::options_description hidden( "Hidden options" );
    hidden.add_options( )
        ( "input-files", po::value< std::vector< std::string > >( ), "input files" );

    // All options
    po::options_description cmdlineOptions;
    cmdlineOptions.add( desc ).add( hidden );

    // Arguments without flags
    po::positional_options_description p;
    p.add( "input-files", -1 );

    // parse line
    po::variables_map vm;
    po::store( po::command_line_parser( argc, argv ).
               options( cmdlineOptions ).positional( p ).run( ), vm );
    po::notify( vm );

    // usage message
    std::stringstream usageMessage;
    usageMessage << "Generates pdat,fix,tag files after comparing two smf files" 
                 << std::endl << std::endl
                 << "Usage: smf2pdat [options] fullMesh.smf region.smf" << std::endl
                 << "       If output file name is not given output is written to <fullMesh>.pdat" << std::endl;

    // help
    if ( vm.count( "help" ) or ( vm.size() == 0 ) ) {
        std::cout << usageMessage.str( ) << std::endl
                  << desc << std::endl;
        return 1;
    }

    // inverse
    inverse = vm.count( "inverse" );

    // tag
    writeTag = vm.count( "tag" );

    // edgetag
    writeETag = vm.count( "edgetag" );

    // inverse
    writeFix = vm.count( "fix" );

    // be quiet
    beQuiet = vm.count( "quiet" );

    // input files
    meshFilename = regionFilename = "";
    if ( vm.count( "input-files" ) ) {
        const std::vector< std::string > inputFiles =
            vm[ "input-files" ].as< std::vector< std::string > >( );
        if ( inputFiles.size( ) >= 2 ) {
            meshFilename = inputFiles[ 0 ];
            regionFilename = inputFiles[ 1 ];
        }
        else {
            std::cerr << usageMessage.str( ) << std::endl
                      << desc << std::endl;
            return -1;
        }
    }

    // output file
    outFilename = "";
    if ( vm.count( "output-file" ) ) {
        outFilename = vm[ "output-file" ].as< std::string >( );
    } else 
        outFilename = meshFilename.substr( 0, meshFilename.find( ".smf" ) ) + ".pdat";

    // tolereance
    tol = 1.e-6;
    if ( vm.count( "delta" ) )
        tol = vm[ "delta" ].as< double >( );

    return 0;
}

#endif
