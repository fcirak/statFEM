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

//! @author Jakub Sistek
//! @date   5/2011

/** \brief Inverse partitioner, which creates single SMF file out of set of PSMF files
 * \info Convertor reads a set of PSMF files containing parallel mesh. It creates 
 * a single SMF file with mesh. The SMF is NOT identical with the original SMF file,
 * because nodes are renumbered and elements are placed linearly with respect to processors.
 */
#include <iostream>
#include <fstream>

#include <boost/lexical_cast.hpp>

#include <corlib/UniqueFilename.hpp>
#include <corlib/Shape.hpp>
#include <corlib/ShapeToVtk.hpp>
#include <corlib/verify.hpp>

#include <tools/input/smf2vtu/readwrite.hpp>

using std::cin;
using std::cout;
using std::endl;

//------------------------------------------------------------------------------
int main(int argc, char **argv)
{
    // a word to the user
    if (argc != 3) {
        std::cerr << " Join PSMF files into single SMF file. " << std::endl
                  << " Usage: psmf2smf P basename " << std::endl
                  << " Arguments: " << std::endl
                  << "       P        -- number of subdomains " << std::endl
                  << "       basename -- will search for files basename.isub.psmf file with mesh in PSMF " << std::endl
                  << std::endl;
        return -1;
    }

    // read command line arguments
    const unsigned    numSub   = boost::lexical_cast<int>( argv[1] );
    const std::string basename = argv[2];

    // prepare memory for storing pressure nodes
    std::set <unsigned> pressureNodesSet;

    // prepare memory for starts of subdomain node indices
    std::vector <unsigned> subNodeStarts( numSub + 1, 0 );
    // prepare memory for starts of subdomain elements indices
    std::vector <unsigned> subElementStarts( numSub + 1, 0 );

    // global mesh description
    std::map< unsigned, boost::numeric::ublas::bounded_vector<double,3> > coordinatesGlobMap;
    std::map< unsigned, std::vector<unsigned> >           connectivityGlobMap;

//    corlib::SmfHead smfHead;

    enum corlib::shape eleShape;

    // run loop over subdomains
    for ( unsigned iSub = 0; iSub < numSub; iSub++ ) {

        // prepare local variables
        std::vector< boost::numeric::ublas::bounded_vector<double,3> > coordinates;
        std::vector< std::vector<unsigned> >           connectivity;
        std::vector<unsigned>                          nodesGlobal4local;
        std::vector<unsigned>                          elementsGlobal4local;
        std::vector<unsigned>                          nodesOwner;
        std::vector<unsigned>                          elementsOwner;
        
        // read the mesh from PSMF file
        corlib::UniqueFilename filenameGen( 3 );
        const std::string psmfFilename = filenameGen( basename, iSub, "psmf" );
        std::ifstream psmf( psmfFilename.c_str() );
        FTL_VERIFY_DESCRIPTIVE( psmf.is_open(),
                                "Could not open psmf file %s\n", psmfFilename.c_str() );
        tools::input::smf2vtu::readSmfFile( psmf, eleShape, coordinates, connectivity, 
                                            &nodesGlobal4local, &elementsGlobal4local, 
                                            &nodesOwner,        &elementsOwner );
        psmf.close();
        
        const unsigned numNodes    = coordinates.size( );
        const unsigned numElements = connectivity.size( );

        const unsigned numNodesUnique    = std::count( nodesOwner.begin(), nodesOwner.end(), 1 );
        const unsigned numElementsUnique = std::count( elementsOwner.begin(), elementsOwner.end(), 1 );

        subNodeStarts[ iSub + 1 ]    = subNodeStarts[ iSub ]    + numNodesUnique;
        subElementStarts[ iSub + 1 ] = subElementStarts[ iSub ] + numElementsUnique;
        //subElementStarts[ iSub + 1 ] = subElementStarts[ iSub ] + numElements;

        // insert corrdinates into the map
        for ( unsigned in = 0; in < numNodes; in++ ) {
            if ( nodesOwner[ in ] == 1 ) {
                coordinatesGlobMap.insert( std::make_pair( nodesGlobal4local[in], coordinates[in] ) );
            }
        }

        // convert elements into global numbering and insert them into the connectivity table
        unsigned ieUnique = 0;
        for ( unsigned ie = 0; ie < numElements; ie++ ) {
            if ( elementsOwner[ ie ] == 1 ) {
                std::vector<unsigned> elementGlob ( connectivity[ie].size(), 0 );
                for ( unsigned n = 0; n < connectivity[ie].size( ); n ++ ) {
                    elementGlob[n] = nodesGlobal4local[ connectivity[ie][n] ];
                }
                connectivityGlobMap.insert( std::make_pair( subElementStarts[ iSub ] + ieUnique, elementGlob ) );
                ieUnique++;
                //connectivityGlobMap.insert( std::make_pair( elementsGlobal4local[ie], elementGlob ) );
            }
        }
    }

    // global number of nodes
    unsigned numNodesGlob    = subNodeStarts[ numSub ];

    // global number of elements
    unsigned numElementsGlob = subElementStarts[ numSub ];

    // check sizes
    FTL_VERIFY_DESCRIPTIVE( numNodesGlob ==  coordinatesGlobMap.size(),
                            "Size mismatch in nodes  %d %d \n", 
                            numNodesGlob, coordinatesGlobMap.size() );
    FTL_VERIFY_DESCRIPTIVE( numElementsGlob ==  connectivityGlobMap.size(),
                            "Size mismatch in elements  %d %d \n", 
                            numElementsGlob, connectivityGlobMap.size() );

    // convert map of coordinates into vector
    std::vector< boost::numeric::ublas::bounded_vector<double,3> > coordinatesGlob;
    coordinatesGlob.reserve( coordinatesGlobMap.size() );
    std::map< unsigned, boost::numeric::ublas::bounded_vector<double,3> >::const_iterator iter1 = coordinatesGlobMap.begin();
    for ( ; iter1 != coordinatesGlobMap.end(); ++iter1 ) {
        coordinatesGlob.push_back( iter1 -> second );
    }
    coordinatesGlobMap.clear( );

    // convert map of elements into vector
    std::vector< std::vector<unsigned> > connectivityGlob;
    connectivityGlob.reserve( connectivityGlobMap.size() );
    std::map< unsigned, std::vector<unsigned> >::const_iterator iter2 = connectivityGlobMap.begin();
    for ( ; iter2 != connectivityGlobMap.end(); ++iter2 ) {
        connectivityGlob.push_back( iter2 -> second );
    }

    // write the joined SMF file
    const std::string smfFilename = basename + "_joined.smf";
    std::ofstream smf( smfFilename.c_str() );
    tools::input::smf2vtu::writeSmfFile( smf, eleShape, coordinatesGlob, connectivityGlob );
    smf.close();

    cout << " Conversion O.K., file " << smfFilename << " created or overwritten." << endl;

    return 0;
}
