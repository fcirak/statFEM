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
//! @date   12/2012

//------------------------------------------------------------------------------
#include <iostream>
#include <fstream>
#include <fstream>
#include <string>
#include <vector>
#include <boost/lexical_cast.hpp>
#include <boost/timer/timer.hpp>

#include <corlib/UniqueFilename.hpp>
#include <corlib/Shape.hpp>
#include <corlib/verify.hpp>

// tools for simple file handling
#include <tools/input/smf2vtu/readwrite.hpp>

#include <tools/input/partitioners/utils/readOptions.hpp>
#include <tools/input/partitioners/rcb/mesh/PartitionerRCB.hpp>
#include <tools/input/partitioners/utils/misc.hpp>

//==============================================================================
/** \brief Partitioner of a mesh using RCB 
*  \details This program calls functions of Partitioner RCB to create partition of mesh into
*  subdomains based on geometry. It can run in two distinct modes:
*    1) Non-overlapping partition ( overlapLayers = 0 ) of elements into subdomains - useful for substructuring.
*    2) Overlapping partition based on locality of nodes  ( overlapLayers = 1 )
*       ( this is in general not the same as adding layers to non-overlapping division )
*       This is useful to avoid costly assembly of MPIAIJ matrix in PETSc. One can produce parts large enough to 
*       contain all data required to assembly of whole rows of a distributed matrix. This can be combined with
*       filtering of non-local rows in SystemSolvePetsc.
*
*    
*/
int main( int argc, char * argv[] )
{
    namespace partitioners = tools::input::partitioners;

    // initialize timing
#ifdef VERBOSE
    boost::timer::cpu_timer totalTimer, partTimer;
    totalTimer.start();
    partTimer.start();
#endif

    int numSub = 0;
    unsigned overlapLayers = 0;
    std::string basename;
    std::vector<std::string> listSuffixes;
    partitioners::utils::readOptionsMeshBasedPartitioner( argc, argv,
                                                          numSub, overlapLayers,
                                                          basename, listSuffixes );

    // read the mesh
    const std::string smfFilename = basename + ".smf";
    std::ifstream smf( smfFilename.c_str() );
    FTL_VERIFY( smf.good( ) );
    
    namespace ublas = boost::numeric::ublas;
    std::vector< ublas::bounded_vector<double,3> > coordinates;
    std::vector< std::vector<unsigned> >           elements;
    enum corlib::shape                             eleShape;

    tools::input::smf2vtu::readSmfFile( smf, eleShape, coordinates, elements );

    smf.close();
    FTL_VERIFY( coordinates.size() );
    FTL_VERIFY( elements.size() );
    
#ifdef VERBOSE
    std::cout << "Time for reading mesh: " << partTimer.elapsed().wall * 1e-9 << std::endl;
    partTimer.start();
#endif

    const unsigned numNodes           = coordinates.size();
    const unsigned numElements        = elements.size();
    const unsigned numNodesPerElement = elements[0].size();

    // initialize partitioner
    partitioners::rcb::mesh::PartitionerRCB partitioner( numNodes, numElements, 
                                                         numNodesPerElement );

#ifdef VERBOSE
    std::cout << "Time for construction of partitioner: " << partTimer.elapsed().wall * 1e-9 << std::endl;
    partTimer.start();
#endif

    // partition the mesh by RCB
    partitioner.createPartition( elements, coordinates, numSub, overlapLayers );

#ifdef VERBOSE
    std::cout << "Time for creating partition: " << partTimer.elapsed().wall * 1e-9 << std::endl;
    partTimer.start();
#endif

    // print partitions - you do not want to for large problems!
    //partitioner.printPartitions();

    // Read in files with lists of global nodes - optional, used for prescribing boundary conditions
    std::vector< std::vector< unsigned > > nodeListsGlobal;
    typename std::vector<std::string>::iterator itFiles  = listSuffixes.begin();
    typename std::vector<std::string>::iterator itFilesE = listSuffixes.end();
    for ( ; itFiles != itFilesE; ++itFiles ) {

        // read in list of global nodes from the file
        const std::string fileName = basename + "." + *itFiles;
#ifdef VERBOSE
        std::cout << " Loading file:  " << fileName << std::endl;
#endif

        std::ifstream ins(fileName.c_str() );
        FTL_VERIFY( ins.good( ) );
        
        unsigned length;
        ins >> length;
        std::vector< unsigned > oneListGlobal;
        for ( unsigned i = 0; i < length; i++ ) {
            unsigned nodeId;
            ins >> nodeId;
            oneListGlobal.push_back( nodeId );
        }

        ins.close();
        
        // update list to new node numbering
        std::vector< unsigned > oneListGlobalNew;
        partitioner.updateListNumbering( oneListGlobal, oneListGlobalNew );

        // make sure the list is sorted
        std::sort( oneListGlobalNew.begin(), oneListGlobalNew.end() );

        // debug
        //std::copy( nodeListGlobal.begin(), nodeListGlobal.end(), 
        //           std::ostream_iterator<unsigned>( std::cout, " " ) );
        //std::cout << " \n ";

        nodeListsGlobal.push_back( oneListGlobalNew );
    }

    // create output files for each subdomain
    for ( unsigned iSub = 0; iSub < static_cast<unsigned>( numSub ); iSub++ ) {


        std::vector< ublas::bounded_vector<double,3> > coordinatesSub;
        std::vector<unsigned>                          nodesGlobal4local;
        std::vector<unsigned>                          nodesOwner;

        std::vector< std::vector<unsigned> >           elementsSub;
        std::vector<unsigned>                          elementsGlobal4local;
        std::vector<unsigned>                          elementsOwner;

        // get data for the subdomain from partitioner and write it into files
        partitioner.getSubMesh( coordinates, elements, iSub, 
                                elementsSub, elementsGlobal4local, elementsOwner,
                                coordinatesSub, nodesGlobal4local, nodesOwner );


        // open file for submesh
        corlib::UniqueFilename filenameGen;
        const std::string psmfFilename = filenameGen( basename, iSub, "psmf" );
        std::ofstream psmf( psmfFilename.c_str() );

        // write the submesh file
        tools::input::smf2vtu::writeSmfFile( psmf, eleShape,
                                             coordinatesSub, elementsSub, 
                                             &nodesGlobal4local, 
                                             &elementsGlobal4local, 
                                             &nodesOwner, 
                                             &elementsOwner );
        psmf.close();


        // make output files with the sub-lists of nodes 
        std::vector< std::vector<unsigned> >::iterator itList  = nodeListsGlobal.begin();
        std::vector< std::vector<unsigned> >::iterator itListE = nodeListsGlobal.end();
        for ( ; itList != itListE; ++itList ) {

            // get corresponding list
            const int listIndex = std::distance( nodeListsGlobal.begin(), itList );

            // prepare array for intersection
            std::vector<unsigned> locNodesIntersect( nodeListsGlobal[ listIndex ].size() + nodesGlobal4local.size() ) ;

            // localize the list of nodes to subdomain by comparison of list ( *itList ) with nodesGlobal4local
            std::vector<unsigned>::iterator posEnd;
            posEnd = std::set_intersection( nodeListsGlobal[ listIndex ].begin(), nodeListsGlobal[ listIndex ].end(), 
                                            nodesGlobal4local.begin(), nodesGlobal4local.end(), 
                                            locNodesIntersect.begin() );

            // export the set of fixed nodes 
            const std::string listFileName = filenameGen( basename, iSub, listSuffixes[listIndex] );
            std::ofstream lst( listFileName.c_str() );
#ifdef VERBOSE
            std::cout << " Exporting file:  " << listFileName << std::endl;
#endif

            FTL_VERIFY( lst.good( ) );
            
            lst << std::distance( locNodesIntersect.begin(), posEnd ) << std::endl;
            std::copy( locNodesIntersect.begin(), posEnd, 
                       std::ostream_iterator<unsigned>( lst, " \n" ) );

            lst.close( );

        }
    }

    std::cout << " Partition O.K. " << std::endl;
#ifdef VERBOSE
    std::cout << " Elapsed time: " << totalTimer.elapsed().wall * 1e-9 << std::endl;
#endif
    
    return 0;
}

