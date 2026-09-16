//
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//
//                 Computational Structural Mechanics Lab
//                         University of Cambridge
//
// This software is copyrighted and all rights are retained by the CSMLab. All
// use, disclosure, and/or reproduction of any part not expressly authorized by
// F Cirak is prohibited. (C) 2012.
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//

//! @file   rcb/grid/partitioner.cpp
//! @author Matija Kecman
//! @date   3/2012
//! @todo   - does everything work for 1-D? writing vtu does not

// standard library includes
#include <iostream>
#include <fstream>
#include <algorithm>
#include <string>

// third-party includes
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/bind.hpp>

#include <corlib/verify.hpp>
#include <corlib/UniqueFilename.hpp>
#include <corlib/Shape.hpp>

#include <tools/input/partitioners/utils/misc.hpp>
#include <tools/input/smf2vtu/readwrite.hpp>

// local includes
#include "PartitionerRCB.hpp"
#include "Parameters.hpp"


//==============================================================================
namespace aux{

    // make ublas namespace available in namespace
    namespace ublas = boost::numeric::ublas;

    //--------------------------------------------------------------------------
    /** Compute the index of an cell (element) given its ID. Using simple modular
     *  arithmetic.
     *
     *  \param[in] numCells Number of cells per direction in the grid
     *  \param[in] id       ID of the cell
     *  \return             Index of the cell
     */
    ublas::bounded_vector<int,3> idToIndex( const ublas::bounded_vector<int,3> & numCells,
                                            const int id )
    {
        ublas::bounded_vector<int,3> idx; idx.clear();
        int s = id;
        for ( int d = 0; d < 3; ++d ) {
            if ( numCells[d] == 0 ) break;
            idx[d] = s % numCells[d];
            s = s / numCells[d];
        }
        return idx;
    }

    //--------------------------------------------------------------------------
    /** Compute the ID of an item based on its index */
    int indexToId( const ublas::bounded_vector<int,3> & numCells,
                   const ublas::bounded_vector<int,3> & index )
    {
        int id = 0, offset = 1;
        for ( int d = 0; d < 3; ++d ) {
            id += offset * index[d];
            if ( numCells[d] > 0 ) offset *= numCells[d];
            else                   break;
        }
        return id;
    }

    //--------------------------------------------------------------------------
    /** Compute the element-wise index of a node based on its element-wise ID */
    ublas::bounded_vector<int,3> localIdToLocalIndex( const unsigned dim,
                                                      const int localId )
    {
        const boost::array<int,2> lineX = {{ 0, 1 }};
        const boost::array<int,4> quadX = {{ 0, 1, 1, 0 }};
        const boost::array<int,4> quadY = {{ 0, 0, 1, 1 }};
        const boost::array<int,8> hexaX = {{ 0, 1, 1, 0,  0, 1, 1, 0 }};
        const boost::array<int,8> hexaY = {{ 0, 0, 1, 1,  0, 0, 1, 1 }};
        const boost::array<int,8> hexaZ = {{ 0, 0, 0, 0,  1, 1, 1, 1 }};

        ublas::bounded_vector<int,3> localIndex;
        localIndex.clear();

        if      ( dim == 1 ) {  // line
            localIndex[0] = lineX[localId];
        }
        else if ( dim == 2 ) {  // quad
            localIndex[0] = quadX[localId];
            localIndex[1] = quadY[localId];
        }
        else if ( dim == 3 ) {  // hex
            localIndex[0] = hexaX[localId];
            localIndex[1] = hexaY[localId];
            localIndex[2] = hexaZ[localId];
        }
        else
            FTL_VERIFY( false );

        return localIndex;
    }

    //--------------------------------------------------------------------------
    /** Print usage instructions of this command
     */
    void printUsage()
    {
        std::cerr << "Recursive Coordinate Bisection (RCB) partitioning of a Cartesian grid\n"
                  << "Usage:\n"
                  << "        ./partitioner <inputfile> <outputfile>\n"
                  << "Output file format:\n"
                  << "   numTotalCellsX [numTotalCellsY [numTotalCellsZ]] numSubDomains \n"
                  << "   subDomain0IndexMinX [subDomain0IndexMinY [subDomain0IndexMinZ]]     \\\n"
                  << "      subDomain0NumCellsX [subDomain0NumCellsY [subDomain0NumCellsZ]]  \\\n"
                  << "      subDomain0BBMinX    [subDomain0BBMinY    [subDomain0BBMinZ]]     \\\n"
                  << "      subDomain0BBMaxX    [subDomain0BBMaxY    [subDomain0BBMaxZ]]\n"
                  << "   subDomain1IndexMinX [subDomain1IndexMinY [subDomain1IndexMinZ]]     \\\n"
                  << "      subDomain1NumCellsX [subDomain1NumCellsY [subDomain1NumCellsZ]]  \\\n"
                  << "      subDomain1BBMinX    [subDomain1BBMinY    [subDomain1BBMinZ]]     \\\n"
                  << "      subDomain1BBMaxX    [subDomain1BBMaxY    [subDomain1BBMaxZ]]\n"
                  << "   ...                 ...                  ...\n";
        return;
    }

    //--------------------------------------------------------------------------
    /** Generate mesh based on grid cells
     *
     *  Example: 3x2 cells grid (with indices) -> 6 element mesh (with IDs)
     *  \verbatim
     *                          bbmax                                  bbmax
     *  0,2-----1,2-----2,2-----3,2            8 ----- 9 ----- 10----- 11
     *   |       |       |       |             |       |       |       |
     *   |  0,1  |  1,1  |  2,1  |             |   3   |   4   |   5   |
     *   |       |       |       |             |       |       |       |
     *  0,1-----1,1-----2,1-----3,1    ->      4 ----- 5 ----- 6 ----- 7 
     *   |       |       |       |             |       |       |       |
     *   |  0,0  |  1,0  |  2,0  |             |   0   |   1   |   2   |
     *   |       |       |       |             |       |       |       |
     *  0,0-----1,0-----2,0-----3,0            0 ----- 1 ----- 2 ----- 3 
     *  bbmin                                  bbmin
     *  \endverbatim
     *
     *  \param[in]    dim           Number of space dimensions
     *  \param[in]    bbmin         Lower corner of cuboid domain
     *  \param[in]    bbmax         Upper corner of cuboid domain
     *  \param[in]    numCells      Number of cells in each space dimension
     *  \param[in]    cellSizes     Cell sizes in each space dimension
     *  \param[out]   coordinates   Coordinates for each mesh node
     *  \param[out]   connectivity  Element connectivities, ie node IDs for each element
     */
    void generateMesh( const unsigned dim,
                       const ublas::bounded_vector<double,3> & bbmin,
                       const ublas::bounded_vector<double,3> & bbmax,
                       const ublas::bounded_vector<int,3> & numCells,
                       const ublas::bounded_vector<double,3> & cellSizes,
                       std::vector<ublas::bounded_vector<double,3> > & coordinates,
                       std::vector<std::vector<unsigned> > & connectivity )
    {
        typedef ublas::bounded_vector<int,3> VecI3;

        // number of nodes in each direction
        VecI3 numNodes = numCells + ublas::scalar_vector<int>( 3, 1 );
        std::fill( numNodes.begin()+dim, numNodes.end(), 0 );

        // total number of cells and nodes
        const int numTotalCells =
            std::accumulate( numCells.begin(), numCells.begin()+dim, 1,
                             std::multiplies<int>() );
        const int numTotalNodes =
            std::accumulate( numNodes.begin(), numNodes.begin()+dim, 1,
                             std::multiplies<int>() );

        // make co-ordinates of nodes
        coordinates.resize( numTotalNodes );
        for ( int p = 0; p < numTotalNodes; ++p ) {
            const VecI3 index = aux::idToIndex( numNodes, p );
            coordinates[p] = bbmin + ublas::element_prod( index, cellSizes );
        }

        // number of nodes per element
        VecI3 numElementNodes = ublas::scalar_vector<int>( 3, 2 );
        std::fill( numElementNodes.begin()+dim, numElementNodes.end(), 0 );
        const int numNodesPerElement =
            std::accumulate( numElementNodes.begin(), numElementNodes.begin()+dim, 1,
                             std::multiplies<int>() );  // 2^dim (generic)

        // element connectivity
        connectivity.resize( numTotalCells );
        for ( int c = 0; c < numTotalCells; ++c ) {
            connectivity[c].resize( numNodesPerElement );
            const VecI3 cIndex = aux::idToIndex( numCells, c );
            for ( int p = 0; p < numNodesPerElement; ++p ) {
                const VecI3 pIndex = cIndex + aux::localIdToLocalIndex( dim, p );
                const int pId = aux::indexToId( numNodes, pIndex );
                connectivity[c][p] = pId;
            }
        }

    }
}

//==============================================================================
/** Command-line command to do partitioning of a grid
 *  using recursive-coordinate-bisectioning (RCB)
 */
int main( int argc, char * argv[] ){

    // for convenience
    namespace ublas = boost::numeric::ublas;
    namespace rcb   = tools::input::partitioners::rcb;

    std::string inputFile, outputFile;
    if( argc == 3 ){ // standard call
        inputFile  = boost::lexical_cast<std::string>( argv[1] );
        outputFile = boost::lexical_cast<std::string>( argv[2] );
    }
    else{
        aux::printUsage();
        return 0;
    }

    //--------------------------------------------------------------------------
    // open input file
    std::ifstream inp( inputFile.c_str() );
    FTL_VERIFY_DESCRIPTIVE( inp.is_open(), "Cannot open input file %s\n", inputFile.c_str() );

    // read grid parameters
    aux::GridParameters gp( inp );

    //--------------------------------------------------------------------------
    // a few typedefs
    typedef std::vector<int>                              VecI;
    typedef VecI::const_iterator                          VecICIter;
    typedef ublas::bounded_vector<unsigned,3>             VecU3;
    typedef ublas::bounded_vector<int,3>                  VecI3;
    typedef ublas::bounded_vector<double,3>               VecD3;
    typedef std::vector<VecD3>                            VecVecD3;
    typedef std::vector<std::vector<unsigned> >           VecVecU;

    //--------------------------------------------------------------------------
    // create mesh
    VecVecD3 coordinates;
    VecVecU  connectivity;
    aux::generateMesh( gp.dim, gp.bbmin, gp.bbmax, gp.numCells, gp.cellSizes,
                       coordinates, connectivity );
    FTL_VERIFY( coordinates.size() );
    FTL_VERIFY( connectivity.size() );

    //--------------------------------------------------------------------------
    // partition the mesh with RCB
    const unsigned numNodes           = coordinates.size();
    const unsigned numElements        = connectivity.size();
    const unsigned numNodesPerElement = connectivity[0].size();
    rcb::grid::PartitionerRCB partitioner( numNodes, numElements, numNodesPerElement );
    partitioner.createPartition( connectivity, coordinates, gp.numParts );

    // extract element to part mapping
    const VecI elem2part = partitioner.giveElement2PartMapping();

    //--------------------------------------------------------------------------
    // open output file
    std::ofstream output( outputFile.c_str() );
    FTL_VERIFY( output.is_open() );

    // write global number of cells and partitions
    std::copy( gp.numCells.begin(), gp.numCells.begin()+gp.dim,
               std::ostream_iterator<int>( output, " " ) );
    output << "  " << gp.numParts << std::endl;

    // collect vtu files for each partition
    std::vector< std::pair<unsigned,std::string> >  vtuCollection;

    // loop over partitions
    for ( int p = 0; p < gp.numParts; ++p ) {

        // collect element IDs for elements in this part
        VecI inPart;
        for ( VecICIter i = elem2part.begin(); i != elem2part.end(); ++i ) {
            if ( *i == p ) inPart.push_back( std::distance( elem2part.begin(), i ) );
        }

        // compute index of first and last element in partition
        const int minElementId = *std::min_element( inPart.begin(), inPart.end() );
        const int maxElementId = *std::max_element( inPart.begin(), inPart.end() );
        VecI3 indexMin = aux::idToIndex( gp.numCells, minElementId );
        VecI3 indexMax = aux::idToIndex( gp.numCells, maxElementId );

        // grow indicies in the case of overlap
        if ( gp.overlapLayers > 0 ) {
            // increment indicies by overlap. Snap value if new index is beyond extents
            for( unsigned d = 0; d < gp.dim; ++ d ){
                indexMin[d] = std::max( indexMin[d]-gp.overlapLayers, 0 );
                indexMax[d] = std::min( indexMax[d]+gp.overlapLayers, gp.numCells[d]-1 );
            }
        }

        // compute extents of partition
        VecI3 partNumCells = indexMax - indexMin + ublas::scalar_vector<int>( 3, 1 );
        std::fill( partNumCells.begin()+gp.dim, partNumCells.end(), 0 );
        const VecD3 partBBMin = gp.bbmin + ublas::element_prod( indexMin, gp.cellSizes );
        const VecD3 partBBMax = gp.bbmin + ublas::element_prod( indexMax, gp.cellSizes ) + gp.cellSizes;

        // write description of this partition
        std::copy( indexMin.begin(), indexMin.begin()+gp.dim,
                   std::ostream_iterator<int>( output, " " ) );
        output << "  ";
        std::copy( partNumCells.begin(), partNumCells.begin()+gp.dim,
                   std::ostream_iterator<int>( output, " " ) );        
        output << "  ";
        std::copy( partBBMin.begin(), partBBMin.begin()+gp.dim,
                   std::ostream_iterator<double>( output, " " ) );
        output << "  ";
        std::copy( partBBMax.begin(), partBBMax.begin()+gp.dim,
                   std::ostream_iterator<double>( output, " " ) );
        output << std::endl;

        if ( gp.writeVTU ) {
            // generate part mesh
            VecVecD3 partCoords;
            VecVecU  partConnec;
            aux::generateMesh( gp.dim, partBBMin, partBBMax, partNumCells, gp.cellSizes,
                               partCoords, partConnec );

            // label points with part id
            typedef std::pair< std::string, ublas::matrix<double> > PointData;
            PointData partData( "partId",
                                ublas::matrix<double>( partCoords.size(), 1 ) );  // one component
            VecVecD3::iterator pcIter = partCoords.begin();
            for ( ; pcIter != partCoords.end(); ++pcIter )
                partData.second( std::distance( partCoords.begin(), pcIter ), 0 ) = p;
            std::vector<PointData> partDataVec;
            partDataVec.push_back( partData );

            // write vtu
            const std::string partVTU = corlib::UniqueFilename(3)( outputFile, p, "vtu" );
            vtuCollection.push_back( std::make_pair( p, partVTU ) );
            const corlib::shape shape = ( gp.dim == 2 ? corlib::QUADRILATERAL : corlib::HEXAHEDRON );
            tools::input::smf2vtu::writeVtuFile( partVTU, shape, partCoords, partConnec, partDataVec );
        }
    }
    
    if ( gp.writeVTU ) {
        // write the pvd file with collection of vtu files to handle them comfortably in Paraview
        const std::string pvdFilename = outputFile + ".pvd";
        tools::input::smf2vtu::writePvdFile( pvdFilename, vtuCollection );
    }

    // close the output file
    output.close();

    return 0;
}

//------------------------------------------------------------------------------
