// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//
//                 Computational Structural Mechanics Lab
//                         University of Cambridge
//
// This software is copyrighted and all rights are retained by the University of
// Cambridge. All use, disclosure, and/or reproduction of any part not expressly
// authorised by the University of Cambridge is prohibited. (C) 2010
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

//! @file readwrite.hpp

/** \brief Set of I/O functions common to smf2vtu and psmf2vtu convertors.
 */
//------------------------------------------------------------------------------


#ifndef tools_input_smf2vtu_readwrite_h
#define tools_input_smf2vtu_readwrite_h

#include <set>
#include <vector>

#include <boost/numeric/ublas/vector.hpp>

#include <corlib/ShapeToVtk.hpp>
#include <corlib/SmfHead.hpp>

//==============================================================================
// declarations
namespace tools{
    namespace input{
        namespace smf2vtu{

            namespace ublas = boost::numeric::ublas;
            
            void readSmfFile( std::istream & smf,
                              corlib::shape & eleShape,
                              std::vector< ublas::bounded_vector<double,3> > & coordinates,
                              std::vector< std::vector<unsigned> > & connectivity,
                              std::vector<unsigned> * nodesGlobal4local = NULL,
                              std::vector<unsigned> * elementsGlobal4local = NULL,
                              std::vector<unsigned> * nodesOwner = NULL,
                              std::vector<unsigned> * elementsOwner = NULL );

            void readPointDataFile( const std::string & pdatFilename,
                                    std::vector< std::pair< std::string, ublas::matrix< double > > > & pointData );

            void writeSmfFile( std::ostream & smf,
                               const corlib::shape eleShape,
                               const std::vector< ublas::bounded_vector<double,3> > & coordinates, 
                               const std::vector< std::vector<unsigned> > & connectivity, 
                               const std::vector<unsigned> * nodesGlobal4local = NULL,
                               const std::vector<unsigned> * elementsGlobal4local = NULL,
                               const std::vector<unsigned> * nodesOwner = NULL,
                               const std::vector<unsigned> * elementsOwner = NULL );
            
            void writeVtuFile( const std::string & vtuFilename, 
                               const corlib::shape eleShape,
                               const std::vector< ublas::bounded_vector<double,3> > & coordinates,
                               const std::vector< std::vector<unsigned> > & connectivity,
                               const std::vector< std::pair< std::string, ublas::matrix< double > > > & pointData,
                               int subIndex = -1 );
            
            void writePvdFile( const std::string & pvdFilename,
                               const std::vector< std::pair< unsigned, std::string > > & vtuCollection );
            
        }
    }
}

//==============================================================================
// definitions

//------------------------------------------------------------------------------
/** Read mesh from an SMF, PSMF file */
void tools::input::smf2vtu::readSmfFile( std::istream & smf,
                                         corlib::shape & eleShape,
                                         std::vector< ublas::bounded_vector<double,3> > & coordinates,
                                         std::vector< std::vector<unsigned> > & connectivity,
                                         std::vector<unsigned> * nodesGlobal4local,
                                         std::vector<unsigned> * elementsGlobal4local,
                                         std::vector<unsigned> * nodesOwner,
                                         std::vector<unsigned> * elementsOwner )
{
    // try to read header
    corlib::SmfHead smfHead;
    smfHead.read( smf );

    if ( not smfHead.foundHeader( ) )
        FTL_VERIFY_DESCRIPTIVE( false, "SMF file without header! Aborting...");        
    
    unsigned numNodesPerElement = smfHead.giveElementNumPoints();    
    eleShape = smfHead.giveElementShape( );

#ifdef VERBOSE
    std::cout << "smf with header." << std::endl;
    std::cout << "    Number of nodes per element = " << numNodesPerElement << std::endl;
    std::cout << "    Shape of element = " << corlib::convertShapeEnumToString( eleShape ) << std::endl;
#endif

    // read node and element sums
    unsigned numNodes, numElements;

    // file SMF without local number of nodes demanded
    smf >> numNodes >> numElements;
    smf.ignore( std::numeric_limits<std::streamsize>::max(), '\n' );

    // read coordinates
    for ( unsigned n = 0; n < numNodes; n ++ ) {
        boost::numeric::ublas::bounded_vector<double,3> node;
        for ( unsigned d = 0; d < 3; d ++ ) {
            smf >> node[d];
        }
        coordinates.push_back( node );

        if ( nodesGlobal4local ) {
            // for PSMF files, global ID and owner are at the end of node entries
            unsigned globalID;
            smf >> globalID;
            (*nodesGlobal4local).push_back( globalID );
        }
        if ( nodesOwner ) {
            unsigned ownerFlag;
            smf >> ownerFlag;
            (*nodesOwner).push_back( ownerFlag );
        }
        smf.ignore( std::numeric_limits<std::streamsize>::max(), '\n' );
    }

    // read connectivity
    for ( unsigned e = 0; e < numElements; e ++ ) {
        std::vector<unsigned> element( numNodesPerElement );
        for ( unsigned v = 0; v < numNodesPerElement; v ++ ) {
            smf >> element[v];
        }
        connectivity.push_back( element );

        // commented out until it is useful
//        if ( elementsGlobal4local ) {
//            // for PSMF files, global ID and owner are at the end of node entries
//            unsigned globalID;
//            smf >> globalID;
//            (*elementsGlobal4local).push_back( globalID );
//        }
        if ( elementsOwner ) {
            unsigned activity;
            smf >> activity;
            (*elementsOwner).push_back( activity );
        }

        smf.ignore( std::numeric_limits<std::streamsize>::max(), '\n' );
    }
    return;
}

//------------------------------------------------------------------------------
/** Read point data from PDAT file
 *
 * Format of PDAT stream
 * \verbatim
 *      <numVertices> <numDataSets>
 *      <nameFirstDataSet>  <numCompoPerPointOfFirstDataSet>
 *      <nameSecondDataSet> <numCompoPerPointOfSecondDataSet>
 *      ...                 ...
 *      <pointFirstDSCompo0> <pointFirstDSCompo1> <pointFirstDSCompo2> ...
 *      <pointFirstDSCompo0> <pointFirstDSCompo1> <pointFirstDSCompo2> ...
 *      ...                  ...                  ...                  ...
 *      <pointSecondDSCompo0> <pointSecondDSCompo1> ...
 *      <pointSecondDSCompo0> <pointSecondDSCompo1> ...
 *      ...                   ...                   ...
 * \endverbatim
 *
 * \param[in]  pdatFilename  Filename of PDAT
 * \param[out] pointData     Vector of point data and their name
 */
void tools::input::smf2vtu::readPointDataFile( const std::string & pdatFilename,
                                               std::vector< std::pair< std::string, ublas::matrix< double > > > & pointData )
{
    // open smf file
    FTL_VERIFY( pdatFilename != "" );
    std::ifstream pdat( pdatFilename.c_str( ) );
    FTL_VERIFY_DESCRIPTIVE( pdat.is_open( ),
                            "Could not open file %s\n", pdatFilename.c_str( ) );

    // read header
    unsigned numNodes = 0, numData = 0;
    pdat >> numNodes >> numData;
    FTL_VERIFY( numNodes > 0 );
    FTL_VERIFY( numData > 0 );

    // read header - continued
    for ( unsigned d = 0; d < numData; ++d ) {
        std::string dataName = "";
        unsigned dataNumComp = 0;
        pdat >> dataName >> dataNumComp;
        pointData.push_back( std::make_pair( dataName,
                                             ublas::matrix< double >( numNodes, dataNumComp ) ) );
    }

    // loop point data
    for ( unsigned d = 0; d < numData; ++d ) {
        const unsigned dataNumComp = pointData[ d ].second.size2( );
        for ( unsigned n = 0; n < numNodes; n ++ ) {
            for ( unsigned c = 0; c < dataNumComp; c ++ ) {
                double comp = 0.;
                pdat >> comp;
                pointData[ d ].second( n, c ) = comp;
            }
        }
    }

    // done
    pdat.close();
    return;
}

//------------------------------------------------------------------------------
/** Function to write a file in the SMF or PSMF format */
void tools::input::smf2vtu::writeSmfFile( std::ostream & smf, 
                                          const corlib::shape eleShape,
                                          const std::vector< ublas::bounded_vector<double,3> > & coordinates, 
                                          const std::vector< std::vector<unsigned> > & connectivity,
                                          const std::vector<unsigned> * nodesGlobal4local,
                                          const std::vector<unsigned> * elementsGlobal4local,
                                          const std::vector<unsigned> * nodesOwner,
                                          const std::vector<unsigned> * elementsOwner )
{

    FTL_VERIFY( not coordinates.empty( ) );
    FTL_VERIFY( not connectivity.empty( ) );
    
    // write header
    const unsigned numNodesPerElement = connectivity[0].size();
    std::string eleShapeStr = corlib::convertShapeEnumToString( eleShape );
    corlib::SmfHead smfHead;
    smfHead.write( eleShapeStr, numNodesPerElement, smf );

    // write node and element sums
    const unsigned numNodes    = coordinates.size( );
    const unsigned numElements = connectivity.size( );

    smf << numNodes << "  " << numElements << std::endl;

    // write coordinates
    for ( unsigned n = 0; n < numNodes; n ++ ) {
        smf << coordinates[n][0] << "  " << coordinates[n][1] << "  " << coordinates[n][2];
        if ( nodesGlobal4local ) 
            smf << "  " << (*nodesGlobal4local)[n];
        if ( nodesOwner ) 
            smf << "  " << (*nodesOwner)[n];
        smf << std::endl;
    }
    // write elements
    for ( unsigned e = 0; e < connectivity.size(); e ++ ) {
        std::vector<unsigned> element = connectivity[e];
        for ( unsigned n = 0; n < element.size(); n ++ ) {
            smf << element[n] << " ";
        }
        // commented out until it is useful
        //if ( elementsGlobal4local ) 
        //    smf << "  " << (*elementsGlobal4local)[e];
        if ( elementsOwner ) 
            smf << "  " << (*elementsOwner)[e];
        smf << std::endl;
    }
    
    return;
}


//------------------------------------------------------------------------------
/** Write mesh into a VTU file for ParaView */
void tools::input::smf2vtu::writeVtuFile( const std::string & vtuFilename, 
                                          const corlib::shape eleShape,
                                          const std::vector< ublas::bounded_vector<double,3> > & coordinates,
                                          const std::vector< std::vector<unsigned> > & connectivity,
                                          const std::vector< std::pair< std::string, ublas::matrix< double > > > & pointData,
                                          int subIndex )
// subIndex is an optional parameter, if it is present, a CellData block is exported with 
// all cells having constant value (subdomain index), if it is not set or is negative,
// this block is not exported.
// This is often not needed since subdomain number appears as the "vtkCompositeIndex" in ParaView.
{
    FTL_VERIFY( not coordinates.empty( ) );
    FTL_VERIFY( not connectivity.empty( ) );
    
    // open VTU file
    std::ofstream vtuStrm( vtuFilename.c_str() );
    FTL_VERIFY_DESCRIPTIVE( vtuStrm.is_open( ),
                            "Could not open %s\n", vtuFilename.c_str() );

    unsigned numNodesPerElement = (connectivity[0]).size();
    const unsigned vtkType = vtkTypeOfShape( eleShape, numNodesPerElement );


    //--------------------------------------------------------------------------
    // write header
    unsigned nNodes = coordinates.size( );
    unsigned nElements = connectivity.size( );

    vtuStrm << "<?xml version=\"1.0\"?>" << std::endl
            << "<VTKFile type=\"UnstructuredGrid\" byte_order=\"LittleEndian\"> " << std::endl
            << "  <UnstructuredGrid> " << std::endl
            << "    <Piece NumberOfPoints=\"" << nNodes 
            << "\"  NumberOfCells=\"" << nElements << "\">" << std::endl
            << "      <Points>" << std::endl
            << "        <DataArray type=\"Float64\" " 
            << "NumberOfComponents=\"3\" Name=\"Coordinates\" format=\"ascii\">" << std::endl;

    //--------------------------------------------------------------------------
    // write nodes
    for ( unsigned n = 0; n < nNodes; n ++ ) {
        for ( unsigned d = 0; d < 3; d ++ ) {
            vtuStrm << coordinates[n][d] << " ";
        }
        vtuStrm << std::endl;
    }

    //--------------------------------------------------------------------------
    // write node end tags and cell begin tags
    vtuStrm << "        </DataArray>" << std::endl
            << "      </Points> " << std::endl
            << "      <Cells>   " << std::endl
            << "        <DataArray type=\"Int32\" "
            << "Name=\"connectivity\" NumberOfComponents=\"1\" format=\"ascii\">" << std::endl;

    //--------------------------------------------------------------------------
    // read and write cell connectivity    
    for ( unsigned e = 0; e < nElements; e ++ ) {
        for ( unsigned v = 0; v < connectivity[e].size(); v ++ ) {
            // note: vtk does not have all element types (e.g. nine-noded elements)
            if ( v < corlib::vtkTypeNumberOfNodes ( vtkType ) ) {
                vtuStrm << connectivity[e][v] << " ";
            }
        }
        vtuStrm << std::endl;
    }

    //--------------------------------------------------------------------------
    // write offsets and types
    vtuStrm << "        </DataArray>" << std::endl
            << "        <DataArray type=\"Int32\" "
            << " Name=\"offsets\" NumberOfComponents=\"1\" format=\"ascii\">" << std::endl;
    for ( unsigned e = 0; e < nElements; e ++ ) {
        vtuStrm << (e+1) * corlib::vtkTypeNumberOfNodes ( vtkType ) << std::endl;
    }
    vtuStrm << "        </DataArray>" << std::endl
            << "        <DataArray type=\"Int32\" " 
            << "Name=\"types\" format=\"ascii\">" << std::endl;
    for ( unsigned e = 0; e < nElements; e ++ ) {
        vtuStrm << vtkType << std::endl;
    }
    vtuStrm << "        </DataArray> " << std::endl
            << "    </Cells> " << std::endl;

    //--------------------------------------------------------------------------
    // export subdomain index if desired  
    if ( subIndex > -1 ) {
        vtuStrm << "    <CellData> " << std::endl;
        vtuStrm << "        <DataArray type=\"Int32\" "
                << " Name=\"subdomain\" NumberOfComponents=\"1\" format=\"ascii\">" << std::endl;
        for ( unsigned e = 0; e < nElements; e ++ ) {
            vtuStrm << subIndex << std::endl;
        }
        vtuStrm << "        </DataArray>" << std::endl
                << "    </CellData> " << std::endl;
    }

    //--------------------------------------------------------------------------
    // write point data if available
    if ( not pointData.empty() ) {
        // start
        vtuStrm << "    <PointData> " << std::endl;
        // loop point data
        for ( unsigned d = 0; d < pointData.size(); d ++ ) {
            const ublas::matrix< double > & dataField = pointData[ d ].second;
            vtuStrm << "        <DataArray type=\"Float64\""
                    << " Name=\"" << pointData[ d ].first << "\""
                    << " NumberOfComponents=\"" << dataField.size2() << "\""
                    << " format=\"ascii\">" << std::endl;
            for ( unsigned n = 0; n < dataField.size1(); n ++ ) {
                for ( unsigned c = 0; c < dataField.size2(); c ++ ) {
                    vtuStrm << dataField( n, c ) << " ";
                }
                vtuStrm << std::endl;
            }
            vtuStrm << "        </DataArray> " << std::endl;
        }
        // finish
        vtuStrm << "    </PointData> " << std::endl;
    }

    //--------------------------------------------------------------------------
    // close xml file
    vtuStrm << "  </Piece> " << std::endl
            << "</UnstructuredGrid> " << std::endl
            << "</VTKFile> " << std::endl;

    //--------------------------------------------------------------------------
    // close VTU file
    vtuStrm.close();

    return;
}

//------------------------------------------------------------------------------
/** Write collection of filenames into a PVD file to be loaded by ParaView */
void tools::input::smf2vtu::writePvdFile( const std::string & pvdFilename, 
                                          const std::vector< std::pair< unsigned, std::string > > & vtuCollection )
{
    // open pvd file
    std::ofstream pvdStrm( pvdFilename.c_str() );

    FTL_VERIFY_DESCRIPTIVE( pvdStrm.is_open( ),
                            "Could not open %s\n", pvdFilename.c_str() );
    
    // write header
    pvdStrm << "<?xml version=\"1.0\"?>" << std::endl
            << "<VTKFile type=\"Collection\" version=\"0.1\" "
            << "byte_order=\"LittleEndian\">" << std::endl
            << "  <Collection>" << std::endl;
    
    // go through vector of subdomain files
    std::vector< std::pair< unsigned, std::string > >::const_iterator siter = vtuCollection.begin( );
    std::vector< std::pair< unsigned, std::string > >::const_iterator send  = vtuCollection.end( );
    for( ; siter != send; ++ siter ) {
        pvdStrm << "    <DataSet part=\"" << siter -> first 
                << "\" file=\"" << siter -> second << "\"/>" << std::endl;
    }
    
    // finish file
    pvdStrm << "  </Collection>" << std::endl
            << "</VTKFile>" << std::endl;

    
    //--------------------------------------------------------------------------
    // close pvd file
    pvdStrm.close();

    return;
}

#endif
