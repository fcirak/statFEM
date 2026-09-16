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
//! @date   4/2012

/** \brief Tool for changing mesh of linear elements into quadratic elements
*  \details Reads an SMF file with linear mesh and converts it into mesh of quadratic elements.
*           Currently only supports hexahedra.
*/
#include <iostream>
#include <fstream>

#include <fstream>

//------------------------------------------------------------------------------
#include <string>
#include <vector>
#include <map>
#include <limits>
#include <boost/lexical_cast.hpp>
#include <boost/numeric/ublas/vector.hpp>
//------------------------------------------------------------------------------
#include <corlib/UniqueFilename.hpp>
#include <corlib/Shape.hpp>
#include <corlib/verify.hpp>
//------------------------------------------------------------------------------
#include <boost/timer/timer.hpp>
//------------------------------------------------------------------------------
// tools for simple SMF file handling
#include <tools/input/smf2vtu/readwrite.hpp>

void printUsage();

//! returns true if first multiindex is smaller than the second
template< unsigned DIM >
class CompareIndices
{
    public:
    typedef boost::numeric::ublas::bounded_vector<unsigned,DIM> Indices;
    bool operator() ( Indices first, Indices second ) const {
    
        for ( unsigned ind = 0; ind < DIM; ind ++ ) {
            if      ( first(ind) < second(ind) ) {
                return true;
            }
            else if ( first(ind) > second(ind) ) {
                return false;
            }
        }
        return false;
    }
};

// Mapping of local vertex indices to edges of a hexahedron
template< typename INDS> 
void getEdgeIndicesHexa(unsigned edgeIndex, INDS & Indices)
{
    switch (edgeIndex) {
        case 0:
            Indices[0] = 0;
            Indices[1] = 1;
            break;
        case 1:
            Indices[0] = 1;
            Indices[1] = 2;
            break;
        case 2:
            Indices[0] = 2;
            Indices[1] = 3;
            break;
        case 3:
            Indices[0] = 3;
            Indices[1] = 0;
            break;
        case 4:
            Indices[0] = 4;
            Indices[1] = 5;
            break;
        case 5:
            Indices[0] = 5;
            Indices[1] = 6;
            break;
        case 6:
            Indices[0] = 6;
            Indices[1] = 7;
            break;
        case 7:
            Indices[0] = 7;
            Indices[1] = 4;
            break;
        case 8:
            Indices[0] = 0;
            Indices[1] = 4;
            break;
        case 9:
            Indices[0] = 1;
            Indices[1] = 5;
            break;
        case 10:
            Indices[0] = 2;
            Indices[1] = 6;
            break;
        case 11:
            Indices[0] = 3;
            Indices[1] = 7;
            break;
        default:
            FTL_VERIFY_DESCRIPTIVE( false, "Edge index out of bounds for hexahedra (min 0, max 11): %d \n", edgeIndex );
    }
}

// Mapping of local vertex indices to faces of a hexahedron
template< typename INDS> 
void getFaceIndicesHexa(unsigned faceIndex, INDS & Indices)
{
    switch (faceIndex) {
        case 0:
            Indices[0] = 0;
            Indices[1] = 1;
            Indices[2] = 2;
            Indices[3] = 3;
            break;
        case 1:
            Indices[0] = 4;
            Indices[1] = 5;
            Indices[2] = 6;
            Indices[3] = 7;
            break;
        case 2:
            Indices[0] = 0;
            Indices[1] = 1;
            Indices[2] = 5;
            Indices[3] = 4;
            break;
        case 3:
            Indices[0] = 1;
            Indices[1] = 2;
            Indices[2] = 6;
            Indices[3] = 5;
            break;
        case 4:
            Indices[0] = 3;
            Indices[1] = 2;
            Indices[2] = 6;
            Indices[3] = 7;
            break;
        case 5:
            Indices[0] = 0;
            Indices[1] = 3;
            Indices[2] = 7;
            Indices[3] = 4;
            break;
        default:
            FTL_VERIFY_DESCRIPTIVE( false, "Face index out of bounds for hexahedra (min 0, max 5): %d \n", faceIndex );
    }
}

//=================================
int main( int argc, char * argv[] )
{
    // number of arguments
    const int numArgs = argc;

    // initialize timing
#ifdef VERBOSE
    boost::timer::cpu_timer totalTimer, partTimer;
    totalTimer.start();
    partTimer.start();
#endif

    std::string smfFilename;
    std::vector<std::string> listFiles;
    if ( numArgs >= 2 ) { // standard call
        smfFilename = boost::lexical_cast<std::string>( argv[1] );
        // load files with lists of global nodes to partition
        const unsigned numListFiles = numArgs - 2;
        for ( unsigned i = 0; i < numListFiles; i++ ) {
            listFiles.push_back( std::string( argv[ i + 2 ] ) );
        }
    }
    else {
        printUsage();
    }

    // read the mesh
    std::ifstream smf( smfFilename.c_str() );
    if ( not smf.good() ) {
        std::cout << "Problem with opening file: " << smfFilename << std::endl;
        FTL_VERIFY( false );
    }
    
    namespace ublas = boost::numeric::ublas;
    typedef ublas::bounded_vector<double,3>                VecDim;
    typedef std::vector< VecDim >                          Coordinates;
    typedef std::vector< std::vector<unsigned> >           Elements;

    Coordinates                                            coordinates;
    Elements                                               elements;
    enum corlib::shape                                     eleShape;


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


    // collect edges and faces
    unsigned numNodesPerElementNew;
    unsigned numEdgesPerElement;
    unsigned numFacesPerElement;
    unsigned numNodesPerFace;
    unsigned localEdgeShift;
    unsigned localFaceShift;
    unsigned localCentreShift;
    if ( eleShape == corlib::HEXAHEDRON and numNodesPerElement == 8 ) {
        numNodesPerElementNew = 27;
        numEdgesPerElement    = 12;
        numFacesPerElement    = 6;
        numNodesPerFace       = 4;
        localEdgeShift        = 8;
        localFaceShift        = 20;
        localCentreShift      = 26;
    }
    else {
        std::cout << "Unsupported element type " << corlib::convertShapeEnumToString( eleShape ) << " with " 
                  << numNodesPerElement << " vertices." << std::endl;
        FTL_VERIFY( false );
    }

    // container for all edges in the mesh
    typedef ublas::bounded_vector<unsigned,2>                       EdgePoints;
    typedef std::vector<std::pair<unsigned,unsigned> >              EdgeElements;
    typedef std::map< EdgePoints, EdgeElements, CompareIndices<2> > EdgesMap;
    EdgesMap edges;

    // container for all faces in the mesh
    typedef ublas::bounded_vector<unsigned,4>                       FacePoints;
    typedef std::vector<std::pair<unsigned,unsigned> >              FaceElements;
    typedef std::map< FacePoints, FaceElements, CompareIndices<4> > FacesMap;
    FacesMap faces;

    std::pair<double,unsigned> maxAspectRatio(1.,0);                                      // worst aspect ratio in the mesh
    std::pair<double,unsigned> maxEdgeLength(0.,0);                                       // longest edge in the mesh
    std::pair<double,unsigned> minEdgeLength(std::numeric_limits<double>().infinity(),0); // shortest edge in the mesh

    // create lists of edges and faces
    for ( unsigned ie = 0; ie < numElements; ie++ ) {

        // insert edges into set
        EdgePoints   ep;
        EdgeElements ee;
        std::pair<EdgesMap::iterator, bool> returnedE;

        double elMaxLength = 0.;
        double elMinLength = std::numeric_limits<double>().infinity();

        for ( unsigned iEdge = 0; iEdge < numEdgesPerElement; iEdge++ ) {
            // get local indices
            getEdgeIndicesHexa( iEdge, ep );

            // translate local indices to global
            for ( unsigned iNodeEdge = 0; iNodeEdge < 2; iNodeEdge++ ) {
                ep(iNodeEdge) = (elements[ie])[ep(iNodeEdge)];
            }
            std::sort( ep.begin(), ep.end() );

            // insert edge into map
            returnedE = edges.insert( std::make_pair(ep, ee) );
            // add index of element to inserted edge
            ((*(returnedE.first)).second).push_back(std::make_pair(ie,iEdge));

            // length of edge
            VecDim X0 = coordinates[ep(0)];
            VecDim X1 = coordinates[ep(1)];
            double length = ublas::norm_2( X1 - X0 );

            if ( length > elMaxLength ) elMaxLength = length;
            if ( length < elMinLength ) elMinLength = length;
        }

        // insert faces into set
        FacePoints   fp;
        FaceElements fe;
        std::pair<FacesMap::iterator, bool> returnedF;

        for ( unsigned iFace = 0; iFace < numFacesPerElement; iFace++ ) {
            // get local indices
            getFaceIndicesHexa( iFace, fp );

            // translate local indices to global
            for ( unsigned iNodeFace = 0; iNodeFace < numNodesPerFace; iNodeFace++ ) {
                fp(iNodeFace) = (elements[ie])[fp(iNodeFace)];
            }
            std::sort( fp.begin(), fp.end() );

            // insert face into map
            returnedF = faces.insert( std::make_pair(fp, fe) );
            // add index of element to inserted face
            ((*(returnedF.first)).second).push_back(std::make_pair(ie,iFace));
        }

        // evaluate mesh statistics
        double aspectRatio = elMaxLength / elMinLength;

        if ( elMaxLength > maxEdgeLength.first ) {
            maxEdgeLength = std::make_pair( elMaxLength, ie );
        }
        if ( elMinLength < minEdgeLength.first ) {
            minEdgeLength = std::make_pair( elMinLength, ie );
        }
        if ( aspectRatio > maxAspectRatio.first ) {
            maxAspectRatio = std::make_pair( aspectRatio, ie );
        }
    } 

    // prepare space for indices of higher order nodes
    for ( unsigned ie = 0; ie < numElements; ie++ ) {
        elements[ie].resize( numNodesPerElementNew );
    }

    // start creating higher order nodes behind existing nodes
    unsigned indNode = numNodes;

    // process edges
    for ( EdgesMap::iterator itEdg = edges.begin(); itEdg != edges.end(); ++itEdg ) {

        // create node coordinates
        VecDim X; X.clear();
        for ( unsigned iNodeEdge = 0; iNodeEdge < 2; iNodeEdge++ ) {
            X = X + coordinates[((*itEdg).first)(iNodeEdge)];
        }
        X  = X / 2.;
        coordinates.push_back( X );

        // loop over elements collected at edge and mark the edge in their connectivity
        for ( EdgeElements::iterator itee = (*itEdg).second.begin(); itee != (*itEdg).second.end(); ++itee ) {
            (elements[(*itee).first])[localEdgeShift + (*itee).second] = indNode;
        }

        indNode++;
    }

    // process faces
    for ( FacesMap::iterator itFace = faces.begin(); itFace != faces.end(); ++itFace ) {

        // create node coordinates
        VecDim X; X.clear();
        for ( unsigned iNodeFace = 0; iNodeFace < numNodesPerFace; iNodeFace++ ) {
            X = X + coordinates[((*itFace).first)(iNodeFace)];
        }
        X  = X / double( numNodesPerFace );
        coordinates.push_back( X );

        // loop over elements collected at edge and mark the edge in their connectivity
        for ( FaceElements::iterator itfe = (*itFace).second.begin(); itfe != (*itFace).second.end(); ++itfe ) {
            (elements[(*itfe).first])[localFaceShift + (*itfe).second] = indNode;
        }

        indNode++;
    }

    // process element centres
    for ( unsigned ie = 0; ie < numElements; ie++ ) {

        // create node coordinates
        VecDim X; X.clear();
        for ( unsigned iNodeEl = 0; iNodeEl < numNodesPerElement; iNodeEl++ ) {
            X = X + coordinates[(elements[ie])[iNodeEl]];
        }
        X  = X / double( numNodesPerElement );
        coordinates.push_back( X );

        // loop over elements collected at edge and mark the edge in their connectivity
        elements[ie][localCentreShift] = indNode;

        indNode++;
    }

    // extract baseFileName from SMF file
    std::size_t pos = smfFilename.find_last_of( "." );
    std::string baseFileName = smfFilename.substr( 0,  pos );
    std::string newSmfFilename = baseFileName + "_midsides.smf";

    // write the high order mesh
    std::ofstream smfOut( newSmfFilename.c_str() );

    tools::input::smf2vtu::writeSmfFile( smfOut, eleShape, coordinates, elements );

    smfOut.close();
    std::cout << " Created file: " << newSmfFilename << std::endl;

    std::vector<std::string>::iterator itFiles  = listFiles.begin();
    std::vector<std::string>::iterator itFilesE = listFiles.end();

    // array of fixed variables
    std::vector<bool> isInList(numNodes);

    for ( ; itFiles != itFilesE; ++itFiles ) {

        // export the set of listed nodes 
        std::size_t posl = (*itFiles).find_first_of( "." );
        std::string baseListFileName = (*itFiles).substr( 0,  posl );
        std::string suffListFileName = (*itFiles).substr( posl, (*itFiles).size() );
        std::string newListFilename  = baseListFileName + "_midsides" + suffListFileName;

        // convert the list of nodes into bool array
        std::fill( isInList.begin(), isInList.end(), false );
        std::ifstream ilist( (*itFiles).c_str() );
        if ( not ilist.good() ) {
            std::cout << "Problem with opening file: " << *itFiles << std::endl;
            FTL_VERIFY( false );
        }
        unsigned numListedNodes;
        ilist >> numListedNodes;
        unsigned indListed;
        for ( unsigned i = 0; i < numListedNodes; i++ ) {
            ilist >> indListed;

            isInList[indListed] = true;
        }
        ilist.close();

        // start testing higher order nodes to find if they are in the list
        indNode = numNodes;

        std::vector<unsigned> additionalList;
        // process edge nodes
        for ( EdgesMap::iterator itEdg = edges.begin(); itEdg != edges.end(); ++itEdg ) {

            // test endpoints - if both are listed, list also the midside
            bool wholeEdge = true;
            for ( unsigned iNodeEdge = 0; iNodeEdge < 2; iNodeEdge++ ) {
                wholeEdge = wholeEdge and isInList[((*itEdg).first)(iNodeEdge)];
            }

            if (wholeEdge) {
                additionalList.push_back(indNode);
            }

            indNode++;
        }

        // process faces
        for ( FacesMap::iterator itFace = faces.begin(); itFace != faces.end(); ++itFace ) {

            // test face vertices - if all are listed, list also the face centre
            bool wholeFace = true;
            for ( unsigned iNodeFace = 0; iNodeFace < numNodesPerFace; iNodeFace++ ) {
                wholeFace = wholeFace and isInList[((*itFace).first)(iNodeFace)];
            }

            if (wholeFace) {
                additionalList.push_back(indNode);
            }

            indNode++;
        }


        std::ofstream olist( newListFilename.c_str() );

        unsigned numListedNodesNew = numListedNodes + additionalList.size();
        olist << numListedNodesNew << std::endl;

        // transcript of the existing list
        ilist.open( (*itFiles).c_str() );
        ilist >> numListedNodes;
        for ( unsigned i = 0; i < numListedNodes; i++ ) {
            ilist >> indListed;
            olist << indListed << std::endl;
        }
        ilist.close();

        for ( unsigned i = 0; i < additionalList.size( ); i++ ) {
            olist << additionalList[i] << std::endl;
        }

        olist.close();
        std::cout << " Created file: " << newListFilename << std::endl;

    }

    std::cout << " ======================================" << std::endl;
    std::cout << " Mesh statistics: " << std::endl;
    std::cout << "  worst aspect ratio is          " << maxAspectRatio.first << " for element " << maxAspectRatio.second << std::endl;
    std::cout << "  length of the LONGEST edge is  " << maxEdgeLength.first <<  " at element  " << maxEdgeLength.second << std::endl;
    std::cout << "  length of the SHORTEST edge is " << minEdgeLength.first <<  " at element  " << minEdgeLength.second << std::endl;
    std::cout << "  i.e. LONGEST/SHORTEST ratio is " << maxEdgeLength.first / minEdgeLength.first << std::endl;
    std::cout << " ======================================" << std::endl;
    std::cout << " Midsides O.K. " << std::endl;
#ifdef VERBOSE
    std::cout << " Elapsed time: " << totalTimer.elapsed().wall * 1e-9 << std::endl;
#endif
    
    return 0;
}

//------------------------------------------------------------------------------
void printUsage()
{
    std::cerr << " SMF - Midsides " << std::endl << std::endl
              << " Usage: " << std::endl
              << "       ./smfMidsides mesh.smf [nodeListFile1 nodeListFile2 ... ]" << std::endl
              << " Arguments: " << std::endl
              << "       mesh.smf -- will search for file mesh.smf with mesh in SMF " << std::endl
              << "  nodeListFileX -- will search for file basename.nodeListFileX with list of " << std::endl
              << "                   global nodes and partition it to new global numbering " << std::endl

              << std::endl;
    exit( -1 );
    return;
}

