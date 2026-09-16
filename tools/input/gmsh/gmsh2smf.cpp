// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//
//                 Computational Structural Mechanics Lab
//                         University of Cambridge
//
// This software is copyrighted and all rights are retained by the University of
// Cambridge. All use, disclosure, and/or reproduction of any part not expressly
// authorised by the University of Cambridge is prohibited. (C) 2010
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

//! @file gmsh2smf.cpp

#include <iostream>
#include <fstream>
#include <map>
#include <vector>
#include <limits>
#include <algorithm>
#include <set>
#include <iomanip>

#include <boost/tuple/tuple.hpp>

#include <corlib/Shape.hpp>
#include <corlib/SmfHead.hpp>
#include <corlib/verify.hpp>
#include <corlib/eigenX.hpp>

// use longer variables
#define unsigned unsigned long

static const unsigned MAXLINE = std::numeric_limits<std::streamsize>::max();

// helper tyepefs
/// Tuple storing Gmsh's element type, face type together with shape
typedef 
boost::tuple<unsigned, unsigned, enum corlib::shape>            GMSHelemType;
typedef 
std::map<unsigned,std::string>                                  PhysicalNamesMap;
/// Map storing dimension and index of nodes/elements with pyshical entities.
/// This is used for 4.1 gmsh file to identity physical tags
typedef
std::map<std::pair<unsigned, unsigned>, unsigned>               EntityMap;
typedef 
std::vector<corlib::eigenX::VectorSd<3> >                       NodesVec;
typedef
std::vector<std::vector<unsigned> >                             ElementsVec;
typedef 
std::multimap<unsigned,std::vector<unsigned> >                  FacesMMap;
typedef 
std::set<unsigned>                                              MonitoredFacesSet;
typedef 
std::set<unsigned>                                              MonitoredNodesSet;
typedef
std::map<unsigned,std::vector<int long> >                       BCdirectionsMap;
typedef 
std::map<unsigned,std::vector<double> >                         BCvaluesMap;


// function declarations
std::istream & readGMSH( std::istream & inp, 
                         PhysicalNamesMap & physicalNames,
                         NodesVec & nodes,
                         ElementsVec & elements,
                         FacesMMap & faces,
                         unsigned numNodesPerElem,
                         unsigned numNodesPerFace);

std::ostream & writeSMF( std::ostream & smfFile,
                         const NodesVec & nodes, 
                         const ElementsVec & elements,
                         unsigned numNodesPerElem,
                         unsigned numNodesPerFace );

void getConstraints( const PhysicalNamesMap & physicalNames,
                     MonitoredFacesSet & sideEdges,
                     MonitoredNodesSet & monitoredNodes,
                     BCdirectionsMap & directions,
                     BCvaluesMap & values );

void writeMonitoredNodes( const std::string & baseFileName,
                          const PhysicalNamesMap & physicalNames,
                          const MonitoredNodesSet & monitoredNodes,
                          const FacesMMap & faces );

void writeMonitoredFaces( const std::string & baseFileName, 
                          const PhysicalNamesMap & physicalNames,
                          const MonitoredFacesSet & monitoredFaces, 
                          const FacesMMap & faces );

void writeConstraintNodes( const std::string & baseFileName, 
                           const FacesMMap & faces,
                           const BCdirectionsMap & directions,
                           const BCvaluesMap & values);


// helper functions
std::istream & readMeshFormat( std::istream & inp, 
                               std::string & meshFormat );

std::istream & readPhysicalNames( std::istream & inp,
                                  PhysicalNamesMap & physicalNames,
                                  const std::string & meshFormat );

std::istream & readEntitiesVersion4( std::istream & inp,
                                     EntityMap & entityMap );

std::istream & readNodes( std::istream & inp, 
                          NodesVec & nodes,
                          const std::string meshFormat );

// for 2.* version gmsh file
std::istream & readElementsVersion2( std::istream & inp,
                                     ElementsVec & elements,
                                     FacesMMap & faces,
                                     unsigned numNodesPerElem,
                                     unsigned numNodesPerFace);
// for 4.* version gmsh file
std::istream & readElementsVersion4( std::istream & inp,
                                     const EntityMap & entityMap,
                                     ElementsVec & elements,
                                     FacesMMap & faces,
                                     unsigned numNodesPerElem,
                                     unsigned numNodesPerFace );

GMSHelemType gmshElement( unsigned numNodesPerElem, 
                          unsigned numNodesPerFace );

void findOrphans( const unsigned numNodes,
                  const ElementsVec & elements,
                  std::vector<unsigned> & orphans );

void removeOrphans( std::vector<unsigned> & orphans,
                    NodesVec & nodes, 
                    ElementsVec & elements ); 


//==============================================================================
int main (int argc, char *argv []) 
{
    std::string errorMsg = 
        "Usage:  gmsh2smf  -e numNodesPerElem  -f numNodesPerFace  -i inputFile.msh";

    if ( argc != 7) {
        std::cout <<  errorMsg << std::endl;
        return 0;
    }


    unsigned numNodesPerElem, numNodesPerFace;
    std::string gmshFileName;

    char c;
    while ((c = getopt(argc, argv, "e:f:i:")) != EOF)
        switch (c) {
        case 'e' :
            numNodesPerElem = atoi(optarg);
            break;
        case 'f' :
            numNodesPerFace = atoi(optarg);
            break;
        case 'i':
            gmshFileName = optarg;
            break;
        case '?':
            std::cout <<  errorMsg << std::endl;
            return 0;
        }
    
    
    // open mesh file 
    std::ifstream gmshFile( gmshFileName.c_str() );
    FTL_VERIFY( gmshFile.is_open( ) );
    
    // extract baseFileName from gmsh file
    std::size_t pos = gmshFileName.find_last_of( "." );
    std::string baseFileName = gmshFileName.substr( 0,  pos );

    // read and write connectivity and coordinates

    PhysicalNamesMap physicalNames;
    NodesVec nodes;
    ElementsVec elements;
    FacesMMap faces;

    // read gmsh file
    readGMSH( gmshFile, physicalNames, nodes, elements, faces,
              numNodesPerElem, numNodesPerFace);

    gmshFile.close();

    // write smf file 
    std::string smfFileName = baseFileName + ".smf";
    std::ofstream smfFile( smfFileName.c_str() );
    
    writeSMF( smfFile, nodes, elements, numNodesPerElem, numNodesPerFace );
    
    smfFile.close();

    // read and write boundary conditions
    
    MonitoredFacesSet monitoredFaces;
    MonitoredNodesSet monitoredNodes;
    BCdirectionsMap directions;
    BCvaluesMap values;

    // get constraints directions, values  and to be extracted entities from t 
    getConstraints( physicalNames, monitoredFaces, monitoredNodes, directions, values );
    
    // write nodes of faces on which no boundary conditions are described but 
    // need to be identified, for example, for force computations
    writeMonitoredNodes( baseFileName, physicalNames, monitoredNodes, faces );

    // write faces on which no boundary conditions are described but need to be 
    // identified, for example, for force computations
    writeMonitoredFaces( baseFileName, physicalNames, monitoredFaces, faces );
    
    // write constraint nodes with constraint directions and values
    writeConstraintNodes(baseFileName, faces, directions, values );
    
    return 0;
}


//------------------------------------------------------------------------------
void writeConstraintNodes( const std::string & baseFileName, 
                           const FacesMMap & faces,
                           const BCdirectionsMap & directions,
                           const BCvaluesMap & values)
// write constraint nodes with constraint directions and values
{    
    // compress constraints (remove multiple vertices)
    std::map<unsigned,std::set<unsigned> > constraintNodes;

    FacesMMap::const_iterator fIter = faces.begin();
    FacesMMap::const_iterator fEnd  = faces.end();
    
    for ( ; fIter != fEnd; ++ fIter ) {
        const unsigned nBoundary = fIter -> first;

        FacesMMap::mapped_type face = fIter -> second;
        
        for ( unsigned v = 0; v < face.size(); ++ v ) {
            ( constraintNodes[nBoundary] ).insert( face[v] );
        }
    }

    // write constraint file 
    std::string cstrFileName = baseFileName + ".cnstr";
    std::ofstream cstrFile( cstrFileName.c_str() );
    unsigned nConstraints = 0;

    std::map<unsigned,std::set<unsigned> >::iterator cIter = constraintNodes.begin();
    std::map<unsigned,std::set<unsigned> >::iterator cEnd  = constraintNodes.end();
    std::set<unsigned>::iterator sIter, sEnd;

    for ( ; cIter != cEnd; ++ cIter ) {
        const unsigned nBoundary = cIter -> first;

        sIter = (cIter -> second).begin();
        sEnd = (cIter -> second).end();

        BCdirectionsMap::const_iterator itd = directions.find( nBoundary );
        if (itd == directions.end() ) continue;

        std::vector<int long> dir = (*itd).second;
        for ( ; sIter != sEnd; ++ sIter ) {
            for ( unsigned d = 0; d < dir.size(); ++ d ) {
                nConstraints ++;
            }
        }
    }
    cstrFile << nConstraints << std::endl;

    cIter = constraintNodes.begin();
    for ( ; cIter != cEnd; ++ cIter ) {
        const unsigned nBoundary = cIter -> first;

        BCdirectionsMap::const_iterator itd = directions.find( nBoundary );
        if ( itd == directions.end() ) continue;        
        std::vector<int long>    dir = ( *itd).second;
        
        BCvaluesMap::const_iterator itv = values.find( nBoundary );
        if ( itv == values.end() ) continue;
        std::vector<double> val = (*itv).second;

        sIter = (cIter -> second).begin();
        sEnd = (cIter -> second).end();
        for ( ; sIter != sEnd; ++ sIter ) {
            for ( unsigned d = 0; d < dir.size(); ++d ) {
                cstrFile << *sIter << "  " << dir[d] << "  " << val[d] << std::endl;
            }
        }
    }
    cstrFile.close();

    return;
}


//------------------------------------------------------------------------------
void writeMonitoredFaces( const std::string & baseFileName, 
                          const PhysicalNamesMap & physicalNames,
                          const MonitoredFacesSet & monitoredFaces, 
                          const FacesMMap & faces )
// write faces on which no boundary conditions are described but need to be 
// identified, for example, for force computations
{
    MonitoredFacesSet::const_iterator seIter = monitoredFaces.begin();
    MonitoredFacesSet::const_iterator seEnd = monitoredFaces.end();

    for ( ; seIter != seEnd; ++ seIter ) {

        PhysicalNamesMap::const_iterator itphm = physicalNames.find( *seIter );
        FTL_VERIFY_DESCRIPTIVE( itphm != physicalNames.end(),
                                "Cannot find physical name for number=%d\n", *seIter );

        std::cout << "extracting faces of \"" + (*itphm).second << "\"" << std::endl;

        std::string sideFileName = baseFileName;
        sideFileName += "." + (*itphm).second;
        sideFileName += ".edges";
        std::ofstream sideFile( sideFileName.c_str() );

        FacesMMap::const_iterator face = faces.find( *seIter );
        if ( face == faces.end() ) {
            std::cout << "No faces in PhysicalNames " << *seIter
                      << " \"" <<  (*itphm).second  
                      << "\" ... continuing" << std::endl;
            continue;
        }
        FTL_VERIFY( face != faces.end() );

        const unsigned nFaces = faces.count( *seIter );
        sideFile << nFaces << std::endl;
        for ( unsigned n = 0; n < nFaces; ++ n , ++ face ) {
            for ( unsigned v = 0; v < face->second.size(); ++ v ) {
                sideFile << (face->second)[v] << "  ";
            }
            sideFile << std::endl;
        }

        sideFile.close();
    }

    return;
}


//------------------------------------------------------------------------------
void writeMonitoredNodes( const std::string & baseFileName,
                          const PhysicalNamesMap & physicalNames,
                          const MonitoredNodesSet & monitoredNodes,
                          const FacesMMap & faces )
// write nodes of faces on which no boundary conditions are described but 
// need to be identified, for example, for force computations
{
    MonitoredNodesSet::const_iterator snIter = monitoredNodes.begin();
    MonitoredNodesSet::const_iterator snEnd  = monitoredNodes.end();
    
    for ( ; snIter != snEnd; ++ snIter ) {
        
        PhysicalNamesMap::const_iterator itphm = physicalNames.find( *snIter );
        FTL_VERIFY( itphm != physicalNames.end() );

        std::cout << "Extracting nodes of " << (*itphm).second << std::endl;
        
        std::string sideFileName = baseFileName;
        sideFileName += "." + (*itphm).second;
        sideFileName += ".nodes";

        FacesMMap::const_iterator face = faces.find( *snIter );
        if ( face == faces.end() ) {
            std::cout << "No faces in PhysicalNames " << *snIter
                      << " \"" << (*itphm).second 
                      << "\" ... continuing" << std::endl;
            continue;
        }
        FTL_VERIFY( face != faces.end() );   

        const unsigned nFaces = faces.count( *snIter );
 
        std::set<unsigned> nodeIds;

        for ( unsigned n = 0; n < nFaces; ++ n , ++ face ) {
            for ( unsigned v = 0; v < face->second.size(); ++ v ) {
                nodeIds.insert( (face->second)[v] );
            }
        }

        std::ofstream sideFile( sideFileName.c_str() );

        sideFile << nodeIds.size() << std::endl;
        std::copy( nodeIds.begin(), nodeIds.end(), std::ostream_iterator<unsigned>( sideFile, "\n " ) );

        sideFile.close();

    }

    return;
}


//------------------------------------------------------------------------------
void getConstraints( const PhysicalNamesMap & physicalNames,
                     MonitoredFacesSet & monitoredFaces,
                     MonitoredNodesSet & monitoredNodes,
                     BCdirectionsMap & directions,
                     BCvaluesMap & values )
// get constraints directions, values  and to be extracted entities from  
{
    PhysicalNamesMap::const_iterator pIter = physicalNames.begin();
    PhysicalNamesMap::const_iterator pEnd  = physicalNames.end();

    std::string dirNames( "xyz" );

    for( ; pIter != pEnd; ++ pIter ) {
        
        std::cout << "Give boundary conditions for surface \"" 
                  << pIter -> second << "\" " << std::endl
                  << "Choose between {x=0|y=1|z=2|all=3|none=4|extractEdges=5|extractNodes=6}" 
                  << std ::endl
                  << "Enter number: ";
        unsigned number;
        std::cin >> number;

        std::vector<int long> dir;
        if ( 0 <= number && number <=2 ) {         
            dir.push_back( number );
        }  else if ( number == 3 ) { 
            dir.push_back( 0 );
            dir.push_back( 1 );
            dir.push_back( 2 );
        }
        else if ( number == 5 ) dir.push_back( -1 );
        else if ( number == 6 ) dir.push_back( -2 );

        if ( dir.empty() ) continue;

        // ask for values
        if ( dir[0] > -1 ) {
            std::vector<double> val;
            std::vector<int long>::iterator dIter = dir.begin();
            std::vector<int long>::iterator dEnd  = dir.end();

            for ( ; dIter != dEnd; ++ dIter ) {
                std::cout << "Give value for the " << dirNames[ *dIter ] 
                          << "-direction: ";
                double value;
                std::cin >> value;
                val.push_back( value );
            }
            
            directions[ pIter -> first ] = dir;
            values[ pIter -> first ] = val;
        }

        else if ( dir[0] == -1 ) {
            monitoredFaces.insert( pIter -> first );
        }
        else if ( dir[0] == -2 ) {
            monitoredNodes.insert( pIter -> first );
        }
    }

    return;
}


//------------------------------------------------------------------------------
std::ostream & writeSMF( std::ostream & smfFile,
                         const NodesVec & nodes, 
                         const ElementsVec & elements,
                         unsigned numNodesPerElem,
                         unsigned numNodesPerFace )
{
    FTL_VERIFY( !nodes.empty() );
    FTL_VERIFY( !elements.empty() );
    
    // write header
    corlib::SmfHead smfHead;
    
    GMSHelemType gmshElemType = gmshElement( numNodesPerElem, numNodesPerFace );
    
    smfHead.setElementShape( boost::get<2>(gmshElemType) );
    smfHead.setElementNumPoints( numNodesPerElem );
    smfHead.write( smfFile );

    // write sum of nodes and elements
    smfFile << nodes.size() << "  " << elements.size() << std::endl;
    NodesVec::const_iterator nIter = nodes.begin();
    NodesVec::const_iterator nEnd  = nodes.end();

    for( ; nIter !=  nEnd; ++ nIter ) {
        smfFile << std::fixed << std::setprecision(10) << (*nIter)[0] << "  ";
        smfFile << std::fixed << std::setprecision(10) << (*nIter)[1] << "  ";
        smfFile << std::fixed << std::setprecision(10) << (*nIter)[2] << std::endl;

//        smfFile << (*nIter)[0] << "  " << (*nIter)[1] << "  " << (*nIter)[2] << std::endl;
    }
    
    ElementsVec::const_iterator eIter = elements.begin();
    ElementsVec::const_iterator eEnd  = elements.end();
    for( ; eIter != eEnd; ++ eIter ) {
        for ( unsigned v = 0; v < numNodesPerElem; ++ v ) 
            smfFile << (*eIter)[v] << "  ";
        smfFile << std::endl;
    }
    
    return smfFile;
}


//------------------------------------------------------------------------------
GMSHelemType gmshElement( unsigned numNodesPerElem, 
                          unsigned numNodesPerFace ) 
//  types of element and corresponding face according to GMSH manual
{
    // combination 2-noded line and 1-noded point
    if (numNodesPerElem == 2 and  numNodesPerFace == 1) {
        return boost::make_tuple(1, 15, corlib::LINE);
    }

    // combination 3-noded line and 1-noded point
    if (numNodesPerElem == 3 and numNodesPerFace == 1) {
        return boost::make_tuple(8, 15, corlib::LINE);
    }

    // combination 3-noded triangle and 2-noded lines
    if (numNodesPerElem == 3 and  numNodesPerFace == 2) {
        return boost::make_tuple(2, 1, corlib::TRIANGLE);
    } 

    // combination 6-noded triangle and 3-noded line
    else if (numNodesPerElem == 6 and  numNodesPerFace == 3) { 
        return boost::make_tuple(9, 8, corlib::TRIANGLE);
    } 

    // combination of 4-noded quadrilateral and 2-noded lines
    else if (numNodesPerElem == 4 and  numNodesPerFace == 2) { 
        return boost::make_tuple(3, 1, corlib::QUADRILATERAL);
    } 

    // combination 10-noded tetrahedron and 6-noded triangle
    else if (numNodesPerElem == 10 and  numNodesPerFace == 6) { 
        return boost::make_tuple(11, 9, corlib::TETRAHEDRON);
    } 
    
    // combination 4-noded tetrahedron and 3-noded triangle
    else if (numNodesPerElem == 4 and  numNodesPerFace == 3) { 
        return boost::make_tuple(4, 2, corlib::TETRAHEDRON);
    } 

    // combination 9-noded quadrangle with 3-noded lines
    else if (numNodesPerElem == 9 and  numNodesPerFace == 3) { 
        return boost::make_tuple(10, 8, corlib::QUADRILATERAL);
    } 

    // combination 8-noded quadrangle with 3-noded lines
    else if (numNodesPerElem == 8 and  numNodesPerFace == 3) { 
        return boost::make_tuple(16, 8, corlib::QUADRILATERAL);        
    } 
    
    // combination 27-noded hexahedron and 9-noded face
    else if (numNodesPerElem == 27 and  numNodesPerFace == 9) { 
        return boost::make_tuple(12, 10, corlib::HEXAHEDRON);
    } 
    
    // combination 20-noded hexahedron and 8-noded face
    else if (numNodesPerElem == 20 and  numNodesPerFace == 8) { 
        return boost::make_tuple(17, 16, corlib::HEXAHEDRON);
    } 
    
    // combination 8-noded hexahedron and 4-noded quadrilateral
    else if (numNodesPerElem == 8 and  numNodesPerFace == 4) { 
        return boost::make_tuple(5, 3, corlib::HEXAHEDRON);
    } 

    // 2-noded line 
    else if (numNodesPerElem == 2 and  numNodesPerFace == 2) { 
        return boost::make_tuple(1, 1, corlib::LINE);
    } else {
        FTL_VERIFY_DESCRIPTIVE(false, "Element version not implemented.");
    }
    
    return GMSHelemType(); // to calm down the compiler
}

//------------------------------------------------------------------------------
std::istream & readGMSH( std::istream & inp,
                         PhysicalNamesMap & physicalNames,
                         NodesVec & nodes,
                         ElementsVec & elements,
                         FacesMMap & faces,
                         unsigned numNodesPerElem,
                         unsigned numNodesPerFace )
{
    std::string tagName;
    
    // read tag 'MeshFormat'
    std::string meshFormat;

    inp >> tagName;
    inp.ignore( MAXLINE , '\n' );
    if ( !tagName.compare( 1, 10, "MeshFormat" ) ) { 
        readMeshFormat( inp, meshFormat );
    }
    else 
        FTL_VERIFY_DESCRIPTIVE( false, "Cannot find tag MeshFormat" );

     // read tag 'PhysicalNames
    inp >> tagName;
    inp.ignore( MAXLINE , '\n' );
    bool physical = false;
    if ( !tagName.compare( 1, 13, "PhysicalNames") ) { 
        readPhysicalNames( inp, physicalNames, meshFormat );
        physical = true;
    }
         
    if ( physical ) {
        // skip tag 'EndPhysicalNames'
        inp >> tagName;
        inp.ignore( MAXLINE , '\n' );
    }

    // in 4.* gmsh file, physical entities are obtained through the tag 'Entities'
    EntityMap entityMap;
    if ( meshFormat == "4.1" ){
        readEntitiesVersion4( inp, entityMap );
        inp >> tagName;
        inp.ignore( MAXLINE , '\n' );

        // skip tag 'EndEntities'
        inp >> tagName;
    }

    // read tag 'Nodes'
    if ( !tagName.compare( 1,  5, "Nodes"        ) ) { 
        readNodes( inp, nodes, meshFormat );
    }
    else 
        FTL_VERIFY_DESCRIPTIVE( false, "Cannot find tag Nodes" );

    // skip tag 'EndNodes'
    inp >> tagName;
    inp.ignore( MAXLINE , '\n' );

    // read tag 'Elements'
    if ( !tagName.compare( 1,  8, "Elements"     ) ) { 
        if ( meshFormat == "2" or meshFormat == "2.1" or meshFormat == "2.2")
            readElementsVersion2( inp, elements, faces, numNodesPerElem, numNodesPerFace );
        else
            readElementsVersion4( inp, entityMap, elements, faces, numNodesPerElem, numNodesPerFace );
    }
    else 
        FTL_VERIFY_DESCRIPTIVE( false, "Cannot find tag Elements" );


    //--------------------------------------------------------------------------
    // check for dangling nodes (orphans)
    std::vector<unsigned> orphans;
    findOrphans( nodes.size(), elements, orphans );
    // in case of orphans found
    if ( orphans.size() ) {
        std::cout << "Found these orphans: " << std::endl;
        std::copy( orphans.begin(), orphans.end(), std::ostream_iterator<unsigned>( std::cout, "\n" ) );

        // remove the orphans from nodes and renumber elements
        removeOrphans( orphans, nodes, elements );
        orphans.clear();
        // check again
        findOrphans( nodes.size(), elements, orphans );
        // if orphan persists --> exit with error
        if ( orphans.size() ) {
            std::cout << "Could not remove those orphans: " << std::endl;
            std::copy( orphans.begin(), orphans.end(), std::ostream_iterator<unsigned>( std::cout, "\n" ) );
            FTL_VERIFY( false );
        }
    }

    return inp;
}


//------------------------------------------------------------------------------
std::istream & readMeshFormat( std::istream & inp, std::string & meshFormat )
{
    double dummy;
    inp >> meshFormat >> dummy >> dummy;
    inp.ignore( MAXLINE , '\n' );
    std::string endtag;
    inp >> endtag;
    inp.ignore( MAXLINE , '\n' );
    return inp;
}


//------------------------------------------------------------------------------
std::istream & readPhysicalNames( std::istream & inp,
                                  PhysicalNamesMap & physicalNames,
                                  const std::string & meshFormat )
// read physical names and corresponding integer tags from gmsh file
{
    unsigned numPN;
    inp >> numPN;
    inp.ignore( MAXLINE , '\n' );
    std::cout << numPN << std::endl;
    
    for ( unsigned i = 0; i < numPN; ++ i ) {
        unsigned num;
        unsigned dim;
        
        std::string name;
        if ( meshFormat == "2" )
            inp >> num >> name;
        else if ( meshFormat == "2.1" or meshFormat == "2.2" or meshFormat == "4.1" )
            inp >> dim >> num >> name;
        else
            FTL_VERIFY_DESCRIPTIVE( 
                false, "Provided Gmsh mesh file version unknown: "
                "Please check $PhysicalNames section and look "
                "at which column the physcal name strings appear." );

        std::cout << name << std::endl;
        inp.ignore( MAXLINE , '\n' );
        const unsigned size = name.length();
        physicalNames[ num ] = name.substr( 1, size-2);
    }
    
    std::string endtag;
    inp >> endtag;
    inp.ignore( MAXLINE , '\n' );

    return inp;
}

//------------------------------------------------------------------------------
/** Read entities between the tags "$Entities" and "$EndEntities" in 4.1 gmsh files.
 *  The entity section is formated as
 *  nodeNum edgeNum surfaceNum VolumeNum
 *  node1Index xMin yMin zMin phyeicalNum physicalEntity1 physicalEntity2 ...
 *  node2Index xMin yMin zMin phyeicalNum physicalEntity1 physicalEntity2 ...
 *  ...
 *  edge1Index xMin yMin zMin xMax yMax zMax phyeicalNum physicalEntity1 physicalEntity2 ...
 *  edge2Index xMin yMin zMin xMax yMax zMax phyeicalNum physicalEntity1 physicalEntity2 ...
 *  ...
 *  surface1Index xMin yMin zMin xMax yMax zMax phyeicalNum physicalEntity1 physicalEntity2 ...
 *  surface2Index xMin yMin zMin xMax yMax zMax phyeicalNum physicalEntity1 physicalEntity2 ...
 *  ...
 *  volume1Index xMin yMin zMin xMax yMax zMax phyeicalNum physicalEntity1 physicalEntity2 ...
 *  volume2Index xMin yMin zMin xMax yMax zMax phyeicalNum physicalEntity1 physicalEntity2 ...
 *  ...
 */
std::istream & readEntitiesVersion4( std::istream & inp, EntityMap & entityMap )
// read between the tags "$Entities" and "$EndEntities" in 4.1 gmsh files.
{
    unsigned nodeNum, edgeNum, surfaceNum, volumeNum;
    // read number of different entities
    inp >> nodeNum >> edgeNum >> surfaceNum >> volumeNum;

    unsigned index, tagNum, dim, physicalNum, dummy;

    // loop over node entities for recording the entity with physical names
    for ( unsigned i = 0; i < nodeNum; i ++ ) {

        inp >> index >> dummy >> dummy >> dummy >> physicalNum;

        // there is no physical entity associated with this node
        if ( physicalNum == 0 ) {
            inp.ignore( MAXLINE , '\n' );
        }
        // if there is one physical entity associated with this node, save it
        else if ( physicalNum == 1 ) {
            inp >> tagNum;
            inp.ignore( MAXLINE , '\n' );
            dim = 0;
            entityMap[std::make_pair(dim, index)] = tagNum;
        }
        // latest version does not support multiple physical entities
        else {
            FTL_VERIFY_DESCRIPTIVE(
                false, "Latest GMSH version does not support multiple physical entities. "
                "Check NODE entities." );
        }
    }

    // loop over edge entities for recording the entity with physical names
    for ( unsigned i = 0; i < edgeNum; i ++ ) {

        inp >> index >> dummy >> dummy >> dummy >> dummy >> dummy >> dummy >> physicalNum;

        // there is no physical entity associated with this edge
        if ( physicalNum == 0 ) {
            inp.ignore( MAXLINE , '\n' );
        }
        // if there is one physical entity associated with this edge, save it
        else if ( physicalNum == 1 ) {
            inp >> tagNum;
            inp.ignore( MAXLINE , '\n' );
            dim = 1;
            entityMap[std::make_pair(dim, index)] = tagNum;
        }
        // latest version does not support multiple physical entities
        else {
            FTL_VERIFY_DESCRIPTIVE(
                false, "Latest GMSH version does not support multiple physical entities. "
                "Check EDGE entities." );
        }
    }

    // loop over surface entities for recording the entity with physical names
    for ( unsigned i = 0; i < surfaceNum; i ++ ) {

        inp >> index >> dummy >> dummy >> dummy >> dummy >> dummy >> dummy >> physicalNum;

        // there is no physical entity associated with this surface
        if ( physicalNum == 0 ) {
            inp.ignore( MAXLINE , '\n' );
        }
        // if there is one physical entity associated with this surface, save it
        else if ( physicalNum == 1 ) {
            inp >> tagNum;
            inp.ignore( MAXLINE , '\n' );
            dim = 2;
            entityMap[std::make_pair(dim, index)] = tagNum;
        }
        // latest version does not support multiple physical entities
        else {
            FTL_VERIFY_DESCRIPTIVE(
                false, "Latest GMSH version does not support multiple physical entities. "
                "Check SURFACE entities." );
        }
    }

    // loop over volume entities for recording the entity with physical names
    for ( unsigned i = 0; i < volumeNum; i ++ ) {

        inp >> index >> dummy >> dummy >> dummy >> dummy >> dummy >> dummy >> physicalNum;

        // there is no physical entity associated with this volume
        if ( physicalNum == 0 ) {
            inp.ignore( MAXLINE , '\n' );
        }
        // if there is one physical entity associated with this volume, save it
        else if ( physicalNum == 1 ) {
            inp >> tagNum;
            inp.ignore( MAXLINE , '\n' );
            dim = 3;
            entityMap[std::make_pair(dim, index)] = tagNum;
        }
        // latest version does not support multiple physical entities
        else {
            FTL_VERIFY_DESCRIPTIVE(
                false, "Latest GMSH version does not support multiple physical entities. "
                "Check VOLUME entities." );
        }
    }

    return inp;
}


//------------------------------------------------------------------------------
/** Read nodal coordinates from gmsh file
 *  In 2.* version gmsh file, section 'Node' is formatted as
 *  nodeNum
 *  node1X node1Y node1Z
 *  node2X node2Y node2Z
 *  ...
 *
 *  In 4.* version gmsh file, section 'Node' is formatted as
 *  blockNum(according to different physical entities) nodeNum nodeIndexMin nodeIndexMax
 *  dim physicalTagOfBlock1(index of the tagged entity in .geo file) parametricOrNot nodeNumInBlock1
 *  nodeIndex1InBlock1
 *  nodeIndex2InBlock1
 *  ...
 *  node1InBlock1X node1InBlock1Y node1InBlock1Z
 *  node1InBlock2X node1InBlock2Y node1InBlock2Z
 *  ...
 *  dim physicalTagOfBlock2(index of the tagged entity in .geo file) parametricOrNot nodeNumInBlock2
 *  nodeIndex1InBlock2
 *  nodeIndex2InBlock2
 *  ...
 *  node1InBlock1X node1InBlock1Y node1InBlock1Z
 *  node1InBlock2X node1InBlock2Y node1InBlock2Z
 *  ...
 */
std::istream & readNodes( std::istream & inp, 
                          NodesVec & nodes,
                          const std::string meshFormat = "2.2")
{
    if( meshFormat == "2" or meshFormat == "2.1" or meshFormat == "2.2" ) {

        unsigned numOfNodes;
        inp >> numOfNodes;
        inp.ignore( MAXLINE , '\n' );

        for ( unsigned i = 0; i < numOfNodes; i ++ ) {
            unsigned ind;
            corlib::eigenX::VectorSd<3> coord;
            inp >> ind >> coord[0] >> coord[1] >> coord[2];
            inp.ignore( MAXLINE , '\n' );
            nodes.push_back( coord );
        }
    }
    else if( meshFormat == "4.1" ) {
        unsigned numOfBlock, dummy;
        inp >> numOfBlock >> dummy >> dummy >> dummy;
        inp.ignore( MAXLINE , '\n' );

        // first loop over entity blocks
        for ( unsigned i = 0; i < numOfBlock; i ++ ) {

            unsigned numOfNodesPerBlock;
            inp >> dummy >> dummy >> dummy >> numOfNodesPerBlock;
            inp.ignore( MAXLINE , '\n' );

            // second loop over nodes indices in the block, which we do not care
            for ( unsigned j = 0; j < numOfNodesPerBlock; j ++ ) {
                inp >> dummy;
                inp.ignore( MAXLINE , '\n' );
            }

            // third loop over nodes coordinates in the block
            for ( unsigned j = 0; j < numOfNodesPerBlock; j ++ ) {
                corlib::eigenX::VectorSd<3> coord;
                inp >> coord[0] >> coord[1] >> coord[2];
                inp.ignore( MAXLINE , '\n' );
                nodes.push_back( coord );
            }
        }
    }
    else
        FTL_VERIFY_DESCRIPTIVE( false, "Provided Gmsh mesh file version unknown." );

   std::string endtag;
   inp >> endtag;
   inp.ignore( MAXLINE , '\n' );
   return inp;
}

//------------------------------------------------------------------------------
std::vector<unsigned> permuteHexaMidsides( const std::vector<unsigned> & element )
// perform permutation of midside nodes generated by GMSH to openFTL 
{
    // initialize new element by element
    std::vector<unsigned> permElement ( element );
    // permute midsides
    permElement[  9 ] = element[ 11 ];
    permElement[ 10 ] = element[ 13 ];
    permElement[ 11 ] = element[  9 ];
    permElement[ 12 ] = element[ 16 ];
    permElement[ 13 ] = element[ 18 ];
    permElement[ 14 ] = element[ 19 ];
    permElement[ 15 ] = element[ 17 ];
    permElement[ 16 ] = element[ 10 ];
    permElement[ 17 ] = element[ 12 ];
    permElement[ 18 ] = element[ 14 ];
    permElement[ 19 ] = element[ 15 ];

    return permElement;
}

//------------------------------------------------------------------------------
std::vector<unsigned> permuteHexaFaces( const std::vector<unsigned> & element )
// perform permutation of nodes at centers of faces generated by GMSH to openFTL 
{
    // initialize new element by element
    std::vector<unsigned> permElement ( element );
    // permute faces
    permElement[ 20 ] = element[ 20 ];
    permElement[ 21 ] = element[ 25 ];
    permElement[ 22 ] = element[ 21 ];
    permElement[ 23 ] = element[ 23 ];
    permElement[ 24 ] = element[ 24 ];
    permElement[ 25 ] = element[ 22 ];
    permElement[ 26 ] = element[ 26 ];

    return permElement;
}

//------------------------------------------------------------------------------
/** Read face and element connectivity from the 2.* gmsh file
 *  Faces and elements in the gmsh file are identified with numNodesPerElem and
 *  numNodesPerFace
 */
std::istream & readElementsVersion2( std::istream & inp,
                                     ElementsVec & elements,
                                     FacesMMap & faces,
                                     unsigned numNodesPerElem,
                                     unsigned numNodesPerFace)
{
    unsigned numOfElements;
    inp >> numOfElements;
    inp.ignore( MAXLINE , '\n' );
    
    GMSHelemType gmshType = gmshElement( numNodesPerElem, numNodesPerFace );
    
    const unsigned eType = boost::get<0>( gmshType );
    const unsigned fType = boost::get<1>( gmshType );
    const unsigned pType = 15;  // type of 1-node point

    for ( unsigned e = 0; e < numOfElements; ++ e  ) {

        unsigned index, type, numTags, tagNum, dummy;
        inp >> index >> type >> numTags >> tagNum;
        for ( unsigned t = 1; t < numTags; ++ t ) inp >> dummy;
        
        if ( type == eType ) {
            std::vector<unsigned> element( numNodesPerElem );
            for ( unsigned v = 0; v < numNodesPerElem; ++ v ) {
                inp >> element[v]; element[v] -= 1;
            }
            
            // adapt gmsh numbering to openFTL for ten-noded tet
            if ( numNodesPerElem == 10 ) {
                std::swap(element[9], element[8]);
            }

            // adapt gmsh numbering to openFTL for twenty-noded hex
            if ( numNodesPerElem == 20 and numNodesPerFace == 8 ) {
                element = permuteHexaMidsides ( element );
            }

            // adapt gmsh numbering to openFTL for twenty seven-noded hex
            if ( numNodesPerElem == 27 and numNodesPerFace == 9 ) {
                element = permuteHexaMidsides ( element );
                element = permuteHexaFaces ( element );
            }

            elements.push_back( element );
        }
        else if ( type == fType ) {
            std::vector<unsigned> face( numNodesPerFace );
            for ( unsigned v = 0; v < numNodesPerFace; ++ v ) {
                inp >> face[v]; face[v] -= 1;
            }
            faces.insert( std::make_pair( tagNum, face ) );
        }
        else if ( type == pType ) {
            std::vector<unsigned> pointElement( 1 );
            inp >> pointElement[0];  pointElement[0] -= 1;
            faces.insert( std::make_pair( tagNum, pointElement ) );
        }
        
        inp.ignore( MAXLINE , '\n' );
    }
    
    if ( elements.size() == 0 ) {
        std::cerr << "(EE) Could not find any element of type " << eType << std::endl;
        FTL_VERIFY( false );
    }
    if ( faces.size() == 0 )
        std::cerr << " (WW) Could not find any face of type " << fType << std::endl;

    return inp;
}

//------------------------------------------------------------------------------
/** Read face and element connectivity from the 4.* gmsh file.
 *  Faces and elements in the gmsh file are identified with numNodesPerElem and
 *  numNodesPerFace.
 *  Section 'Element' is formatted as
 *  blockNum(according to different physical entities) eleNum eleIndexMin eleIndexMax
 *  dim physicalTagOfBlock1(index of the tagged entity in .geo file) eleType nodeNumInBlock1
 *  ele1InBlock1Index node1InEle1 node2InEle1 ...
 *  ele2InBlock1Index node1InEle2 node2InEle2 ...
 *  ...
 *  dim physicalTagOfBlock2 eleType nodeNumInBlock2
 *  ele1InBlock2Index node1InEle1 node2InEle1 ...
 *  ele2InBlock2Index node1InEle2 node2InEle2 ...
 *  ...
 */
std::istream & readElementsVersion4( std::istream & inp,
                                     const EntityMap & entityMap,
                                     ElementsVec & elements,
                                     FacesMMap & faces,
                                     unsigned numNodesPerElem,
                                     unsigned numNodesPerFace )
{
    unsigned dim, numOfBlock, type, tagNum, dummy;
    inp >> numOfBlock >> dummy >> dummy >> dummy;
    inp.ignore( MAXLINE , '\n' );

    GMSHelemType gmshType = gmshElement( numNodesPerElem, numNodesPerFace );
       
    const unsigned eType = boost::get<0>( gmshType );
    const unsigned fType = boost::get<1>( gmshType );
    const unsigned pType = 15;  // type of 1-node point

    // first loop over entity blocks
    for ( unsigned i = 0; i < numOfBlock; i ++ ){

        unsigned numOfElementsPerBlock, physicalIndex;
        inp >> dim >> tagNum >> type >> numOfElementsPerBlock;
        inp.ignore( MAXLINE , '\n' );

        // search the physical name index using dim and tagNum
        std::pair<unsigned, unsigned> entityMapKey = std::make_pair( dim, tagNum );
        EntityMap::const_iterator pos = entityMap.find( entityMapKey );
        if ( pos != entityMap.end() )
        // there is a physical entity
            physicalIndex = pos->second;
        else
            physicalIndex = 0;

        // second loop over elements in the block
        for ( unsigned j = 0; j < numOfElementsPerBlock; j ++ ) {
            inp >> dummy;
            if ( type == eType ) {
                std::vector<unsigned> element( numNodesPerElem );
                for ( unsigned v = 0; v < numNodesPerElem; ++ v ) {
                    inp >> element[v]; element[v] -= 1;
                }

                // adapt gmsh numbering to openFTL for ten-noded tet
                if ( numNodesPerElem == 10 ) {
                    std::swap(element[9], element[8]);
                }

                // adapt gmsh numbering to openFTL for twenty-noded hex
                if ( numNodesPerElem == 20 and numNodesPerFace == 8 ) {
                    element = permuteHexaMidsides ( element );
                }

                // adapt gmsh numbering to openFTL for twenty seven-noded hex
                if ( numNodesPerElem == 27 and numNodesPerFace == 9 ) {
                    element = permuteHexaMidsides ( element );
                    element = permuteHexaFaces ( element );
                }

                elements.push_back( element );
            }
            else if ( type == fType ) {
                std::vector<unsigned> face( numNodesPerFace );
                for ( unsigned v = 0; v < numNodesPerFace; ++ v ) {
                    inp >> face[v]; face[v] -= 1;
                }
                faces.insert( std::make_pair( physicalIndex, face ) );
            }
            else if ( type == pType ) {
                std::vector<unsigned> pointElement( 1 );
                inp >> pointElement[0];  pointElement[0] -= 1;
                faces.insert( std::make_pair( physicalIndex, pointElement ) );
            }
        }
    }

    if ( elements.size() == 0 ) {
        std::cerr << "(EE) Could not find any element of type " << eType << std::endl;
        FTL_VERIFY( false );
    }
    if ( faces.size() == 0 )
        std::cerr << " (WW) Could not find any face of type " << fType << std::endl;

    return inp;
}

//------------------------------------------------------------------------------
/** Find orphans in the nodes. An orphan is a node which is not connected to
 *  any element. Such nodes will not carry dofs ans ought to be removed.
 *  \param[in]  numNodes  size of node vector
 *  \param[in]  elements  vector of element connectivities
 *  \param[out] orphans   vector of node numbers not connected to an element
 */
void findOrphans( const unsigned numNodes,
                  const ElementsVec & elements,
                  std::vector<unsigned> & orphans )
{
    // markers
    std::vector<bool> isNodeOfElement( numNodes, false );

    // check elements / faces for references to nodes
    ElementsVec::const_iterator eIter = elements.begin();
    ElementsVec::const_iterator eEnd  = elements.end();
    unsigned numNodesPerElem = ( elements[0] ).size();
    // go through elements
    for ( ; eIter != eEnd; eIter++ ) {
        // go through connected node numbers
        for ( unsigned v = 0; v < numNodesPerElem; v++ ) {
            // mark connected node numbers
            isNodeOfElement[ (*eIter)[v] ] = true;
        }
    }

    // collect indices of orphans
    for ( unsigned n = 0; n < numNodes; n ++ ) {
        if ( not isNodeOfElement[n] ) orphans.push_back(n);
    }

    return;
}

//------------------------------------------------------------------------------
/** Remove orphans from node vector and renumber connectivity accordingly
 *  \param[in]      orphans  vector of node indices which are orphans
 *  \param[in,out]  nodes    vector of node coordinates 
 *  \param[in,out]  elements vector of element connectivities
 */
void removeOrphans( std::vector<unsigned> & orphans,
                    NodesVec & nodes, 
                    ElementsVec & elements ) 
{
    std::cout << "Removing orphans " << std::endl;

    const unsigned numNodesPerElement = (elements[0]).size();

    // sort in order to operate from highest to lowest number
    std::sort( orphans.begin(), orphans.end() );

    // go from highest to lowest orphan number
    for ( int i = orphans.size() -1; i >= 0; i -- ) {

        const unsigned orphan = orphans[i];

        //----------------------------------------------------------------------
        // renumber element connectivity:
        // every node index larger than the orphan's number has
        // to be decreased by one
        for (unsigned e = 0; e < elements.size(); e++ ) {
            for (unsigned v = 0; v < numNodesPerElement; v ++ ) {
                if ( elements[e][v] > orphan )
                    elements[e][v] --;
            }
        }

        //----------------------------------------------------------------------
        // remove orphan from vector of nodes (important: have orphan numbers
        //  sorted first)
        nodes.erase( nodes.begin() + orphan );
    }

    return;
}
                    
