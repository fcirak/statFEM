#include <vtkSmartPointer.h>
#include <vtkXMLPolyDataReader.h>
#include <vtkPolyData.h>
#include <vtkUnstructuredGrid.h>
#include <vtkPoints.h>
#include <vtkPointData.h>
#include <vtkIdList.h>
//------------------------------------------------------------------------------
#include <string>
#include <iostream>
#include <vector>
#include <cassert>
#include <sstream>
//------------------------------------------------------------------------------
#include "TECIO.h"


// return the name of the element type for TECPLOT
std::string elementName( const unsigned vtkNum )
{
    switch( vtkNum ) {
    case 22: // TRI6
    case  5: // TRI3
        return "FETRIANGLE"; break;
    case 23: // QUAD8/9
    case 8:  // PIXEL
    case 9:  // QUAD4
        return "FEQUADRILATERAL"; break;
    case 11: // VOXEL
    case 12: // HEX8
        return "FEBRICK"; break;
    default: return "UNDEFINED";
    }
    return "UNDEFINED";
}

// number of nodes  per element
unsigned numNodesPerElement( const unsigned vtkNum )
{
    switch( vtkNum ) {
    case 22:
    case  5: return 3; break;
    case 23:
    case 8: 
    case 9: return 4; break;
    case 11: 
    case 12: return 8; break;
    default: return 0;
    }
    return 0;
}

// number of nodes  per element
unsigned numNodesPerCell( const unsigned vtkNum )
{
    switch( vtkNum ) {
    case 22: return 6; break;
    case  5: return 3; break;
    case 23: return 8; break;
    case 8: 
    case 9:  return 4; break;
    case 11: 
    case 12: return 8; break;
    default: return 0;
    }
    return 0;
}


// number of sub-elements per element (quadratic elements are partitioned)
unsigned numElementsPerCell( const unsigned vtkNum )
{
    switch( vtkNum ) {
    case 22 : 
    case 23 : return 4; break;
    default : return 1; break;
    }
    return 0;
}

// return TECPLOT element type number
INTEGER4 zoneType( const unsigned vtkNum )
{
    switch( vtkNum ) {
    case 22:
    case  5: return 2; break;
    case 23:
    case 8: 
    case 9: return 3; break;
    case 11:
    case 12: return 5; break;
    }
    return -1;
}

// swap some nodes for PIXEL elements
vtkIdList * reverse( const unsigned vtkNum, vtkIdList * original )
{
    if ( vtkNum == 8 ) { // PIXEL element -> swap last two nodes
        const unsigned tmp = original -> GetId( 2 );
        original -> SetId( 2, original -> GetId( 3 ) );
        original -> SetId( 3, tmp );
    }
    
    if ( vtkNum == 11 ) { // VOXEL element -> swap some nodes
        const unsigned tmp1 = original -> GetId( 2 );
        original -> SetId( 2, original -> GetId( 3 ) );
        original -> SetId( 3, tmp1 );
        const unsigned tmp2 = original -> GetId( 6 );
        original -> SetId( 6, original -> GetId( 7 ) );
        original -> SetId( 7, tmp2 );
    }

    return original;
}

// lines get too long in BLOCK mode for ascii output
char lineBreak( const unsigned num )
{
    if ( num % 100 == 0 ) return '\n';
    return ' ';
}

// make a name out of first character and component number
std::string makeName( std::string name, unsigned comp )
{
    std::stringstream aux;
    aux << name[0] << name[1] << name[2] << comp ;
    return aux.str();
}

//------------------------------------------------------------------------------
int main(int argc, char *argv[])
{
    //parse command line arguments
    if( argc != 2 ) {
        std::cerr << "Usage: " << argv[0]
                  << " filename.vtp" << std::endl;
        return EXIT_FAILURE;
    }
 
    // file names
    std::string filename       = argv[1];
    std::string basename       = filename.substr( 0, filename.find( ".vtp" ) );
    std::string binaryOut      = basename + ".plt";
 
    //read XML vtp file
    vtkSmartPointer<vtkXMLPolyDataReader> reader =
        vtkSmartPointer<vtkXMLPolyDataReader>::New();
    reader -> SetFileName( filename.c_str() );
    reader -> Update();

    // get access to the unstructured grid
    vtkSmartPointer<vtkPolyData> polyData = reader -> GetOutput();

    // get header data
    const unsigned numPoints   = polyData -> GetNumberOfPoints();
    const unsigned cellType    = polyData -> GetCellType( 0 );
    const unsigned numCells    = polyData -> GetNumberOfCells();
    const unsigned numSubElem  = numElementsPerCell( cellType );
    const unsigned numElements = numCells * numSubElem;

    // point data
    vtkSmartPointer<vtkPointData> pointData =  polyData -> GetPointData();
    const unsigned numPointDataComponents = pointData -> GetNumberOfComponents();
    const unsigned numPointDataArrays     = pointData -> GetNumberOfArrays();

    // access to points
    vtkSmartPointer<vtkPoints> points = polyData -> GetPoints();


    // generate data header string
    std::string dataHeader = "X  Y  Z ";
    for ( unsigned na = 0; na < numPointDataArrays; na ++ ) {
        std::string dataName = ( pointData -> GetArray( na ) ) -> GetName();
        const unsigned numArrayComponents = ( pointData -> GetArray( na ) ) -> GetNumberOfComponents();
        for ( unsigned nc = 0; nc < numArrayComponents; nc ++ ) {
            std::string appendix = makeName( dataName, nc ) + "  ";
            dataHeader.append( appendix );
        }
    }
    std::cout << dataHeader << std::endl;
    
    INTEGER4 fieldFileType = 0;
    INTEGER4 debug = 0;
    INTEGER4 vIsDouble = 0;

    // open and initialize binary file
    INTEGER4 I = TECINI112( (char*)"Generated by vtu2tec",
                            (char*)dataHeader.c_str(), 
                            (char*)binaryOut.c_str(),
                            (char*)".",
                            &fieldFileType,
                            &debug,
                            &vIsDouble );

    // write a zone
    INTEGER4 zT = zoneType( cellType );
    INTEGER4 dummyInt   = 0;
    double   dummyFloat = 0.0;
    INTEGER4 isBlock    = 1;
    INTEGER4 nP = numPoints;
    INTEGER4 nE = numElements;
    I = TECZNE112((char*)"Zone 1",
                  &zT,
                  &nP,
                  &nE,
                  &dummyInt,
                  &dummyInt,
                  &dummyInt,
                  &dummyInt,
                  &dummyFloat,
                  &dummyInt,
                  &dummyInt,
                  &isBlock,
                  &dummyInt,
                  &dummyInt,
                  0,         /* TotalNumFaceNodes */
                  0,         /* NumConnectedBoundaryFaces */
                  0,         /* TotalNumBoundaryConnections */
                  NULL,      /* PassiveVarList */
                  NULL,      /* ValueLocation */
                  NULL,      /* ShareVarFromZone */
                  &dummyInt);

    // write coordinates
    std::vector<double> X, Y, Z;
    for ( unsigned np = 0; np < numPoints; np ++ ) {
        // coordinates
        double xyz[3];
        points -> GetPoint( np, xyz );
        X.push_back( xyz[0] );
        Y.push_back( xyz[1] );
        Z.push_back( xyz[2] );
    }
    INTEGER4 isDouble = 1;
    I   = TECDAT112(&nP, &(X[0]), &isDouble);
    I   = TECDAT112(&nP, &(Y[0]), &isDouble);
    I   = TECDAT112(&nP, &(Z[0]), &isDouble);


    // write point data
    for ( unsigned na = 0; na < numPointDataArrays; na ++ ) {
        vtkSmartPointer<vtkDataArray> array = pointData -> GetArray( na );
        const unsigned numArrayComponents = array -> GetNumberOfComponents();
        for ( unsigned nc = 0; nc < numArrayComponents; nc ++ ) {
            std::vector<double> data;
            for ( unsigned np = 0; np < numPoints; np ++ ) {
                double datum[ numArrayComponents ];
                array -> GetTuple( np, datum );
                data.push_back( datum[nc] );
            }
            I = TECDAT112( &nP, &(data[0]), &isDouble );
        }
    }

    // write connectivity
    const unsigned numNPE = numNodesPerElement( cellType );

    //INTEGER4 connec[numElements][ numNPE ];
    INTEGER4 * connec = new INTEGER4 [ numElements * numNPE ];
    
    for ( unsigned c = 0; c < numCells; c ++ ) {
        const unsigned otherCT = polyData -> GetCellType( c );
        assert( polyData -> GetCellType( c ) == cellType );
        vtkSmartPointer<vtkIdList> idList = vtkSmartPointer<vtkIdList>::New();
        idList -> SetNumberOfIds( numNPE );
        polyData -> GetCellPoints( c, idList );
        idList = reverse( cellType, idList );
        for ( unsigned p = 0; p < numNPE; p ++ ) {
            //connec[ c ][ p ] = idList -> GetId( p ) + 1;
            connec[ c * numNPE + p ] = idList -> GetId( p ) + 1;
        }
    }
    I = TECNOD112((INTEGER4 *)connec);

    // terminate file io
    I = TECEND112();

    delete[] connec;

    return EXIT_SUCCESS;
}

 
