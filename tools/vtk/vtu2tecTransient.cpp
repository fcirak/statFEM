#include <vtkSmartPointer.h>
#include <vtkXMLUnstructuredGridReader.h>
#include <vtkXMLUnstructuredGridWriter.h>
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
#include <boost/lexical_cast.hpp>

#include <corlib/UniqueFilename.hpp>
//------------------------------------------------------------------------------
#include "TECIO.h"


// return the name of the element type for TECPLOT
std::string elementName( const unsigned vtkNum )
{
    switch( vtkNum ) {
    case 22:
    case  5: return "FETRIANGLE"; break;
    case 23:
    case 8:  
    case 9:  return "FEQUADRILATERAL"; break;
    case 11: 
    case 12: return "FEBRICK"; break;
    default: return "UNDEFINED";
    }
    return "UNDEFINED";
}

// number of nodes  per element
unsigned numNodes( const unsigned vtkNum )
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

// swap some nodes for PIXEL/VOXEL elements
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



// make a name out of first character and component number
std::string makeName( std::string name, unsigned comp )
{
    std::stringstream aux;
    aux << name[0] << comp ;
    return aux.str();
}

//------------------------------------------------------------------------------
int main(int argc, char *argv[])
{
    //parse command line arguments
    if( argc != 6 ) {
        std::cerr << "Usage: " << argv[0]
                  << " basename " 
                  << " XXXX (begin number string)" 
                  << " YYYY (end number string)"
                  << " shift "
                  << " deltaT "
                  << std::endl;
        return EXIT_FAILURE;
    }
 
    // file names
    std::string basename       = argv[1];
    std::string binaryOut      = basename + ".plt";

    INTEGER4 I;

    // analyze time stepping
    const unsigned begin  = boost::lexical_cast<unsigned>( argv[2] );
    const unsigned end    = boost::lexical_cast<unsigned>( argv[3] );
    const unsigned inc    = boost::lexical_cast<unsigned>( argv[4] );
    const double   deltaT = boost::lexical_cast<double>(   argv[5] );
    const std::string endNumber = std::string( argv[2] );
    const unsigned numeralWidth = endNumber.size();

    // count the number of variables (at least three for the coordinates)
    unsigned numVars = 3;


    // go through time steps
    for ( unsigned step = begin; step <= end; step += inc ) {

        // generate the filename
        std::string filename = basename + corlib::UniqueFilename::convertNumberToString( step, numeralWidth ) + ".vtu";

        //read XML vtu file
        vtkSmartPointer<vtkXMLUnstructuredGridReader> reader =
            vtkSmartPointer<vtkXMLUnstructuredGridReader>::New();
        reader -> SetFileName( filename.c_str() );
        reader -> Update();

        // get access to the unstructured grid
        vtkSmartPointer<vtkUnstructuredGrid> unstructuredGrid = reader -> GetOutput();

        // get header data
        const unsigned numPoints = unstructuredGrid -> GetNumberOfPoints();
        const unsigned numCells  = unstructuredGrid -> GetNumberOfCells();
        const unsigned cellType  = unstructuredGrid -> GetCellType( 0 );

        // point data
        vtkSmartPointer<vtkPointData> pointData =  unstructuredGrid -> GetPointData();
        const unsigned numPointDataComponents = pointData -> GetNumberOfComponents();
        const unsigned numPointDataArrays     = pointData -> GetNumberOfArrays();

        // access to points
        vtkSmartPointer<vtkPoints> points = unstructuredGrid -> GetPoints();

        // open and initialize binary file ( first step only )
        if ( step == begin ) {
            INTEGER4 fieldFileType = 0;
            INTEGER4 debug = 0;
            INTEGER4 vIsDouble = 0;

            // generate data header string
            std::string dataHeader = "X  Y  Z ";
            for ( unsigned na = 0; na < numPointDataArrays; na ++ ) {
                std::string dataName = ( pointData -> GetArray( na ) ) -> GetName();
                const unsigned numArrayComponents = ( pointData -> GetArray( na ) ) -> GetNumberOfComponents();
                for ( unsigned nc = 0; nc < numArrayComponents; nc ++ ) {
                    std::string appendix = makeName( dataName, nc ) + "  ";
                    dataHeader.append( appendix );
                    numVars ++;
                }
            }

            std::cout << "Header: " << dataHeader << std::endl;

            I = TECINI112( (char*)"Generated by vtu2tecTransient",
                           (char*)dataHeader.c_str(), 
                           (char*)binaryOut.c_str(),
                           (char*)".",
                           &fieldFileType,
                           &debug,
                           &vIsDouble );
        }

        // write a zone
        INTEGER4 zT = zoneType( cellType );
        INTEGER4 dummyInt   = 0;
        INTEGER4 isBlock    = 1;
        INTEGER4 nP = numPoints;
        INTEGER4 nE = numCells;
        double   solTime  = step * deltaT;
        INTEGER4 strandID = 1;

        INTEGER4 shareVarFromZone    = ( step == begin ? 0 : 1 );
        INTEGER4 shareConnecFromZone = ( step == begin ? 0 : 1 );

        INTEGER4 sharedVars[numVars];
        if ( step == begin ) { std::fill( sharedVars, sharedVars + numVars, 0 ); }
        else {
            sharedVars[0] = 1; sharedVars[1] = 1; sharedVars[2] = 1; 
            std::fill( sharedVars + 3, sharedVars + numVars, 0 );
        }

        // generate a name of the zone 
        const std::string zoneTitle =  "Zone" + makeName( " ", step );

        std::cout << "Writing " << zoneTitle << " to file " << std::endl;

        I = TECZNE112( (char *)zoneTitle.c_str(),
                       &zT,
                       &nP,
                       &nE,
                       &dummyInt,
                       &dummyInt,
                       &dummyInt,
                       &dummyInt,
                       &solTime,
                       &strandID,
                       &dummyInt,
                       &isBlock,
                       &dummyInt,
                       &dummyInt,
                       0,         /* TotalNumFaceNodes */
                       0,         /* NumConnectedBoundaryFaces */
                       0,         /* TotalNumBoundaryConnections */
                       NULL,      /* PassiveVarList */
                       NULL,      /* ValueLocation */
                       sharedVars,
                       &shareConnecFromZone);

        INTEGER4 isDouble = 1;

        // write coordinates; only in the first step, since they are shared
        if ( step == begin ) {
            std::vector<double> X, Y, Z;
            for ( unsigned np = 0; np < numPoints; np ++ ) {
                // coordinates
                double xyz[3];
                points -> GetPoint( np, xyz );
                X.push_back( xyz[0] );
                Y.push_back( xyz[1] );
                Z.push_back( xyz[2] );
            }
            I   = TECDAT112(&nP, &(X[0]), &isDouble);
            I   = TECDAT112(&nP, &(Y[0]), &isDouble);
            I   = TECDAT112(&nP, &(Z[0]), &isDouble);
        }


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

        // write connectivity; only in the first step since it is shared
        if ( step == begin ) {
            const unsigned numNPE = numNodes( cellType );
            INTEGER4 connec[numCells][ numNPE ];
            for ( unsigned c = 0; c < numCells; c ++ ) {
                assert( unstructuredGrid -> GetCellType( c ) == cellType );
                vtkSmartPointer<vtkIdList> idList = vtkSmartPointer<vtkIdList>::New();
                idList -> SetNumberOfIds( numNPE );
                unstructuredGrid -> GetCellPoints( c, idList );
                idList = reverse( cellType, idList );
                for ( unsigned p = 0; p < numNPE; p ++ ) {
                    connec[ c ][ p ] = idList -> GetId( p ) + 1;
                }
            }
            I = TECNOD112((INTEGER4 *)connec);
        }

    
    }
    
    // terminate file io
    I = TECEND112();
    
    return EXIT_SUCCESS;
}

 
