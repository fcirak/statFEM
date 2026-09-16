#include <vtkSmartPointer.h>
#include <vtkXMLUnstructuredGridReader.h>
#include <vtkUnstructuredGrid.h>
#include <vtkXMLPolyDataWriter.h>
#include <vtkGeometryFilter.h>
#include <vtkPointData.h>
#include <vtkXMLUnstructuredGridWriter.h>

#include <boost/lexical_cast.hpp>
//------------------------------------------------------------------------------
// Convert a vtu to vtp file
//------------------------------------------------------------------------------
int main( int argc, char *argv[] )
{
    //parse command line arguments
    if( argc != 4 ) {
        std::cerr << "Usage: " << argv[0]
                  << " filename.vtu  delta1  delta2" << std::endl;
        return EXIT_FAILURE;
    }

    const double delta1 = boost::lexical_cast<double>( argv[2] );
    const double delta2 = boost::lexical_cast<double>( argv[3] );

    // input and output file names
    vtkstd::string iFilename = argv[1];
    vtkstd::string oFilename = iFilename.substr( 0, iFilename.find( ".vtu" ) ) + "_.vtu";
 
    //read all the data from the file
    vtkSmartPointer<vtkXMLUnstructuredGridReader> reader = vtkSmartPointer<vtkXMLUnstructuredGridReader>::New();
    reader->SetFileName( iFilename.c_str() );
    reader->Update();
 
    // get the pointdata
    vtkSmartPointer<vtkPointData> pointData = (reader -> GetOutput()) -> GetPointData();
    const unsigned nArrays = pointData -> GetNumberOfArrays();
    bool found = false;
    vtkstd::string vectorName;
    
    // looking for point data field with 3 components
    unsigned array = 0;
    for ( ; array < nArrays; array ++ ) {

        // found first field with 3 components -> the vector field
        if (pointData -> GetArray( array ) -> GetNumberOfComponents() == 3 ) {
            vectorName = pointData -> GetArrayName( array );
            found = true;
            break;
        }
    }

    vtkstd::cout << "Adding (" << delta1 << ", " << delta2 << ") to "
                 << vectorName << vtkstd::endl;


    // extract field
    vtkSmartPointer<vtkDataArray> field = pointData -> GetArray( array );

    // number of entries in that field
    const unsigned numVectors = field -> GetNumberOfTuples();
    
    // add number to each tuple
    for ( unsigned n = 0; n < numVectors; n ++ ) {
        double * tuple = field -> GetTuple( n );
        tuple[0] += delta1;
        tuple[1] += delta2;
        field -> SetTuple( n, tuple );
    }

    // write XML binary vtu file
    vtkSmartPointer<vtkXMLUnstructuredGridWriter> vtuWriter = 
        vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
    vtuWriter -> SetFileName( oFilename.c_str() );
    vtuWriter -> SetInputConnection( reader->GetOutputPort() );
    vtuWriter -> SetDataModeToAppended();
    vtuWriter -> EncodeAppendedDataOff();
    vtuWriter -> Write();


    return 0;
}
