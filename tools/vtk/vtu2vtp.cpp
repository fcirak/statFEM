#include <vtkSmartPointer.h>
#include <vtkXMLUnstructuredGridReader.h>
#include <vtkUnstructuredGrid.h>
#include <vtkXMLPolyDataWriter.h>
#include <vtkGeometryFilter.h>

//------------------------------------------------------------------------------
// Convert a vtu to vtp file
//------------------------------------------------------------------------------
int main(int argc, char *argv[])
{
    //parse command line arguments
    if( argc != 2 ) {
        std::cerr << "Usage: " << argv[0]
                  << " filename.vtu" << std::endl;
        return EXIT_FAILURE;
    }

    // input and output file names
    vtkstd::string iFilename = argv[1];
    vtkstd::string oFilename = iFilename.substr( 0, iFilename.find( ".vtu" ) ) + ".vtp";
 
    //read all the data from the file
    vtkSmartPointer<vtkXMLUnstructuredGridReader> reader = vtkSmartPointer<vtkXMLUnstructuredGridReader>::New();
    reader->SetFileName( iFilename.c_str() );
    reader->Update();
 
    // extract geometry
    vtkSmartPointer<vtkGeometryFilter> extract = vtkSmartPointer<vtkGeometryFilter>::New();
    extract -> SetInput( reader -> GetOutput() );
    extract -> Update();

    // write VTP file
    vtkSmartPointer<vtkXMLPolyDataWriter> writer = vtkSmartPointer<vtkXMLPolyDataWriter>::New();
    writer -> SetFileName( oFilename.c_str() );
    writer -> SetInput( extract -> GetOutput( ) );
    writer -> SetDataModeToAppended();
    writer -> EncodeAppendedDataOff();
    writer -> Write();

    return 0;
}
