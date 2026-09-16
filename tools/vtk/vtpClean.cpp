#include <vtkSmartPointer.h>
#include <vtkXMLPolyDataReader.h>
#include <vtkPolyData.h>
#include <vtkCleanPolyData.h>
#include <vtkXMLPolyDataWriter.h>

#include <boost/lexical_cast.hpp>

int main( int argc, char * argv[] )
{
    if( argc != 2 ) {

        std::cout << "Usage: " << argv[0]  
                  << "  file.vtp  " << std::endl << std::endl
                  << " where: " << std::endl
                  << "  file.vtp  - the data source (converted by vtu2vtp) " << std::endl;

        return 1;
    }

    const std::string filename = boost::lexical_cast<std::string>( argv[1] );

    // Read XML vtp file
    vtkSmartPointer<vtkXMLPolyDataReader> reader =
        vtkSmartPointer<vtkXMLPolyDataReader>::New();
    reader -> SetFileName( filename.c_str() );
    reader -> Update();

    // Clean the poly data
    vtkSmartPointer<vtkCleanPolyData> cleaner = vtkSmartPointer<vtkCleanPolyData>::New();
    cleaner -> SetInputConnection( reader -> GetOutputPort() );

    // write data back to same file
    vtkSmartPointer<vtkXMLPolyDataWriter> writer = vtkSmartPointer<vtkXMLPolyDataWriter>::New();
    writer->SetFileName( filename.c_str() );
    writer->SetInput( cleaner -> GetOutput() );
    writer -> SetDataModeToAppended();
    writer -> EncodeAppendedDataOff();
    writer->Write();


    return 0;
}
