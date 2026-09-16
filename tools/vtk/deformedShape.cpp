#include <vtkSmartPointer.h>
#include <vtkXMLUnstructuredGridReader.h>
#include <vtkUnstructuredGrid.h>
#include <vtkPointData.h>
#include <vtkXMLUnstructuredGridWriter.h>

#include <boost/lexical_cast.hpp>
#include <corlib/UniqueFilename.hpp>

int main( int argc, char * argv[] )
{
    if( argc != 2 ) {
        std::cerr << "Usage: " << argv[0]
                  << " filename " 
                  << std::endl;
        return EXIT_FAILURE;
    }

    std::string filename    = boost::lexical_cast<std::string>( argv[1] );
    std::string basename    = filename.substr( 0, filename.find( ".vtu" ) );
    std::string outFileName = basename + ".def";

    //read XML vtu file
    vtkSmartPointer<vtkXMLUnstructuredGridReader> reader =
        vtkSmartPointer<vtkXMLUnstructuredGridReader>::New();
    reader -> SetFileName( filename.c_str() );
    reader -> Update();

    // get access to the unstructured grid
    vtkSmartPointer<vtkUnstructuredGrid> unstructuredGrid = reader -> GetOutput();

    // point data
    vtkSmartPointer<vtkPointData> pointData =  unstructuredGrid -> GetPointData();

    
    const unsigned numPoints = unstructuredGrid -> GetNumberOfPoints();

    // choose datum
    std::cout << "Choose point datum from " << std::endl;
    const unsigned numPointDataArrays = pointData -> GetNumberOfArrays();
    for ( unsigned na = 0; na < numPointDataArrays; na ++ ) {
        const std::string dataName = ( pointData -> GetArray( na ) ) -> GetName();
        std::cout << "  " << na << ": " << dataName << std::endl;
    }
    unsigned pointDatum;
    std::cin >> pointDatum;

    vtkSmartPointer<vtkDataArray> field = pointData -> GetArray( pointDatum );


    std::ofstream out( outFileName.c_str() );

    for ( unsigned n = 0; n < numPoints; n ++ ) {
        double * point = unstructuredGrid -> GetPoint( n );
        double * disp  = field -> GetTuple( n );
        for ( unsigned d = 0; d < 3; d ++ )
            out << point[d] + disp[d] << "  ";
        out << std::endl;
    }

    out.close();

    return 0;
}
