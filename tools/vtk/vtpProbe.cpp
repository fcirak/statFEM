#include <vtkSmartPointer.h>
#include <vtkXMLPolyDataReader.h>
#include <vtkProbeFilter.h>
#include <vtkPlaneSource.h>
#include <vtkPolyData.h>
#include <vtkXMLPolyDataWriter.h>
#include <vtkPointData.h>

#include <vector>
#include <fstream>
#include <boost/lexical_cast.hpp>
#include <corlib/UniqueFilename.hpp>

int main( int argc, char * argv[] )
{
    if( argc != 3 ) {

        std::cout << "Usage: " << argv[0]  
                  << "  file.vtp  grid " << std::endl << std::endl
                  << " where: " << std::endl
                  << "  file.vtp  - the data source (converted by vtu2vtp) " << std::endl
                  << "  grid      - the grid defined by: " << std::endl << std::endl
                  << "     orX  orY  orZ " << std::endl
                  << "     p1X  p1Y  p1Z " << std::endl
                  << "     p2X  p2Y  p2Z " << std::endl
                  << "     res1  res2  "   << std::endl << std::endl
                  << "  which denote the origin (first line), the two "
                  << "end points " << std::endl
                  << "  spanning the grid (second and third line), and " << std::endl
                  << "  the resolutions along these directions " << std::endl;
        return 1;
    }

    const std::string filename = boost::lexical_cast<std::string>( argv[1] );
    const std::string gridFile = boost::lexical_cast<std::string>( argv[2] );

    const std::string basename = filename.substr( 0, filename.find( ".vtp" ) );
   
    // Read XML vtp file
    vtkSmartPointer<vtkXMLPolyDataReader> reader =
        vtkSmartPointer<vtkXMLPolyDataReader>::New();
    reader -> SetFileName( filename.c_str() );
    reader -> Update();

    // read grid file do define a plane
    std::ifstream grid( gridFile.c_str() );
    double orX, orY, orZ, p1X, p1Y, p1Z, p2X, p2Y, p2Z;
    unsigned res1, res2;
    grid >> orX >> orY >> orZ;
    grid >> p1X >> p1Y >> p1Z;
    grid >> p2X >> p2Y >> p2Z;
    grid >> res1 >> res2;
    grid.close();

    // set up a plane source 
    vtkSmartPointer<vtkPlaneSource> plane = vtkSmartPointer<vtkPlaneSource>::New();
    plane -> SetOrigin( orX, orY, orZ );
    plane -> SetPoint1( p1X, p1Y, p1Z );
    plane -> SetPoint2( p2X, p2Y, p2Z );
    plane -> SetXResolution( res1 );
    plane -> SetYResolution( res2 );

    // Set up a probe filter
    vtkSmartPointer<vtkProbeFilter> probeFilter = vtkSmartPointer<vtkProbeFilter>::New();
    probeFilter -> SetInputConnection(  plane -> GetOutputPort() );
    probeFilter -> SetSourceConnection( reader -> GetOutputPort() );
    probeFilter -> Update();

    // get output from probe
    vtkSmartPointer<vtkPolyData> probeData = probeFilter -> GetPolyDataOutput();

    // write VTP file
    vtkSmartPointer<vtkXMLPolyDataWriter> writer = vtkSmartPointer<vtkXMLPolyDataWriter>::New();
    const std::string vtpOutFilename = basename + "_probe.vtp";
    writer -> SetFileName( vtpOutFilename.c_str() );
    writer -> SetInput( probeData );
    writer -> SetDataModeToAppended();
    writer -> EncodeAppendedDataOff();
    writer -> Write();

    //--------------------------------------------------------------------------
    // write ASCII file

    // extract  point data
    vtkSmartPointer<vtkPointData> pointData = probeData -> GetPointData();
    const unsigned nArrays = pointData -> GetNumberOfArrays();
    const std::string asciiFile = basename + ".dat";
    std::ofstream out( asciiFile.c_str() );
    // write header 
    out << "# coordinates(3) ";
    for ( unsigned a =0; a < nArrays; a ++ )
        out << pointData -> GetArrayName( a ) << "("
            << (pointData -> GetArray(a)) -> GetNumberOfComponents()
            << ")  ";
    out << std::endl;

    // write point data
    for ( unsigned p = 0; p < probeData -> GetNumberOfPoints(); p ++ ) {
        // write coordinates
        double * point = probeData -> GetPoint( p );
        for ( unsigned d = 0; d < 3; d ++ )
            out << point[d] << "   ";
        // write data
        for ( unsigned a = 0; a < nArrays; a ++ ) {
            // number of components
            const unsigned nComp = (pointData -> GetArray(a)) -> GetNumberOfComponents(); 
            // write all components
            double * tuple = (pointData -> GetArray(a)) -> GetTuple(p);
            for ( unsigned c = 0; c < nComp; c ++ ) 
                out << tuple[c] << " ";
        }
        out << std::endl;
    }
    out.close();

    return 0;
}
