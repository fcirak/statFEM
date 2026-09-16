#include <vtkSmartPointer.h>
#include <vtkXMLUnstructuredGridReader.h>
#include <vtkUnstructuredGrid.h>
#include <vtkClipDataSet.h>
#include <vtkPolyData.h>
#include <vtkClipPolyData.h>
#include <vtkGeometryFilter.h>
#include <vtkBox.h>
#include <vtkXMLUnstructuredGridWriter.h>
#include <vtkXMLPolyDataWriter.h>

#include <boost/lexical_cast.hpp>
#include <corlib/UniqueFilename.hpp>

//------------------------------------------------------------------------------
int main( int argc, char * argv[] )
{
    if( argc < 8  or argc > 9 ) {
        std::cerr << "Usage: " << argv[0] 
                  << " filename.vtu xmin xmax ymin ymax zmin zmax [2D]"  <<  std::endl
                  << "  with the input \'filename\' and the bounding box. " << std::endl
                  << "  2D is an optional flag (default=0) which enforces " << std::endl
                  << "  a triangulated poly-data output. " << std::endl;
        return EXIT_FAILURE;
    }

    // read input file name
    std::string filename    = boost::lexical_cast<std::string>( argv[1] );
    std::string basename    = filename.substr( 0, filename.find( ".vtu" ) );

    // read bounding box
    const double xmin = boost::lexical_cast<double>( argv[2] );
    const double xmax = boost::lexical_cast<double>( argv[3] );
    const double ymin = boost::lexical_cast<double>( argv[4] );
    const double ymax = boost::lexical_cast<double>( argv[5] );
    const double zmin = boost::lexical_cast<double>( argv[6] );
    const double zmax = boost::lexical_cast<double>( argv[7] );

    // flag for 2D conversion: generate a vtp-surface file which is triangulated
    const bool twoDee = ( argc == 9 ? boost::lexical_cast<bool>( argv[8] ) : false );

    //read XML vtu file
    vtkSmartPointer<vtkXMLUnstructuredGridReader> reader =
        vtkSmartPointer<vtkXMLUnstructuredGridReader>::New();
    reader -> SetFileName( filename.c_str() );
    reader -> Update();

    // get unstructured grid
    vtkSmartPointer<vtkUnstructuredGrid> grid = reader -> GetOutput();

    // set an implicit function
    vtkSmartPointer<vtkBox> box = vtkSmartPointer<vtkBox>::New();
    box -> SetBounds( xmin, xmax, ymin, ymax, zmin, zmax );

    std::string outFileName = basename + ".clipped.vt" + ( twoDee ? "p" : "u" );

    if ( twoDee ) {
        // go via poly data and write a polydata file

        // extract the surface
        vtkSmartPointer<vtkGeometryFilter> surf = vtkSmartPointer<vtkGeometryFilter>::New();
        surf -> SetInput( grid );
        surf -> Update();
        vtkSmartPointer<vtkPolyData>       pd   = surf -> GetOutput();

        // clip the poly data
        vtkSmartPointer<vtkClipPolyData> clipper = vtkSmartPointer<vtkClipPolyData>::New();
        clipper -> SetInput( pd );
        clipper -> InsideOutOn();
        clipper -> SetClipFunction( box );
        clipper -> Update();

        // write new output file
        vtkSmartPointer<vtkPolyData> clippedGrid = clipper -> GetOutput();
        vtkSmartPointer<vtkXMLPolyDataWriter> writer = vtkSmartPointer<vtkXMLPolyDataWriter>::New();
        writer -> SetFileName( outFileName.c_str() );
        writer -> SetInput( clippedGrid );
        writer -> SetDataModeToAppended();
        writer -> Write();

    }
    else {
        // go directly via the unstructured grid

        // clip the grid
        vtkSmartPointer<vtkClipDataSet> clipper = vtkSmartPointer<vtkClipDataSet>::New();
        clipper -> SetInput( grid );
        clipper -> InsideOutOn();
        clipper -> SetClipFunction( box );
        clipper -> Update();

        // write new output file
        vtkSmartPointer<vtkUnstructuredGrid> clippedGrid = clipper -> GetOutput();
        vtkSmartPointer<vtkXMLUnstructuredGridWriter> writer = vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();    
        writer -> SetFileName( outFileName.c_str() );
        writer -> SetInput( clippedGrid );
        writer -> SetDataModeToAppended();
        writer -> Write();
    }
    
    return 0;
}
