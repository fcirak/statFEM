#include <vtkSmartPointer.h>
#include <vtkXMLUnstructuredGridReader.h>
#include <vtkUnstructuredGrid.h>
#include <vtkPointLocator.h>

#include <boost/lexical_cast.hpp>

#include <corlib/fuzzyEqual.hpp>

//------------------------------------------------------------------------------
void check( const bool expression, const std::string & message )
{
    if ( not expression ) {
        std::cerr << "(EE) " << message << std::endl;
        exit( -1 );
    }
}

//------------------------------------------------------------------------------
int main( int argc, char *argv[] )
{
    //parse command line arguments
    if( argc != 5 ) {
        std::cerr << "Usage: " << argv[0]
                  << " file.vtu   x  y  z" << std::endl;
        return EXIT_FAILURE;
    }

    // files to compare
    const std::string fileName = boost::lexical_cast<std::string>( argv[1] );
    const double      x        = boost::lexical_cast<double>(      argv[2] );
    const double      y        = boost::lexical_cast<double>(      argv[3] );
    const double      z        = boost::lexical_cast<double>(      argv[4] );

    // read data from files
    vtkSmartPointer<vtkXMLUnstructuredGridReader> reader = 
        vtkSmartPointer<vtkXMLUnstructuredGridReader>::New();
    reader -> SetFileName( fileName.c_str() );
    reader -> Update();

    // get the grids
    vtkSmartPointer<vtkUnstructuredGrid> grid = reader -> GetOutput();

    // 3-D point
    const double point[3] = {x, y, z};

    // set up a locator
    vtkSmartPointer<vtkPointLocator> locator = vtkSmartPointer<vtkPointLocator>::New();
    locator -> SetDataSet( grid );
    std::cout << locator -> FindClosestPoint( point ) << std::endl;

    return 0;
}
