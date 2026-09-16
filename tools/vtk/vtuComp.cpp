#include <vtkSmartPointer.h>
#include <vtkXMLUnstructuredGridReader.h>
#include <vtkUnstructuredGrid.h>
#include <vtkPointData.h>
#include <vtkCellData.h>
#include <vtkCell.h>

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
    if( argc != 4 ) {
        std::cerr << "Usage: " << argv[0]
                  << " file1.vtu  file2.vtu  tol " << std::endl;
        return EXIT_FAILURE;
    }

    // files to compare
    const std::string fileName1 = boost::lexical_cast<std::string>( argv[1] );
    const std::string fileName2 = boost::lexical_cast<std::string>( argv[2] );

    // comparison tolerance
    const double tol = boost::lexical_cast<double>( argv[3] );

    // read data from files
    vtkSmartPointer<vtkXMLUnstructuredGridReader> reader1 = 
        vtkSmartPointer<vtkXMLUnstructuredGridReader>::New();
    reader1 -> SetFileName( fileName1.c_str() );
    reader1 -> Update();

    vtkSmartPointer<vtkXMLUnstructuredGridReader> reader2 = 
        vtkSmartPointer<vtkXMLUnstructuredGridReader>::New();
    reader2 -> SetFileName( fileName2.c_str() );
    reader2 -> Update();
 
    // get the grids
    vtkSmartPointer<vtkUnstructuredGrid> grid1 = reader1 -> GetOutput();
    vtkSmartPointer<vtkUnstructuredGrid> grid2 = reader2 -> GetOutput();

    //grid1 -> PrintSelf( std::cout, vtkIndent(2) );

    //--------------------------------------------------------------------------
    // GRID COMPARISON

    // number of points
    const unsigned nPoints1 = grid1 -> GetNumberOfPoints();
    const unsigned nPoints2 = grid2 -> GetNumberOfPoints();
    check( nPoints1 == nPoints2, "Number of points does not match" );

    // point coordinates
    for ( unsigned n = 0; n < nPoints1; n ++ ) {
        double * p1 = grid1 -> GetPoint( n );
        double * p2 = grid2 -> GetPoint( n );
        for ( unsigned d = 0; d < 3; d ++ ) {
            check( corlib::fuzzyEqual( p1[d], p2[d] ), "Coordinates do not match" );
        }
    }

    // number of cells
    const unsigned nCells1 = grid1 -> GetNumberOfCells();
    const unsigned nCells2 = grid2 -> GetNumberOfCells();
    check( nCells1 == nCells2, "Number of cells does not match" );

    // cell structure
    for ( unsigned c = 0; c < nCells1; c ++ ) {
        // cell types
        check( (grid1 -> GetCellType( c )) == (grid2 -> GetCellType( c )), 
               "Cell types do not match" );
        // individual cells
        vtkCell * cell1 = grid1 -> GetCell( c );
        vtkCell * cell2 = grid2 -> GetCell( c );
        // number of points
        const unsigned nPointsPerCell1 = cell1 -> GetNumberOfPoints();
        const unsigned nPointsPerCell2 = cell2 -> GetNumberOfPoints();
        check( nPointsPerCell1 == nPointsPerCell2, 
               "Number of points per cell does not match" );
        // individual points
        for ( unsigned n = 0; n < nPointsPerCell1; n ++ ) {
            check( (cell1 -> GetPointId(n)) == (cell2 -> GetPointId(n) ),
                   "Connectivity does not match");
        }
    }
        
    //--------------------------------------------------------------------------
    // CHECK POINTDATA
    {
        // get point data
        vtkSmartPointer<vtkPointData> pointData1 = grid1 -> GetPointData();
        vtkSmartPointer<vtkPointData> pointData2 = grid2 -> GetPointData();

        // number of arrays
        const unsigned nArrays1 = pointData1 -> GetNumberOfArrays();
        const unsigned nArrays2 = pointData2 -> GetNumberOfArrays();
        check( nArrays1 == nArrays2, 
               "Number of arrays in point data does not match" );
    
        // go through arrays
        for ( unsigned n = 0; n < nArrays1; n ++ ) {
            // get arrays
            vtkDataArray * array1 = pointData1 -> GetArray( n );
            vtkDataArray * array2 = pointData2 -> GetArray( n );
            // array names
            const std::string name1( pointData1 -> GetArrayName( n ) );
            const std::string name2( pointData2 -> GetArrayName( n ) );
        
            // get number of components
            const unsigned nComponents1 = array1 -> GetNumberOfComponents();
            const unsigned nComponents2 = array2 -> GetNumberOfComponents();
            check(  nComponents1 == nComponents2, 
                    "Number of components in point data does not match in field "
                    + name1 + " / " + name2 );
            // array size
            const unsigned nTuples1 = array1 -> GetNumberOfTuples();
            const unsigned nTuples2 = array2 -> GetNumberOfTuples();
            check( nTuples1 == nTuples2, 
                   "Number of tuples in point point does not match in field "
                   + name1 + " / " + name2 );
            // go through array
            for ( unsigned t = 0; t < nTuples1; t ++ ) {
                double * tuple1 = array1 -> GetTuple( t );
                double * tuple2 = array2 -> GetTuple( t );
                // check each component
                for ( unsigned c = 0; c < nComponents1; c ++ ) {
                    check( corlib::fuzzyEqual( tuple1[c], tuple2[c], tol ), 
                           "Entries of point point do not match for field "
                           + name1 + " / " + name2 );
                }
            }
        }
    }    
    //--------------------------------------------------------------------------
    // CHECK CELLDATA
    {
        vtkSmartPointer<vtkCellData> cellData1 = grid1 -> GetCellData();
        vtkSmartPointer<vtkCellData> cellData2 = grid2 -> GetCellData();

        // number of arrays
        const unsigned nArrays1 = cellData1 -> GetNumberOfArrays();
        const unsigned nArrays2 = cellData2 -> GetNumberOfArrays();
        check( nArrays1 == nArrays2, 
               "Number of arrays in cell data does not match" );
    
        // go through arrays
        for ( unsigned n = 0; n < nArrays1; n ++ ) {
            // get arrays
            vtkDataArray * array1 = cellData1 -> GetArray( n );
            vtkDataArray * array2 = cellData2 -> GetArray( n );
            // array names
            const std::string name1( cellData1 -> GetArrayName( n ) );
            const std::string name2( cellData2 -> GetArrayName( n ) );
        
            // get number of components
            const unsigned nComponents1 = array1 -> GetNumberOfComponents();
            const unsigned nComponents2 = array2 -> GetNumberOfComponents();
            check(  nComponents1 == nComponents2, 
                    "Number of components in cell data does not match in field "
                    + name1 + " / " + name2 );
            // array size
            const unsigned nTuples1 = array1 -> GetNumberOfTuples();
            const unsigned nTuples2 = array2 -> GetNumberOfTuples();
            check( nTuples1 == nTuples2, 
                   "Number of tuples in cell data does not match in field "
                   + name1 + " / " + name2 );
            // go through array
            for ( unsigned t = 0; t < nTuples1; t ++ ) {
                double * tuple1 = array1 -> GetTuple( t );
                double * tuple2 = array2 -> GetTuple( t );
                // check each component
                for ( unsigned c = 0; c < nComponents1; c ++ ) {
                    check( corlib::fuzzyEqual( tuple1[c], tuple2[c], tol ), 
                           "Entries of cell data do not match for field "
                           + name1 + " / " + name2 );
                }
            }
        }
    }    
    
    return 0;
}
