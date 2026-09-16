//! System includes
#include <string>

//! Boost includes
#include <boost/lexical_cast.hpp>

//! VTK includes
#include <vtkSmartPointer.h>
#include <vtkXMLDataReader.h>
#include <vtkXMLWriter.h>
#include <vtkXMLUnstructuredGridReader.h>
#include <vtkXMLUnstructuredGridWriter.h>
#include <vtkXMLStructuredGridReader.h>
#include <vtkXMLStructuredGridWriter.h>
#include <vtkUnstructuredGrid.h>
#include <vtkStructuredGrid.h>
#include <vtkPointSet.h>
#include <vtkPointData.h>
#include <vtkCellData.h>
#include <vtkCell.h>

//! Corlib includes
#include <corlib/fuzzyEqual.hpp>
#include <corlib/verify.hpp>

//------------------------------------------------------------------------------
/** Tool to compare two vtu/vts files and write the difference to a new file.
 */
int main( int argc, char *argv[] )
{
    //! Check proper call
    if( argc != 3 ) {
        std::cerr << "Usage: " << argv[0]
                  << " file1.vt{u|s}  file2.vt{u|s} " << std::endl;
        return 0;
    }

    //! files to compare
    const std::string filename1 = boost::lexical_cast<std::string>( argv[1] );
    const std::string filename2 = boost::lexical_cast<std::string>( argv[2] );

    //! check if the type of files is VTS
    bool isVTS;
    if (      filename1.find( ".vts" ) != std::string::npos ) isVTS = true;
    else if ( filename1.find( ".vtu" ) != std::string::npos ) isVTS = false;
    else {
        std::cerr << filename1 << " does not have the right suffix " << std::endl;
        FTL_VERIFY( false );
    }

    //! Sanity check concerning other file
    if (     isVTS ) FTL_VERIFY( filename2.find(".vts") != std::string::npos );
    if ( not isVTS ) FTL_VERIFY( filename2.find(".vtu") != std::string::npos );

    //! output file name
    const std::string suffix = ( isVTS ? ".vts" : ".vtu" );
    const std::string outFilename = "diff" + suffix;

    //! Set up readers and writer according to VTU/VTS choice
    vtkSmartPointer<vtkXMLDataReader> reader1, reader2;
    vtkSmartPointer<vtkXMLWriter> writer;
    if ( isVTS ) {
        reader1 = vtkSmartPointer<vtkXMLStructuredGridReader>::New();
        reader2 = vtkSmartPointer<vtkXMLStructuredGridReader>::New();
        writer  = vtkSmartPointer<vtkXMLStructuredGridWriter>::New();
    }
    else {
        reader1 = vtkSmartPointer<vtkXMLUnstructuredGridReader>::New();
        reader2 = vtkSmartPointer<vtkXMLUnstructuredGridReader>::New();
        writer  = vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
    }

    //! Read
    {
        reader1 -> SetFileName( filename1.c_str() );
        reader1 -> Update();

        reader2 -> SetFileName( filename2.c_str() );
        reader2 -> Update();
    }
 
    //! Get the grids from the reader
    vtkSmartPointer<vtkPointSet> grid1, grid2;
    if ( isVTS ) {
        vtkXMLStructuredGridReader * tmp1 
            = static_cast<vtkXMLStructuredGridReader*>( reader1.GetPointer() );
        grid1 = tmp1 -> GetOutput();
        vtkXMLStructuredGridReader * tmp2 
            = static_cast<vtkXMLStructuredGridReader*>( reader2.GetPointer() );
        grid2 = tmp2 -> GetOutput();
    }
    else {
        vtkXMLUnstructuredGridReader * tmp1 
            = static_cast<vtkXMLUnstructuredGridReader*>( reader1.GetPointer() );
        grid1 = tmp1 -> GetOutput();
        vtkXMLUnstructuredGridReader * tmp2 
            = static_cast<vtkXMLUnstructuredGridReader*>( reader2.GetPointer() );
        grid2 = tmp2 -> GetOutput();
    }

    //--------------------------------------------------------------------------
    // GRID COMPARISON

    //! number of points
    const unsigned nPoints1 = grid1 -> GetNumberOfPoints();
    const unsigned nPoints2 = grid2 -> GetNumberOfPoints();
    FTL_VERIFY( nPoints1 == nPoints2 );

    //! point coordinates
    for ( unsigned n = 0; n < nPoints1; n ++ ) {
        double * p1 = grid1 -> GetPoint( n );
        double * p2 = grid2 -> GetPoint( n );
        for ( unsigned d = 0; d < 3; d ++ ) {
            FTL_VERIFY( corlib::fuzzyEqual( p1[d], p2[d] ) );
        }
    }

    //! number of cells
    const unsigned nCells1 = grid1 -> GetNumberOfCells();
    const unsigned nCells2 = grid2 -> GetNumberOfCells();
    FTL_VERIFY( nCells1 == nCells2 );

    //! cell structure
    for ( unsigned c = 0; c < nCells1; c ++ ) {
        // cell types
        FTL_VERIFY( (grid1 -> GetCellType( c )) == (grid2 -> GetCellType( c ) ) ); 
        // individual cells
        vtkCell * cell1 = grid1 -> GetCell( c );
        vtkCell * cell2 = grid2 -> GetCell( c );
        // number of points
        const unsigned nPointsPerCell1 = cell1 -> GetNumberOfPoints();
        const unsigned nPointsPerCell2 = cell2 -> GetNumberOfPoints();
        FTL_VERIFY( nPointsPerCell1 == nPointsPerCell2 );
        // individual points
        for ( unsigned n = 0; n < nPointsPerCell1; n ++ ) {
            FTL_VERIFY( (cell1 -> GetPointId(n)) == (cell2 -> GetPointId(n) ) );
        }
    }
        
    //--------------------------------------------------------------------------
    // DIFF POINTDATA
    {
        // get point data
        vtkSmartPointer<vtkPointData> pointData1 = grid1 -> GetPointData();
        vtkSmartPointer<vtkPointData> pointData2 = grid2 -> GetPointData();

        // number of arrays
        const unsigned nArrays1 = pointData1 -> GetNumberOfArrays();
        const unsigned nArrays2 = pointData2 -> GetNumberOfArrays();
        FTL_VERIFY( nArrays1 == nArrays2 );

        std::cout << "Comparing " << nArrays1 << " arrays of point data: "
                  << std::endl;
    
        // go through arrays
        for ( unsigned n = 0; n < nArrays1; n ++ ) {
            // get arrays
            vtkDataArray * array1 = pointData1 -> GetArray( n );
            vtkDataArray * array2 = pointData2 -> GetArray( n );
            // array names
            const std::string name1( pointData1 -> GetArrayName( n ) );
            const std::string name2( pointData2 -> GetArrayName( n ) );
            std::cout << n << ": " << name1 << "  vs. " << name2 << std::endl;
        
            // get number of components
            const unsigned nComponents1 = array1 -> GetNumberOfComponents();
            const unsigned nComponents2 = array2 -> GetNumberOfComponents();
            FTL_VERIFY(  nComponents1 == nComponents2 );
            // array size
            const unsigned nTuples1 = array1 -> GetNumberOfTuples();
            const unsigned nTuples2 = array2 -> GetNumberOfTuples();
            FTL_VERIFY( nTuples1 == nTuples2 );
            // go through array
            for ( unsigned t = 0; t < nTuples1; t ++ ) {
                double * tuple1 = array1 -> GetTuple( t );
                double * tuple2 = array2 -> GetTuple( t );
                // diff each component
                for ( unsigned c = 0; c < nComponents1; c ++ ) {
                    tuple1[c] = tuple1[c] - tuple2[c];
                }
                array1 -> SetTuple( t, tuple1 );
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
        FTL_VERIFY( nArrays1 == nArrays2 );
    
        std::cout << "Comparing " << nArrays1 << " arrays of cell data: "
                  << std::endl;

        // go through arrays
        for ( unsigned n = 0; n < nArrays1; n ++ ) {
            // get arrays
            vtkDataArray * array1 = cellData1 -> GetArray( n );
            vtkDataArray * array2 = cellData2 -> GetArray( n );
            // array names
            const std::string name1( cellData1 -> GetArrayName( n ) );
            const std::string name2( cellData2 -> GetArrayName( n ) );
            std::cout << n << ": " << name1 << "  vs. " << name2 << std::endl;
        
            // get number of components
            const unsigned nComponents1 = array1 -> GetNumberOfComponents();
            const unsigned nComponents2 = array2 -> GetNumberOfComponents();
            FTL_VERIFY( nComponents1 == nComponents2 );
            // array size
            const unsigned nTuples1 = array1 -> GetNumberOfTuples();
            const unsigned nTuples2 = array2 -> GetNumberOfTuples();
            FTL_VERIFY( nTuples1 == nTuples2 );
            // go through array
            for ( unsigned t = 0; t < nTuples1; t ++ ) {
                double * tuple1 = array1 -> GetTuple( t );
                double * tuple2 = array2 -> GetTuple( t );
                // check each component
                for ( unsigned c = 0; c < nComponents1; c ++ ) {
                    tuple1[c] = tuple1[c] - tuple2[c];
                }
                array1 -> SetTuple( t, tuple1 );
            }
        }
    }    

    //--------------------------------------------------------------------------
    // WRITE NEW FILE
    writer -> SetFileName( outFilename.c_str() );
    writer -> SetInputConnection( reader1 -> GetOutputPort() );
    writer -> SetDataModeToAppended();
    writer -> EncodeAppendedDataOff();
    writer -> Write();
    
    return 0;
}
