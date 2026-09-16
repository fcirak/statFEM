#include <vtkSmartPointer.h>
#include <vtkXMLUnstructuredGridReader.h>
#include <vtkXMLUnstructuredGridWriter.h>
#include <vtkUnstructuredGrid.h>
#include <vtkMergeCells.h>

#include <vector>
#include <boost/lexical_cast.hpp>
#include <corlib/UniqueFilename.hpp>

//------------------------------------------------------------------------------

int main( int argc, char * argv[] )
{
    if ( argc != 5 ) {
        std::cout << "Usage: " << argv[0] 
                  << " string1  lower  upper  string2 " << std::endl
                  << " with: " << std::endl
                  << "  string1  - begin of file name until number range " << std::endl
                  << "  lower    - lower string of number range " << std::endl
                  << "  upper    - upper string of number range " << std::endl
                  << "  string2  - end of file name following the number range " << std::endl
                  << std::endl
                  << " Example call:  " << argv[0] 
                  << "  basename.  000  010  -0340.vtu " << std::endl
                  << "  will merge the 11 files  " << std::endl
                  << "  basename.000-0340.vtu ... basename.010-0340.vtu " << std::endl
                  << "  into:  basename.merge-0340.vtu " << std::endl;

        return 1;
    }


    // read from command line
    const std::string prefix = boost::lexical_cast<std::string>( argv[1] );
    const std::string suffix = boost::lexical_cast<std::string>( argv[4] );
    const unsigned lowerNum  = boost::lexical_cast<unsigned>(    argv[2] );
    const unsigned upperNum  = boost::lexical_cast<unsigned>(    argv[3] );
    const std::string upperString = boost::lexical_cast<std::string>( argv[3] );
    const unsigned numeralWidth = upperString.size();

    // container for grid pointers
    std::vector<vtkSmartPointer<vtkUnstructuredGrid> > grids;

    // iterate over all files in order to get total numbers
    unsigned numGrids  = 0;
    unsigned numPoints = 0;
    unsigned numCells  = 0;

    for ( unsigned gridNum = lowerNum; gridNum <= upperNum; ++ gridNum ) {
        
        // generate the filename
        const std::string numString = corlib::UniqueFilename::convertNumberToString( gridNum, numeralWidth );
        const std::string filename  = prefix + numString + suffix;

        std::cout << " Reading " << filename << std::endl;

        //read XML vtu file
        vtkSmartPointer<vtkXMLUnstructuredGridReader> reader =
            vtkSmartPointer<vtkXMLUnstructuredGridReader>::New();
        reader -> SetFileName( filename.c_str() );
        reader -> Update();

        // get access to the unstructured grid
        vtkSmartPointer<vtkUnstructuredGrid> unstructuredGrid = reader -> GetOutput();

        // get header data
        numGrids  ++;
        numPoints += unstructuredGrid -> GetNumberOfPoints();
        numCells  += unstructuredGrid -> GetNumberOfCells();

        // store grid
        grids.push_back( unstructuredGrid );

    }

    // use merge cells
    vtkSmartPointer<vtkMergeCells> cellMerge = vtkSmartPointer<vtkMergeCells>::New();

    // strange hack
    cellMerge -> SetPointMergeTolerance( std::numeric_limits<double>::min() );

    // allocate data for the merged grid
    vtkSmartPointer<vtkUnstructuredGrid> gridUnion = vtkSmartPointer<vtkUnstructuredGrid>::New();

    // first set total data upper bounds
    cellMerge -> SetUnstructuredGrid( gridUnion );
    cellMerge -> SetTotalNumberOfDataSets( numGrids );
    cellMerge -> SetTotalNumberOfPoints(   numPoints );
    cellMerge -> SetTotalNumberOfCells(    numCells );

    // do the merging
    std::cout << " Merging ... " << std::endl;
    for ( unsigned g = 0; g < numGrids; g ++ ) {
        cellMerge -> MergeDataSet( grids[g] );
    }

    // finalize
    cellMerge -> Finish();

    // write output
    const std::string outputFilename = prefix + "merge" + suffix;

    gridUnion = cellMerge -> GetUnstructuredGrid();

    // write XML binary vtu file
    vtkSmartPointer<vtkXMLUnstructuredGridWriter> vtuWriter = 
        vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
    vtuWriter -> SetFileName( outputFilename.c_str() );
    vtuWriter -> SetInput( gridUnion );
    vtuWriter -> SetDataModeToAppended();
    vtuWriter -> EncodeAppendedDataOff();
    vtuWriter -> Write();


    return 0;
}
