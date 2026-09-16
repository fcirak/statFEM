//! System includes
#include <string>
#include <iostream>
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
//! Corlib includes
#include <corlib/verify.hpp>

//------------------------------------------------------------------------------
// Converts ASCII-XML vtu/s-files to binary 
//------------------------------------------------------------------------------

int main(int argc, char *argv[])
{
    // print usage message
    if( argc != 2 ) {
        std::cout << "Usage: " << argv[0] << " filename.vt{u|s}" << std::endl
                  << std::endl
                  << "Output: filename.bin.vt{u|s}" << std::endl;
            
        return 0;
    }

    // get input file name
    const std::string filename = boost::lexical_cast<std::string>( argv[1] );
 
    bool isVTS;
    if (      filename.find( ".vts" ) != std::string::npos ) isVTS = true;
    else if ( filename.find( ".vtu" ) != std::string::npos ) isVTS = false;
    else {
        std::cerr << filename << " does not have the right suffix " << std::endl;
        FTL_VERIFY( false );
    }

    // get basename
    const std::string suffix   = ( isVTS ? ".vts" : ".vtu" );
    const std::string basename = filename.substr( 0, filename.find( suffix ) );
    const std::string outputFilename = basename + ".bin" + suffix;

    // create reader and writer
    vtkSmartPointer<vtkXMLDataReader> reader;
    vtkSmartPointer<vtkXMLWriter> writer;
    if ( isVTS ) {
        reader = vtkSmartPointer<vtkXMLStructuredGridReader>::New();
        writer = vtkSmartPointer<vtkXMLStructuredGridWriter>::New();
    }
    else {
        reader = vtkSmartPointer<vtkXMLUnstructuredGridReader>::New();
        writer = vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
    }

    // read file
    reader -> SetFileName( filename.c_str() );
    reader -> Update();

    // write XML binary file
    writer -> SetFileName( outputFilename.c_str() );
    writer -> SetInputConnection( reader->GetOutputPort() );
    writer -> SetDataModeToAppended();
    writer -> EncodeAppendedDataOff();
    writer -> Write();
 
    return 0;
}

 
