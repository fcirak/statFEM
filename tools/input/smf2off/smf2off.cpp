//
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//
//                 Computational Structural Mechanics Lab
//                         University of Cambridge
//
// This software is copyrighted and all rights are retained by the CSMLab. All
// use, disclosure, and/or reproduction of any part not expressly authorized by
// F Cirak is prohibited. (C) 2010.
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//

//------------------------------------------------------------------------------
/// small converter for files in the smf-format to off-files                  //
/// Compile:  make -B                                                         //
///                                                                           //
/// Usage:  ./smf2off < file.smf > file.off                                   //
///                                                                           //
/// Note:  - this tool can handle only meshes with the same element type      //
///        - the element type must be either                                  //
///            *  3-noded (linear) triangles                                  //
///            *  4-noded (linear) quadrilaterals                             //
//------------------------------------------------------------------------------

#include <iostream>
#include <corlib/Shape.hpp>
#include <corlib/SmfHead.hpp>
#include <corlib/verify.hpp>

//------------------------------------------------------------------------------
int main(int argc, char **argv)
{
    if ( argc != 1 ) {
        std::cerr << "Usage: smf2off < input.smf > output.off" << std::endl;
        return 1;
    }

    // search header section in SMF file to identify element type
    corlib::SmfHead smfHead;
    smfHead.read( std::cin );
    FTL_VERIFY_DESCRIPTIVE( smfHead.foundHeader(),
                            "SMF file without header. Please add one. Thank you.\n" );
    const enum corlib::shape eleShape = smfHead.giveElementShape();
    const unsigned nPointsPerElement = smfHead.giveElementNumPoints();
    FTL_VERIFY( ( nPointsPerElement == 3 and eleShape == corlib::TRIANGLE ) or
                ( nPointsPerElement == 4 and eleShape == corlib::QUADRILATERAL ) );

    unsigned nNodes, nElements;
    std::cin >> nNodes >> nElements;

    // write header
    std::cout << "OFF" << std::endl
              << nNodes << " " << nElements << " 0 " << std::endl
              << std::endl;

    // read and write nodes
    const unsigned dim = 3;
    for ( unsigned n = 0; n < nNodes; n ++ ) {
        for ( unsigned d = 0; d < dim; d ++ ) {
            double x;
            std::cin >> x;
            std::cout << x << " ";
        }
        std::cout << std::endl;
    }
	
    // read and write cell connectivity
    for ( unsigned e = 0; e < nElements; e ++ ) {
    	std::cout << nPointsPerElement << " ";
        for ( unsigned v = 0; v < nPointsPerElement; v ++ ) {
            int p;
            std::cin >> p;
            std::cout << p << " ";
        }
        std::cout << std::endl;
    }

    return 0;
}
