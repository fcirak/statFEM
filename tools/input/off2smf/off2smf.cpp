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

#include <iostream>
#include <vector>

#include <corlib/verify.hpp>

//==============================================================================
int main( int argc, char **argv )
{
    if ( argc > 1 ) {
        std::cerr << "Usage: off2smf < input.off > output.smf" << std::endl;
        return 1;
    }

    // read "OFF" key word
    std::string name;    
    std::cin >> name;
    FTL_VERIFY( name == "OFF" );

    // read mesh magnitudes
    unsigned numVertices, numElements, numFacets;
    std::cin >> numVertices >> numElements >> numFacets;

    // read coordinates
    const unsigned dim = 3;
    unsigned numCoor = numVertices * dim;
    std::vector<double> coordinates;
    coordinates.reserve( numCoor );
    for ( unsigned i = 0; i < numCoor; i ++ ) {
	double xyz;
	std::cin >> xyz;
	coordinates.push_back( xyz );
    }
   
    // determine number of nodes per element
    unsigned numVtxPerElem;
    std::cin >> numVtxPerElem;
    FTL_VERIFY( ( numVtxPerElem == 3 ) or ( numVtxPerElem == 4 ) );

    // read connectivity
    std::vector<int> connectivity;
    connectivity.reserve( numElements * numVtxPerElem );
    for ( unsigned e = 0; e < numElements; e ++ ) {
        // skim number of nodes; was already done for very first element
        if ( e > 0 ) {
            unsigned numVtxPerElemCurrent;
            std::cin >> numVtxPerElemCurrent;
        }

        // read element's nodes
        for ( unsigned j = 0; j < numVtxPerElem; j ++ ) {
            int numVtx;
            std::cin >> numVtx;
            connectivity.push_back( numVtx );
        }
    }
    
    // write smf header
    if ( numVtxPerElem == 3 ) {
	std::cout << "! elementShape triangle " << std::endl;
	std::cout << "! elementNumPoints 3 " << std::endl;
    }
    else if ( numVtxPerElem == 4 ) {
	std::cout << "! elementShape quadrilateral " << std::endl;
	std::cout << "! elementNumPoints 4 " << std::endl;
    }

    // write mesh magnitudes
    std::cout << numVertices << " " << numElements << std::endl;

    // write vertex coordinates 
    std::vector<double>::iterator itc = coordinates.begin();
    for ( int i = 0; i < numVertices; ++i ) {
	for ( int j = 0; j < dim; ++j ) {
	    std::cout << *(itc++) << " ";
	}
	std::cout << std::endl;
    }

    // write connectivity
    std::vector<int>::iterator itcon = connectivity.begin();
    for ( int i = 0; i < numElements; ++i ) {
	for ( int j = 0; j < numVtxPerElem; ++j ) { 
	    std::cout << *(itcon++) << " ";	
	}
	std::cout << std::endl;	
    }
    
    return 0;
}
