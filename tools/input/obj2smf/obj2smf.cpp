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
/// obj file to smf file converter. assumes that there is information         //
/// about texture and normal coordinates     		                          //
/// Compile:  make -B                                                         //
///                                                                           //
/// Usage:  ./obj2smf < file.obj > file.smf                                   //
///                                                                           //
/// Note:  - this tool can handle meshes with tri or quad elements 	          //
//------------------------------------------------------------------------------

#include <iostream>
#include <fstream>
#include <sstream>
#include <corlib/Shape.hpp>
#include <corlib/SmfHead.hpp>
#include <corlib/verify.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <vector>

namespace ublas = boost::numeric::ublas;

typedef ublas::bounded_vector<double, 3>			VecDim;
typedef std::vector<VecDim>							Nodes;
typedef std::vector< std::vector<int> >				Facet;

static const unsigned MAXLINE = std::numeric_limits<std::streamsize>::max();

void ignoreLines(std::ifstream & inp, const char delimiter)
{
	std::string garbage;
	while (inp.peek() != EOF and inp.peek() != delimiter)
	{
		inp >> garbage;
	    inp.ignore(MAXLINE, '\n');
	}
}

void readVertices(std::ifstream & inp, Nodes & nodes, const char delimiter)
{
	VecDim  coord;
	std::string word;

	while (inp.peek() != EOF and inp.peek() == delimiter)
	{
		inp >> word;
		inp >> coord[0] >> coord[1] >> coord[2];
		inp.ignore(MAXLINE, '\n');
		nodes.push_back(coord);
	}
}

void readFaces(std::ifstream & inp, Facet & faces, const char delimiter)
{
	while ( true ) {
		std::string line;
		std::getline(inp, line);

		std::istringstream lines( line );

		std::string faceIndicator;
		lines >> faceIndicator;

		if ( faceIndicator == "f" ) {
			std::vector<int> connects;

			int vtxNumber;

			while ( lines >> vtxNumber ) {
				--vtxNumber;
				connects.push_back( vtxNumber );
			}

			faces.push_back( connects );

			char c = inp.peek( );
			if ( c != 'f' ) break;
		}
	}

}

void writeSmf(std::ofstream & out, const Facet & faces, const Nodes & nodes)
{
	//! write numNodes and numfaces
	out << nodes.size() << "  " << faces.size() << std::endl;

	//! write nodal coordinates
	for (int i = 0; i < nodes.size(); ++i) {
		out << nodes[i][0] << "  " << nodes[i][1] << "  " << nodes[i][2] << std::endl;
	}

	//! write faces coordinates
	for (int i = 0; i < faces.size(); ++i) {
		std::vector<int> face = faces[i];
		int numVertices = face.size();

		for (int j = 0; j < numVertices; ++j) {
			out << faces[i][j] << " ";
		}
		out << std::endl;
	}
}


//------------------------------------------------------------------------------
int main(int argc, char **argv)
{
	//! if appropriate arguments are not provided stop
    if ( argc != 3 ) {
        std::cerr << "Usage: obj2smf < input.obj > <output.obj>" << std::endl;
        return 1;
    }

	//-----------------read obj file--------------------------------
    std::string inputData = argv[1];
    std::string outFile   = argv[2];
	std::ifstream inp( inputData.c_str() );
	FTL_VERIFY( inp.is_open() );

	ignoreLines( inp, 'v' );

	Nodes nodes;
	readVertices( inp, nodes, 'v' );

	ignoreLines( inp, 'f' );

	Facet faces;
	readFaces( inp, faces, 'f' );

	inp.close();

	int numVerticesPelement = faces[0].size();

	std::ofstream file( outFile.c_str() );
	FTL_VERIFY( file.is_open() );

	if ( numVerticesPelement == 3 ) {
		file << "! " << "elementShape" << " " << "triangle" << std::endl;
		file << "! " << "elementNumPoints" << " " << 3 << std::endl;
	} else if ( numVerticesPelement == 4 ) {
		file << "! " << "elementShape" << " " << "quadrilateral" << std::endl;
		file << "! " << "elementNumPoints" << " " << 4 << std::endl;
	} else {
		FTL_VERIFY( false );
	}

	writeSmf( file, faces, nodes );

	file.close();
	file.clear();

    return 0;
}
