// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//
//                 Computational Structural Mechanics Lab
//                         University of Cambridge
//
// This software is copyrighted and all rights are retained by the University of
// Cambridge. All use, disclosure, and/or reproduction of any part not expressly
// authorised by the University of Cambridge is prohibited. (C) 2010 --
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

//! @file ioHelpers.hpp

#ifndef del2_apps_reference_laplaceRHSinfer_iohelpers_h
#define del2_apps_reference_laplaceRHSinfer_iohelpers_h

#include <vector>
//#include <string>
#include <fstream>
#include <istream>
#include <algorithm>

#include <Eigen/Core>

#include <corlib/verify.hpp>

namespace apps{
	//! read data from a given input .dat file and store it as Eigen's vector
	//------------------------------------------------------------------------------
    void readVecData ( const std::string & fileName, Eigen::VectorXd & data )
    {
        std::ifstream fileStream( fileName );
        FTL_VERIFY( fileStream.is_open( ) );
        unsigned dataSize = data.rows( );

        std::vector<double> tmp;
        double x;
        while ( fileStream >> x ) {
            tmp.push_back( x );
        }
        fileStream.close( );

        FTL_VERIFY( tmp.size( ) == dataSize );

        for ( unsigned i = 0; i < dataSize; ++ i ) {
            data( i ) = tmp[ i ];
        }

        return;
    }

	//! read data from a given input .dat file and store it as Eigen's vector
	//------------------------------------------------------------------------------
    template <typename VECDIM>
    void readCoordData ( const std::string & fileName, std::vector<VECDIM> & data )
    {
        std::ifstream fileStream( fileName );
        FTL_VERIFY( fileStream.is_open() );

        VECDIM coord;
        const auto dim = VECDIM::RowsAtCompileTime;
        FTL_VERIFY( dim > 0 );

        // read the lines from the file
        std::string line;
        while ( std::getline( fileStream, line ) ) {
            // split the line into dimension coordinates
            std::stringstream lineStream( line );
            std::istream_iterator<std::string> start( lineStream ), end;
            std::vector<std::string> coorStr( start, end );
            FTL_VERIFY( coorStr.size() == dim );

            // assign the coordinates according to the dimension
            if ( dim >= 1 ) coord[0] = atof(coorStr[0].c_str());
            if ( dim >= 2 ) coord[1] = atof(coorStr[1].c_str());
            if ( dim >= 3 ) coord[2] = atof(coorStr[2].c_str());
            data.push_back( coord );
        }
        fileStream.close();

        return;
    }

	//------------------------------------------------------------------------------
	// read non-square matrix: Phi
    void readMatrixPhi ( const std::string & fileName, Eigen::MatrixXd & matPhi,
                         const unsigned & numDofs )
    {
        std::ifstream inputStream( fileName.c_str( ) );
        FTL_VERIFY( inputStream.is_open( ) );

        //! create a triplet list
        std::vector < Eigen::Triplet<double> > tripletList;
        std::vector<int> rowIndices, colIndices;

        int rowIndex, colIndex;
        double dataVal;

        while ( inputStream >> rowIndex >> colIndex >> dataVal ) {
            rowIndices.push_back( rowIndex );
            colIndices.push_back( colIndex );
            tripletList.push_back(
                    Eigen::Triplet<double>( rowIndex, colIndex, dataVal ) );
        }

        inputStream.close( );

        int rowSize = *std::max_element( rowIndices.begin( ),
                                         rowIndices.end( ) ) + 1;
        int colSize = *std::max_element( colIndices.begin( ),
                                         colIndices.end( ) ) + 1;

        //! construct sparse matrix
        Eigen::SparseMatrix<double> matSparse;
        matSparse.resize( rowSize, numDofs );
        matSparse.setFromTriplets( tripletList.begin( ), tripletList.end( ) );
        matPhi.resize( rowSize, numDofs );
        matPhi = Eigen::MatrixXd( matSparse );

        return;
    }

    //! read data from a given fileName.dat file and store it as Eigen's matrix
    //------------------------------------------------------------------------------
    void readMatrixData ( const std::string & fileName, Eigen::MatrixXd & data )
    {
        std::ifstream fileStream( fileName );
        FTL_VERIFY( fileStream.is_open( ) );
        assert( fileStream.good( ) );

        std::vector<double> arrayData;
        unsigned numRow, numCol;
        double x;

        fileStream >> numRow >> numCol;
        while ( fileStream >> x ) {
            arrayData.push_back( x );
        }
        fileStream.close( );
        FTL_VERIFY( arrayData.size( ) == numRow * numCol );

        data.resize( numRow, numCol );
        for ( unsigned i = 0; i < numRow; ++ i ) {
            for ( unsigned j = 0; j < numCol; ++ j ) {
                data( i, j ) = arrayData[ ( i * numCol ) + j ];
            }
        }

        return;
    }

	//! read parameters and insert them in a std::map
    //------------------------------------------------------------------------------
    void readParamToMap ( const std::string & fileName,
                          std::map<std::string, double> & hParam )
    {
        std::ifstream fileStream( fileName );
        FTL_VERIFY( fileStream.is_open( ) );
        assert( fileStream.good( ) );

        std::string key, line;
        double value;

        while( std::getline( fileStream, line ) ) {
            std::istringstream iss( line, std::istringstream::in );

            if ( !line.length( ) ) continue;
            if ( line[0] == '#' ) continue; // Ignore the line starts with #

            while ( iss >> key >> value ) {
                hParam[ key ] = value; // input them into the map
            }

        }

        fileStream.close( );

        return;
    }
}

#endif
