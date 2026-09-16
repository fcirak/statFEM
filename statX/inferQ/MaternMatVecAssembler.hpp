// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//
//                 Computational Structural Mechanics Lab
//                         University of Cambridge
//
// This software is copyrighted and all rights are retained by the University of
// Cambridge. All use, disclosure, and/or reproduction of any part not expressly
// authorised by the University of Cambridge is prohibited. (C) 2010 --
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

//! @file MaternMatVecAssembler.hpp

#ifndef statX_inferQ_maternmatvecassembler_h
#define statX_inferQ_maternmatvecassembler_h

//------------------------------------------------------------------------------
// system includes
#include <vector>
#include <set>
// Eigen includes
#include <Eigen/Core>
#include <Eigen/Sparse>
// corlib includes
#include <corlib/Triplet.hpp>
#include <corlib/verify.hpp>

//------------------------------------------------------------------------------
namespace statX {
    namespace inferQ {
        class MaternMatVecAssembler;
    }
}

//------------------------------------------------------------------------------
/**
 * \brief System assembler for global Matern matrix and vector.
 * \details Provides storage and access to global Matern matrix and vector.
 */
class statX::inferQ::MaternMatVecAssembler
{
public:

    //! Construct with number of points (matrix size)
    MaternMatVecAssembler( const unsigned & numPoints )
      : numPoints_( numPoints ),
        massVec_( numPoints )
    {
        massVec_.setZero( numPoints );
    }

    //! @name Insertion routines
    ///@{

    void insertToMatrix( const Eigen::MatrixXd & eMatern,
                         const std::vector<unsigned> & rowIndices,
                         const std::vector<unsigned> & colIndices );
    
    void insertToRhs( const Eigen::VectorXd & eLumpedMass,
                      const std::vector<unsigned> & dofIndices );
    ///@}

    //! @name Storage management routines
    ///@{

    //! Clear matrix storage (i.e., fill with zeros)
    void clearMaternMatrix( ) { maternMat_.erase( ); return; }

    //! Clear vector storage (i.e., fill with zeros)
    void clearLumpedMass( )   { massVec_.setZero( ); return; }
    ///@}

    //! @name Accessor routines
    ///@{

    //! Give global Matern stiffness matrix
    void giveMaternMatrix( Eigen::SparseMatrix<double> & maternMat );

    //! Give global unit lumped mass vector
    void giveLumpedMassVector( Eigen::VectorXd & massVec ) const {
        massVec = massVec_;
        return;
    }
    ///@}


private:

    void createEigenTriplet_( std::vector<Eigen::Triplet<double> > & triplet,
                              const std::vector<unsigned> & colIndices,
                              const std::vector<unsigned> & rowIndices,
                              const std::vector<double> & val );

private:

    const unsigned                   numPoints_; //!< number of points = matrix size
    corlib::Triplet<unsigned,double> maternMat_; //!< global Matern stiffness matrix
    Eigen::VectorXd                  massVec_;   //!< global lumped mass vector

};

//------------------------------------------------------------------------------
//! Insert element Matern stiffness matrix into global Matern stiffness matrix
void statX::inferQ::MaternMatVecAssembler::
insertToMatrix( const Eigen::MatrixXd & eMatern,
                const std::vector<unsigned> & rowIndices,
                const std::vector<unsigned> & colIndices )
{
    for ( unsigned i = 0; i < rowIndices.size( ); i ++ ) {
        for ( unsigned j = 0; j < colIndices.size( ); j ++ ) {
            maternMat_.insert( rowIndices[ i ], colIndices[ j ],
                               eMatern( i, j ) );
        }
    }

    return;
}

//------------------------------------------------------------------------------
//! Insert element lumped mass vector into global lumped mass vector
void statX::inferQ::MaternMatVecAssembler::
insertToRhs( const Eigen::VectorXd & eLumpedMass,
             const std::vector<unsigned> & dofIndices )
{
    for ( unsigned i = 0; i < dofIndices.size( ); i ++ ) {
        massVec_[ dofIndices[ i ] ] += eLumpedMass( i );
    }

    return;
}

//------------------------------------------------------------------------------
void statX::inferQ::MaternMatVecAssembler::
giveMaternMatrix( Eigen::SparseMatrix<double> & maternMat )
{
    // resize also initialises the matrix to zeros
    maternMat.resize( numPoints_, numPoints_ );

    std::vector<double> Anz;
    std::vector<unsigned> colInd, rowInd;

    // let triplet object fill the arrays
    maternMat_.extractArrays( rowInd, colInd, Anz );

    // create Eigen triplets from corlib triplets
    std::vector< Eigen::Triplet<double> >  tripletVec;
    createEigenTriplet_( tripletVec, colInd, rowInd, Anz );

    // clear the vectors
    std::vector<unsigned>( ).swap( colInd );
    std::vector<unsigned>( ).swap( rowInd );
    std::vector<double>( ).swap( Anz );

    // set Eigen sparse matrix from triplet vector
    maternMat.setFromTriplets( tripletVec.begin( ), tripletVec.end( ) );

    return;
}

//------------------------------------------------------------------------------
//! Obtain a vector of triplets according to Eigen library from corlib::triplet
void statX::inferQ::MaternMatVecAssembler::
createEigenTriplet_ ( std::vector< Eigen::Triplet<double> > & tripletVec,
                      const std::vector<unsigned>  & colIndices,
                      const std::vector<unsigned>  & rowIndices,
                      const std::vector<double>    & val )
{
    auto itr = rowIndices.begin();
    auto itc = colIndices.begin();

    FTL_VERIFY( colIndices.size() == val.size() );
    FTL_VERIFY( rowIndices.size() == val.size() );

    for ( auto v : val ) {
        const Eigen::Triplet<double> entry( *itr ++, *itc ++, v );
        tripletVec.push_back( entry );
    }

    return;
}

#endif
