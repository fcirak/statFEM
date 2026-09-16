// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//
//                 Computational Structural Mechanics Lab
//                         University of Cambridge
//
// This software is copyrighted and all rights are retained by the University of
// Cambridge. All use, disclosure, and/or reproduction of any part not expressly
// authorised by the University of Cambridge is prohibited. (C) 2010
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

//! @file SystemSolveEigenSparse.hpp

#ifndef corlib_systemsolveeigensparse_h
#define corlib_systemsolveeigensparse_h

//------------------------------------------------------------------------------
//! system includes
#include <ostream>
#include <vector>
#include <set>
//! Eigen includes
#include <Eigen/Core>
#include <Eigen/Sparse>
#include <Eigen/SparseLU>
//! corlib includes
#include <corlib/Triplet.hpp>
#include <corlib/linalg.hpp>
#include <corlib/verify.hpp>
#include <corlib/eigenX.hpp>

//------------------------------------------------------------------------------
namespace corlib {
    class SystemSolveEigenSparse;
}

//------------------------------------------------------------------------------
/** \brief Solution of A x = b by means of Eigen sparse solver
 *  \details The system A x = b, with A possibly non-symmetric is solved for x.
 */
class corlib::SystemSolveEigenSparse
{
private:
    //! @name Convenience typedefs
    //@{
    typedef Eigen::MatrixXd                               MatrixLocal_;
    typedef Eigen::SparseMatrix< double >                 SparseMatrix_;
    typedef Eigen::VectorXd                               VectorLocal_;
    typedef std::vector< unsigned >                       VectorUInt_;
    typedef std::vector< std::pair<unsigned,double> >     VectorConstraint_;
    //@}

public:
    //! Construct with number of dofs (matrix size)
    SystemSolveEigenSparse( const unsigned & numDofs )
      : numDofs_( numDofs ),
        A_( numDofs, numDofs ),
        F_( numDofs )
    {
        F_.setZero( numDofs );
    }

    //! Destructor
    ~SystemSolveEigenSparse( ) { this->clearMatrix( ); }

    //! @name Insertion routines
    //@{
    void insertToMatrix( const MatrixLocal_ & eStiff, 
                         const VectorUInt_  & rowIndices,
                         const VectorUInt_  & colIndices );
    
    void insertToRhs(  const VectorLocal_ & eForce, 
                       const VectorUInt_ & dofIndices );
    //@}

    //! @name Constraint routines
    //@{
    void applyConstraints( const VectorConstraint_ & constraints, const double factor, 
                           const double scalar );

    void fixDOF( const unsigned index, const double scalar = 1. );
    //@}

    //! Debug writing of matrix
    void writeMatrix( std::ostream & out ) { K_.write( out ); }
    //! Debug writing of right hand side
    void writeRhs(    std::ostream & out ) { corlib::writeVector( out, F_ ); }

    //! @name Norms
    //@{
    double normMatrix( ) { return K_.frobeniusNorm( ); }
    double normRhs( ) const { return  F_.lpNorm<2>( ); }
    double norm1Rhs( ) const { return F_.lpNorm<1>( ); }
    double normSol( ) const { return  X_.lpNorm<2>( ); }
    //@}

    //! To be compliant with SystemSolvePetsc (does nothing here)
    void finishRhsAssembly( ) { }
    //! To be compliant with SystemSolvePetsc (does nothing here)
    void finishMatAssembly( ) { }
    //! To be compliant with SystemSolvePetsc (does nothing here)
    void finishAssembly( ) { }

    //! Solve
    void solveSystem( );
    void solveSystem( const bool reuseFactoredMatrix );
    void solveSystemTranspose( const Eigen::VectorXd & rhsVec,
                               Eigen::VectorXd & uVec,
                               const bool reuseFactoredMatrix );

    //! Clear matrix storage (i.e., fill with zeros)
    void clearMatrix( ) { K_.erase( ); A_.setZero( ); }
    //! Clear RHS storage (i.e., fill with zeros)
    void clearRhs( )    { F_.setZero( ); }

    //! Solution accessor
    void giveSolution( const VectorUInt_ & dofIndices, VectorLocal_ & result ) const;
    //! Force accessor
    void giveForce( const VectorUInt_ & dofIndices, VectorLocal_ & result ) const;
    //@}

    //! deprecated?
    unsigned numNonZeros() const { return K_.nnz(); }

protected:
    //! Check if the constraints have changed
    //!
    //! \return True if constraint DOFs did not change
    bool updateConstraintDofs_()
    {
        std::vector<unsigned> diffDof;
        std::set_symmetric_difference( newConstrDofs_.begin(), newConstrDofs_.end(),
                                       oldConstrDofs_.begin(), oldConstrDofs_.end(),
                                       std::back_inserter( diffDof ) );
        const bool areSameConstraints = diffDof.empty();

        // update constraint dofs sets
        std::swap( newConstrDofs_, oldConstrDofs_ );
        newConstrDofs_.clear();

        return areSameConstraints;
    }
    //@}

private:
    const unsigned                   numDofs_;  //!< number of dofs = matrix size
    corlib::Triplet<unsigned,double> K_;        //!< global stiffness matrix (stored in a triplet)
    VectorLocal_                     F_;        //!< global force vector

    SparseMatrix_     A_;          //!< sparse matrix on LHS
    VectorLocal_      rhs_;        //!< right-hand-side
    VectorLocal_      X_;          //!< solution vector

    Eigen::SparseLU< SparseMatrix_ >     solver_;   //!< LU decomposition

    std::set< unsigned >     oldConstrDofs_;  //!< previously constrained dofs
    std::set< unsigned >     newConstrDofs_;  //!< currently constrained dofs
    
    void createEigenTriplet_( std::vector<Eigen::Triplet<double> > &triplet,
                              const VectorUInt_ &colIndices,
                              const VectorUInt_ &rowIndices,
                              const std::vector<double> &val );

    void createEigenTripletTranspose_( std::vector<Eigen::Triplet<double> > &triplet,
                                       const VectorUInt_ &colIndices,
                                       const VectorUInt_ &rowIndices,
                                       const std::vector<double> &val );
};

//------------------------------------------------------------------------------
//! insert element stiffness matrix into global system matrix
void corlib::SystemSolveEigenSparse::insertToMatrix( const MatrixLocal_ & eStiff,
                                                     const VectorUInt_ & rowIndices,
                                                     const VectorUInt_ & colIndices )
{
    for ( unsigned i = 0; i < rowIndices.size( ); i ++ ) {
        for ( unsigned j = 0; j < colIndices.size( ); j ++ ) {
            K_.insert( rowIndices[i], colIndices[j], eStiff(i,j) );
        }
    }

    return;
}

//------------------------------------------------------------------------------
//! insert element stiffness matrix into global system matrix
void corlib::SystemSolveEigenSparse::insertToRhs( const VectorLocal_ & eForce,
                                                  const VectorUInt_ & dofIndices )
{
    for ( unsigned i = 0; i < dofIndices.size( ); i ++ ) {
        F_[ dofIndices[ i ] ] += eForce( i );
    }
    return;
}

//------------------------------------------------------------------------------
//! apply constraints by a sort of a penalty method
void corlib::SystemSolveEigenSparse::applyConstraints( const VectorConstraint_ & constraints,
                                                       const double factor,
                                                       const double scalar )
{
    VectorConstraint_::const_iterator cIter = constraints.begin( );
    VectorConstraint_::const_iterator cEnd  = constraints.end( );

    for ( ; cIter != cEnd; ++ cIter ) {

        const unsigned index = cIter -> first;
        const double   value = factor * (cIter -> second);

        K_.zeroRow( index, scalar  );
        F_[ index ] = value * scalar;
        newConstrDofs_.insert( index );
    }

    return;
}

//------------------------------------------------------------------------------
//! Just set the DOF given by 'index' to zero.
void corlib::SystemSolveEigenSparse::fixDOF( const unsigned index, const double scalar )
{
    K_.zeroRow( index, scalar );
    F_[ index ] = 0.;
    newConstrDofs_.insert( index );
    return;
}

//------------------------------------------------------------------------------
//! Obtain a vector of triplets according to Eigen library from corlib::triplet
void corlib::SystemSolveEigenSparse::createEigenTriplet_( std::vector<Eigen::Triplet< double > > & triplet,
                                                          const VectorUInt_ & colIndices,
                                                          const VectorUInt_ & rowIndices,
                                                          const std::vector< double > & val )
{
    auto itr = rowIndices.begin();
    auto itc = colIndices.begin();

    FTL_VERIFY( colIndices.size() == val.size() );
    FTL_VERIFY( rowIndices.size() == val.size() );

    for (auto v : val)
    {
        const Eigen::Triplet<double> entry(*itr++, *itc++, v);
        triplet.push_back(entry);
    }

    return;
}

//------------------------------------------------------------------------------
//! Obtain a vector of triplets according to Eigen library from corlib::triplet
void corlib::SystemSolveEigenSparse::createEigenTripletTranspose_( std::vector<Eigen::Triplet< double > > & triplet,
                                                                   const VectorUInt_ & colIndices,
                                                                   const VectorUInt_ & rowIndices,
                                                                   const std::vector< double > & val )
{
    auto itr = rowIndices.begin();
    auto itc = colIndices.begin();

    FTL_VERIFY( colIndices.size() == val.size() );
    FTL_VERIFY( rowIndices.size() == val.size() );

    for (auto v : val)
    {
        const Eigen::Triplet<double> entry( *itc++, *itr++, v );
        triplet.push_back( entry );
    }

    return;
}

//------------------------------------------------------------------------------
//! solve the system A x = b directly by means of an LU decomposition
void corlib::SystemSolveEigenSparse::solveSystem( )
{
    std::vector<double> Anz;
    VectorUInt_ colInd, rowInd;

    // let triplet object fill the arrays
    K_.extractArrays( rowInd, colInd, Anz );

    // create Eigen triplets from corlib triplets
    std::vector< Eigen::Triplet<double> >  triplets;
    createEigenTriplet_( triplets, colInd, rowInd, Anz );

    // clear the arrays
    std::vector<unsigned> colIndTmp, rowIndTmp;
    colInd.swap( colIndTmp );
    rowInd.swap( rowIndTmp );
    std::vector<double> anzTmp;
    Anz.swap( anzTmp );

    // set Eigen sparse matrix from triplets
    A_.setFromTriplets( triplets.begin(), triplets.end() );

    // clear Eigen triplets
    triplets.resize( 0 );
    triplets.shrink_to_fit( );

    solver_.compute( A_ );
    
    X_ = solver_.solve ( F_ );

    FTL_VERIFY_DESCRIPTIVE( (solver_.info() == Eigen::Success ),
                             "Sparse solver failed to solve with info=%s\n",
                              solver_.lastErrorMessage().c_str( ) );

    oldConstrDofs_.clear();
    newConstrDofs_.clear();

    return;
}

//------------------------------------------------------------------------------
//! solve the system A x = b by means of reusing the factored matrix
void corlib::SystemSolveEigenSparse::solveSystem( const bool reuseFactoredMatrix )
{
    // not reusing factored matrix
    if ( not reuseFactoredMatrix ) {
        
        std::vector<double> Anz;
        VectorUInt_ colInd, rowInd;
        
        // let triplet object fill the arrays
        K_.extractArrays( rowInd, colInd, Anz );

        // create Eigen triplets from corlib triplets
        std::vector< Eigen::Triplet<double> >  triplets;
        createEigenTriplet_( triplets, colInd, rowInd, Anz );

        // clear the arrays
        std::vector<unsigned> colIndTmp, rowIndTmp;
        colInd.swap( colIndTmp );
        rowInd.swap( rowIndTmp );
        std::vector<double> anzTmp;
        Anz.swap( anzTmp );

        // set Eigen sparse matrix from triplets
        A_.setFromTriplets( triplets.begin(), triplets.end() );

        // clear Eigen triplets
        triplets.resize( 0 );
        triplets.shrink_to_fit( );

        // factorise
        solver_.compute( A_ );

    }

    // check constraint DOFs
    const bool hasSameConstr = this -> updateConstraintDofs_();
    if ( reuseFactoredMatrix ) FTL_VERIFY( hasSameConstr );

    X_ = solver_.solve( F_ );

    FTL_VERIFY_DESCRIPTIVE( (solver_.info() == Eigen::Success ),
                             "Sparse solver failed to solve with info=%s\n",
                              solver_.lastErrorMessage().c_str( ) );

    return;
}

//------------------------------------------------------------------------------
//! solve the system A^T uVec = rhsVec by means of reusing the factored matrix
void corlib::SystemSolveEigenSparse::solveSystemTranspose( const Eigen::VectorXd & rhsVec, Eigen::VectorXd & uVec, const bool reuseFactoredMatrix )
{
    // not reusing factored matrix
    if ( not reuseFactoredMatrix ) {

        std::vector<double> Anz;
        VectorUInt_ colInd, rowInd;

        // let triplet object fill the arrays
        K_.extractArrays( rowInd, colInd, Anz );

        // create Eigen triplets from corlib triplets
        std::vector< Eigen::Triplet<double> >  triplets;
        createEigenTripletTranspose_( triplets, colInd, rowInd, Anz );

        // clear the arrays
        std::vector<unsigned> colIndTmp, rowIndTmp;
        colInd.swap( colIndTmp );
        rowInd.swap( rowIndTmp );
        std::vector<double> anzTmp;
        Anz.swap( anzTmp );

        // set Eigen sparse matrix from triplets
        A_.setFromTriplets( triplets.begin(), triplets.end() );

        // clear Eigen triplets
        triplets.resize( 0 );
        triplets.shrink_to_fit( );

        // factorise
        solver_.compute( A_ );
    }

    uVec = solver_.solve( rhsVec );

    FTL_VERIFY_DESCRIPTIVE( (solver_.info() == Eigen::Success ),
                             "Sparse solver failed to solve with info=%s\n",
                              solver_.lastErrorMessage().c_str( ) );

    return;
}

//------------------------------------------------------------------------------
//! return solution values, access with vectors
void corlib::SystemSolveEigenSparse::giveSolution( const VectorUInt_ & dofIndices,
                                                   VectorLocal_      & result ) const
{
    for ( unsigned i = 0; i < dofIndices.size( ); i ++ ) {
        result[ i ] = X_[ dofIndices[i] ];
    }
    return;
}

//------------------------------------------------------------------------------
//! return values of the force vector
void corlib::SystemSolveEigenSparse::giveForce( const VectorUInt_ & dofIndices,
                                                VectorLocal_      & force ) const
{
    for ( unsigned i = 0; i < dofIndices.size( ); i ++ ) {
        force[ i ] = F_[ dofIndices[i] ];
    }
    return;
}

#endif
