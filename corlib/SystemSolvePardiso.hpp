// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//
//                 Computational Structural Mechanics Lab
//                         University of Cambridge
//
// This software is copyrighted and all rights are retained by the University of
// Cambridge. All use, disclosure, and/or reproduction of any part not expressly
// authorised by the University of Cambridge is prohibited. (C) 2010
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

//! @file SystemSolvePardiso.hpp

#ifndef corlib_systemsolvepardiso_h
#define corlib_systemsolvepardiso_h

//------------------------------------------------------------------------------
//! system includes
#include <ostream>
#include <vector>
#include <set>
#include <cmath>
//! boost includes
#include <boost/numeric/ublas/vector.hpp>
//! corlib includes
#include <corlib/Triplet.hpp>
#include <corlib/linalg.hpp>
#include <corlib/verify.hpp>
#include <corlib/linalg.hpp>

//! MKL header
#include <mkl.h>

#if defined(MKL_ILP64)
#define MKL_INT long long
#else
#define MKL_INT int
#endif

//------------------------------------------------------------------------------
namespace corlib{

    class SystemSolvePardiso;

    namespace ublas = boost::numeric::ublas;
}

//------------------------------------------------------------------------------
/** \brief Solve A x = b with Pardiso
 *  \details This solver interfaces with the shared-memory parallel solver
 *  Pardiso (in the Intel MKL).
 */
class corlib::SystemSolvePardiso
{
private:
    typedef ublas::matrix< double >                         MatrixLocal_;
    typedef ublas::vector< double >                         VectorLocal_;
    typedef std::vector< unsigned >                         VectorUInt_;
    typedef std::vector< std::pair<unsigned,double> >       VectorConstraint_;

    //! Solution mode for Pardiso solver
    enum SolveMode_ {
        COMPLETE,
        FACTORISE,
        BACKSUBSTITUTE,
        RELEASE
    };
    
public:
    //! Construct with number of dofs (= matrix size) 
    SystemSolvePardiso( const unsigned numDofs )
        : numDofs_( numDofs )
    {
        F_.resize( numDofs_ );
        F_.clear();
    }


    //! @name Insertion routines
    //@{
    void insertToMatrix( const MatrixLocal_ & eStiff, 
                         const VectorUInt_ & rowIndices,
                         const VectorUInt_ & colIndices );
    
    void insertToRhs( const VectorLocal_ & eForce, const VectorUInt_ & dofIndices );
    //@}

    //! Debug write of the matrix
    void writeMatrix( std::ostream & out ) { K_.write( out ); }
    //! Debug write of the right hand side
    void writeRhs(    std::ostream & out ) { corlib::writeVector( out, F_ ); }

    //! @name norms
    //@{
    double normRhs( ) const { return norm_2( F_ ); }
    double normSol( ) const { return norm_2( F_ ); }
    //@}

    //! @name Constraints routines
    //@{
    void applyConstraints( VectorConstraint_ & constraints, const double factor, 
                           const double scalar );

    void fixDOF( const unsigned index, const double scalar = 1. );

protected:
    //! Check if the constraints have changed
    //!
    //! \return True if constraint DOFs did not change
    bool updateConstraintDofs_();
    //@}

public:
    //! To be compliant with SystemSolvePetsc (does nothing here)
    void finishRhsAssembly( ) { }
    //! To be compliant with SystemSolvePetsc (does nothing here)
    void finishMatAssembly( ) { }
    //! To be compliant with SystemSolvePetsc (does nothing here)
    void finishAssembly( ) { }

    //! SOLVE
    void solveSystem( );
    void solveSystem( const bool reuseFactoredMatrix );

    //! Fill matrix with zeros
    void clearMatrix( )
    { 
        K_.clear( );
        
        if ( not Anz_.empty() )
            pardisosv( RELEASE, colInd_, rowPtr_, Anz_ );
        colInd_.clear();
        rowPtr_.clear();
        Anz_.clear();
    }
    //! Fill RHS with zeros
    void clearRhs( ) {   F_.clear( ); }
  
    //! Accessor to solution
    void giveSolution( const VectorUInt_ & dofIndices, VectorLocal_ & result ) const;
    
    //! deprecated
    unsigned numNonZeros() const { return K_.nnz(); }


private:
    //! Interface to linear solution with Pardiso
    void pardisosv( const SolveMode_ solveMode,
                    std::vector<int> &, 
                    std::vector<int> &, 
                    std::vector<double> & );

private:
    const unsigned numDofs_; //!< number of dofs (= matrix size)
    corlib::Triplet<int,double> K_; //!< global stiffness matrix (stored in a triplet)
    ublas::vector<double>  F_;      //!< global force vector

private:
    //! @name These variables are only necessary here for applying #solveSystem(bool)
    //@{
    std::set<unsigned>   oldConstrDofs_;  //!< previously constrained dofs
    std::set<unsigned>   newConstrDofs_;  //!< currently constrained dofs

    MKL_INT iparm_[64];   //!< Pardiso control parameters
    void * pt_[64];       //!< Internal solver memory pointer, see #pardisosv() for details

    std::vector<int>      colInd_;   //!< column indices
    std::vector<int>      rowPtr_;   //!< row-ends in #Anz_
    std::vector<double>   Anz_;      //!< before factorisation the non-zeros in matrix
    //@}

};

//------------------------------------------------------------------------------
//! insert element stiffness matrix into global system matrix
void corlib::SystemSolvePardiso::insertToMatrix( const MatrixLocal_ & eStiff, 
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
//! insert element force vector into global force vector
void corlib::SystemSolvePardiso::insertToRhs( const VectorLocal_ & eForce, 
                                              const VectorUInt_ & dofIndices )
{
    for ( unsigned i = 0; i < dofIndices.size( ); i ++ ) {
        F_[ dofIndices[ i ] ] += eForce( i );
    }

    return;
}

//------------------------------------------------------------------------------
//! apply constraints by a sort of a penalty method
void corlib::SystemSolvePardiso::applyConstraints( VectorConstraint_ & constraints,
                                                   const double factor,
                                                   const double scalar )
{
    VectorConstraint_::const_iterator cIter = constraints.begin( );
    VectorConstraint_::const_iterator cEnd  = constraints.end( );

    for ( ; cIter != cEnd; ++ cIter ) {
        const unsigned index = cIter -> first;
        const double value   = factor * (cIter -> second);

        K_.zeroRow( index, scalar );
        F_[ index ] = value * scalar;

        newConstrDofs_.insert( index );
    }

    return;
}

//------------------------------------------------------------------------------
//! Just set the DOF given by 'index' to zero.
void corlib::SystemSolvePardiso::fixDOF( const unsigned index, const double scalar )
{
    K_.zeroRow( index, scalar );
    F_[ index ] = 0.;
    newConstrDofs_.insert( index );
    return;
}

//------------------------------------------------------------------------------
//! Check if the constraints have changed
//!
//! \return True if constraint DOFs did not change
bool corlib::SystemSolvePardiso::updateConstraintDofs_()
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

//------------------------------------------------------------------------------
//! solve the system A x = b directly by Pardiso 
void corlib::SystemSolvePardiso::solveSystem( )
{
    FTL_VERIFY( K_.nnz() > 0 );
    FTL_VERIFY( F_.size() > 0 );

    std::vector<double> Anz;
    std::vector<int> colInd, rowPtr;

    // let triplet object fill the arrays
    K_.setCompressedRowArrays( colInd, rowPtr, Anz );
    
    // clear triplet object
    K_.clear();

    // solver
    pardisosv( COMPLETE, colInd, rowPtr, Anz );

    // clean-up
    oldConstrDofs_.clear();
    newConstrDofs_.clear();
    
    return;
}

//------------------------------------------------------------------------------
//! solve the system A x = b directly by Pardiso allow to reuse factored matrix
void corlib::SystemSolvePardiso::solveSystem( const bool reuseFactoredMatrix )
{
    FTL_VERIFY( F_.size() > 0 );

    if ( not reuseFactoredMatrix ) {
        FTL_VERIFY( K_.nnz() > 0 );

        // let triplet object fill the arrays
        K_.setCompressedRowArrays( colInd_, rowPtr_, Anz_ );
    
        // clear triplet object
        K_.clear();
        
        // solver
        pardisosv( FACTORISE, colInd_, rowPtr_, Anz_ );
    }
    else {
        FTL_VERIFY( not colInd_.empty() );
        FTL_VERIFY( not rowPtr_.empty() );
        FTL_VERIFY( not Anz_.empty() );
    }

    // check pattern of constraints
    const bool hasSameConstr = this -> updateConstraintDofs_();
    if ( reuseFactoredMatrix ) FTL_VERIFY( hasSameConstr );

    // back-substitute
    pardisosv( BACKSUBSTITUTE, colInd_, rowPtr_, Anz_ );
    
    return;
}

//------------------------------------------------------------------------------
//! return solution values (access with vectors)
void corlib::SystemSolvePardiso::giveSolution( const VectorUInt_ & dofIndices, 
                                               VectorLocal_    & result ) const
{
    for ( unsigned i = 0; i < dofIndices.size( ); i ++ ) {
        result[ i ] = F_[ dofIndices[i] ];
    }
    return;
}

//-----------------------------------------------------------------------------

#include "SystemSolvePardiso.ipp"

//-----------------------------------------------------------------------------

#endif
