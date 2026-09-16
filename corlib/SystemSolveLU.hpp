// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//
//                 Computational Structural Mechanics Lab
//                         University of Cambridge
//
// This software is copyrighted and all rights are retained by the University of
// Cambridge. All use, disclosure, and/or reproduction of any part not expressly
// authorised by the University of Cambridge is prohibited. (C) 2010
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

//! @file SystemSolveLU.hpp

#ifndef corlib_systemsolvelu_h
#define corlib_systemsolvelu_h

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

//! SuperLU header
#include <slu_ddefs.h>

// Following prototype is missing in slu_ddefs.h but available in SuperLU library.
// This work-around may be resolved in the future.
//
// Tested installs:
//    SuperLU 3.1, 4.2 from source
//    SuperLU 3.0 package managed on Ubuntu 10.04
extern "C" {
    extern double dlangs( char *, SuperMatrix * );
}

//------------------------------------------------------------------------------
namespace corlib{

    class SystemSolveLU;

    namespace ublas = boost::numeric::ublas;
}

//------------------------------------------------------------------------------
/** \brief Solution of A x = b by means of SuperLU
 *  \details Using the single-threaded version of the library SuperLU, the
 *  system A x = b, with A possibly non-symmetric is solved for x. 
 */
class corlib::SystemSolveLU
{
private:
    //! @name Convenience typedefs
    //@{
    typedef ublas::matrix< double >                         MatrixLocal_;
    typedef ublas::vector< double >                         VectorLocal_;
    typedef std::vector< unsigned >                         VectorUInt_;
    typedef std::vector< std::pair<unsigned,double> >       VectorConstraint_;
    //@}

public:
    //! Construct with number of dofs (matrix size)
    SystemSolveLU( const unsigned & numDofs  )
      : numDofs_( numDofs ),
        conditionNumberNorm_( '#' ), 
        conditionNumber_( 0.0 )
    {
        F_.resize( numDofs_ ),
        F_.clear();

        this -> initSuperMatrices_();
    }

    //! Destructor
    ~SystemSolveLU() { this->clearSuperMatrices_(); }

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
    void applyConstraints( VectorConstraint_ & constraints, const double factor, 
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
    double normRhs( ) const { return ublas::norm_2( F_ ); }
    double norm1Rhs( ) const { return ublas::norm_1( F_ ); }
    double normSol( ) const { return ublas::norm_2( F_ ); }
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

    //! Clear matrix storage (i.e., fill with zeros)
    void clearMatrix( ) { K_.erase( ); this->clearSuperMatrices_(); }
    //! Clear RHS storage (i.e., fill with zeros)
    void clearRhs( )    { F_.clear( ); }

    //! Solution accessor
    void giveSolution( const VectorUInt_ & dofIndices, VectorLocal_ & result ) const;
    //@}

    //! deprecated?
    unsigned numNonZeros() const { return K_.nnz(); }

    //! Set estimate condition number
    //!
    //! Note: Must be set prior to inversion with #solveSystem().
    //!
    //! \param[in]  norm  If '#'        then do not estimate condition number,
    //!                   if '1' or 'O' then use 1-norm,
    //!                   if 'I'        then infinity-norm
    void setConditionNumberNorm( const char norm ) { conditionNumberNorm_ = norm; }
    
    //! Get estimate of condition number
    double getConditionNumber() const
    {
        FTL_VERIFY( conditionNumberNorm_ != '#' );  // check if computed
        return conditionNumber_;
    }

protected:
    //! Initialise a SuperLU SuperMatrix such that its data pointer is NULL
    static void initSuperMatrix_( SuperMatrix & superMatrix )
    {
        superMatrix.nrow  = 0;
        superMatrix.ncol  = 0;
        superMatrix.Store = NULL;
    }

    //! Initialise all SuperMatrix'es which are used with SystemSolveLU#systemSolve(bool)
    void initSuperMatrices_()
    {
        initSuperMatrix_( A_ );
        initSuperMatrix_( U_ );
        initSuperMatrix_( L_ );
        initSuperMatrix_( rhs_ );
        initSuperMatrix_( X_ );
    }

    //! Clear all data used with SystemSolveLU#systemSolve(bool)
    void clearSuperMatrices_()
    {
        std::vector<int> colInd_tmp, rowPtr_tmp;
        colInd_.swap( colInd_tmp );
        rowPtr_.swap( rowPtr_tmp );
        std::vector<double> Anz_tmp;
        Anz_.swap( Anz_tmp );

        std::vector<int> perm_c_tmp, perm_r_tmp, etree_tmp;
        perm_c_.swap( perm_c_tmp );
        perm_r_.swap( perm_r_tmp );
        etree_.swap( etree_tmp );

        std::vector<double> R_tmp, C_tmp, sol_tmp;
        R_.swap( R_tmp );
        C_.swap( C_tmp );
        sol_.swap( sol_tmp );

        if ( A_.Store ) { Destroy_SuperMatrix_Store( &A_ ); initSuperMatrix_( A_ ); }
        if ( L_.Store ) { Destroy_SuperNode_Matrix( &L_ );  initSuperMatrix_( L_ ); }
        if ( U_.Store ) { Destroy_CompCol_Matrix( &U_ );    initSuperMatrix_( U_ ); }
    }

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
    const unsigned              numDofs_;  //!< number of dofs = matrix size
    corlib::Triplet<int,double> K_;        //!< global stiffness matrix (stored in a triplet)
    ublas::vector<double>       F_;        //!< global force vector

    char    conditionNumberNorm_;  //!< norm for estimating condition number
    double  conditionNumber_;      //!< condition number estimate of matrix

private:
    //! @name Variables for reusing factored matrix call #systemSolve(bool)
    //@{
    SuperMatrix          A_;          //!< sparse matrix on LHS
    SuperMatrix          L_;          //!< lower triangle of LU decomposition
    SuperMatrix          U_;          //!< upper triangle of LU decomposition
    SuperMatrix          rhs_;        //!< right-hand-side
    SuperMatrix          X_;          //!< solution vector

    std::set<unsigned>   oldConstrDofs_;  //!< previously constrained dofs
    std::set<unsigned>   newConstrDofs_;  //!< currently constrained dofs

    std::vector<double>  Anz_;        //!< non-zeros in #A_
    std::vector<int>     colInd_;     //!< column indices in of content in #A_
    std::vector<int>     rowPtr_;     //!< row offsets of content in #A_

    std::vector<int>     perm_c_;     //!< column permutation vector
    std::vector<int>     perm_r_;     //!< row permutation vector

    std::vector<int>     etree_;      //!< Elimiation tree
    std::vector<double>  sol_;        //!< Data vector for secondary solution #X_
    std::vector<double>  R_;          //!< row scale factors
    std::vector<double>  C_;          //!< column scale factors
    
    char                 equed_;      //!< Form of equilibrium done
    //@}
};

//------------------------------------------------------------------------------
//! insert element stiffness matrix into global system matrix
void corlib::SystemSolveLU::insertToMatrix( const MatrixLocal_ & eStiff, 
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
void corlib::SystemSolveLU::insertToRhs( const VectorLocal_ & eForce, 
                                         const VectorUInt_ & dofIndices )
{
    for ( unsigned i = 0; i < dofIndices.size( ); i ++ ) {
        F_[ dofIndices[ i ] ] += eForce( i );
    }
    return;
}

//------------------------------------------------------------------------------
//! apply constraints by a sort of a penalty method
void corlib::SystemSolveLU::applyConstraints( VectorConstraint_ & constraints,
                                              const double   factor, 
                                              const double   scalar )
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
void corlib::SystemSolveLU::fixDOF( const unsigned index, const double scalar )
{
    K_.zeroRow( index, scalar );
    F_[ index ] = 0.;
    newConstrDofs_.insert( index );
    return;
}

//------------------------------------------------------------------------------
//! solve the system A x = b directly by means of an LU decomposition
void corlib::SystemSolveLU::solveSystem( )
{
    // set up superlu matrix
    SuperMatrix A, L, U, rhs;

    std::vector<double> Anz;
    std::vector<int> colInd, rowPtr;

    // let triplet object fill the arrays
    K_.setCompressedRowArrays( colInd, rowPtr, Anz );

    dCreate_CompRow_Matrix( &A, static_cast<int> (numDofs_),
    						static_cast<int> (numDofs_),
    						static_cast<int> (Anz.size()),
                            &(Anz[0]), &(colInd[0]), &(rowPtr[0]), 
                            SLU_NR, SLU_D, SLU_GE );

    //--------------------------------------------------------------------------
    // set up rhs
    const int numRhs = 1;
    dCreate_Dense_Matrix( &rhs, numDofs_, numRhs, &(F_[0]), numDofs_, 
                          SLU_DN, SLU_D, SLU_GE );

    //--------------------------------------------------------------------------
    // prepare solver
    superlu_options_t options;
    set_default_options( &options );
    //options.DiagPivotThresh = 0.0;
    //options.PrintStat = YES;

    SuperLUStat_t stat;
    StatInit( &stat );

    std::vector<int> perm_c( numDofs_, 0 );
    std::vector<int> perm_r( numDofs_, 0 );
    int info = 0;

    //--------------------------------------------------------------------------
    // solver
    dgssv( &options, &A, &(perm_c[0]), &(perm_r[0]), &L, &U, &rhs, &stat, &info );
    
    FTL_VERIFY_DESCRIPTIVE( (info==0), "SuperLU failed to solve with info=%d\n", info );

    //StatPrint( &stat );

    // estimate condition number
    if ( conditionNumberNorm_ != '#' ) {
        char norm    = conditionNumberNorm_;
        double anorm = dlangs( &norm, &A );  // regarding this function look at top of file
        double rcond = 0.;
        dgscon( &norm, &L, &U, anorm, &rcond, &stat, &info );
        conditionNumber_ = 1. / rcond;
    }

    //--------------------------------------------------------------------------
    // cleanup
    Destroy_SuperMatrix_Store( &A );
    Destroy_SuperMatrix_Store( &rhs ); 
    Destroy_SuperNode_Matrix( &L );
    Destroy_CompCol_Matrix( &U );
    StatFree( &stat );

    oldConstrDofs_.clear();
    newConstrDofs_.clear();

    return;
}

//------------------------------------------------------------------------------
//! solve the system A x = b directly by means of an LU decomposition
void corlib::SystemSolveLU::solveSystem( const bool reuseFactoredMatrix )
{
    // fill matrices
    if ( not reuseFactoredMatrix ) {
        // let triplet object fill the arrays
        K_.setCompressedRowArrays( colInd_, rowPtr_, Anz_ );

        dCreate_CompRow_Matrix( &A_, static_cast<int> (numDofs_),
        						static_cast<int> (numDofs_),
        						static_cast<int> (Anz_.size()),
                                &(Anz_[0]), &(colInd_[0]), &(rowPtr_[0]), 
                                SLU_NR, SLU_D, SLU_GE );
    }

    //--------------------------------------------------------------------------
    // check constraint DOFs
    const bool hasSameConstr = this -> updateConstraintDofs_();
    if ( reuseFactoredMatrix ) FTL_VERIFY( hasSameConstr );

    //--------------------------------------------------------------------------
    // set up RHS : link #rhs_ to #F_ ==> solution will turn up on F_
    const int numRhs = 1;
    dCreate_Dense_Matrix( &rhs_, numDofs_, numRhs, &(F_[0]), numDofs_, 
                          SLU_DN, SLU_D, SLU_GE );

    //--------------------------------------------------------------------------
    // prepare solver
    superlu_options_t options;
    set_default_options( &options );
    //options.DiagPivotThresh = 0.0;
    //options.PrintStat = YES;
    options.Fact = ( reuseFactoredMatrix ) ? FACTORED : DOFACT;

    // allocate a few intermediate vectors
    if ( not reuseFactoredMatrix ) {
        perm_c_.resize( numDofs_, 0 );
        perm_r_.resize( numDofs_, 0 );
        etree_.resize( numDofs_ );
        R_.resize( numDofs_ );
        C_.resize( numDofs_ );
    }
    sol_.resize( numDofs_, 0. );
    dCreate_Dense_Matrix( &X_, numDofs_, numRhs, &(sol_[0]), numDofs_, 
                          SLU_DN, SLU_D, SLU_GE );
    
    // more intermediates for solver
    void *        work  = NULL;
    int           lwork = 0;   // implies internal allocation
    double        rpg   = 0.;
    double        rcond = 0.;  // computation depends options
    double        ferr;        // for certain options this needs to be a vector
    double        berr;        // for certain options this needs to be a vector
    mem_usage_t   mem_usage;
    int           info  = 0;
    SuperLUStat_t stat;
    StatInit( &stat );

    //--------------------------------------------------------------------------
    // solver
#if (OPENFTL_SUPERLU_VERSION == 5)
    GlobalLU_t Glu;
    dgssvx( &options, &A_, &(perm_c_[0]), &(perm_r_[0]), &(etree_[0]), &equed_,
            &(R_[0]), &(C_[0]), &L_, &U_, work, lwork, &rhs_, &X_,
            &rpg, &rcond, &ferr, &berr, &Glu, &mem_usage, &stat, &info );
#else
    dgssvx( &options, &A_, &(perm_c_[0]), &(perm_r_[0]), &(etree_[0]), &equed_,
               &(R_[0]), &(C_[0]), &L_, &U_, work, lwork, &rhs_, &X_,
               &rpg, &rcond, &ferr, &berr, &mem_usage, &stat, &info );
#endif
    
    FTL_VERIFY_DESCRIPTIVE( (info==0), "SuperLU failed to solve with info=%d\n", info );
    std::copy( sol_.begin(), sol_.end(), F_.begin() );

    //StatPrint( &stat );

    // estimate condition number
    if ( conditionNumberNorm_ != '#' ) {
        char norm    = conditionNumberNorm_;
        double anorm = dlangs( &norm, &A_ );  // regarding this function look at top of file
        double rcond = 0.;
        dgscon( &norm, &L_, &U_, anorm, &rcond, &stat, &info );
        conditionNumber_ = 1. / rcond;
    }

    //--------------------------------------------------------------------------
    // cleanup --- clean-up of matrices is called along #clearMatrix()
    std::vector<double>().swap( sol_ );
    Destroy_SuperMatrix_Store( &rhs_ );
    initSuperMatrix_( rhs_ );
    Destroy_SuperMatrix_Store( &X_ );
    initSuperMatrix_( X_ );
    StatFree( &stat );

    return;
}


//------------------------------------------------------------------------------
//! return solution values, access with vectors
void corlib::SystemSolveLU::giveSolution( const VectorUInt_ & dofIndices, 
                                          VectorLocal_      & result ) const
{
    for ( unsigned i = 0; i < dofIndices.size( ); i ++ ) {
        result[ i ] = F_[ dofIndices[i] ];
    }
    return;
}

#endif
