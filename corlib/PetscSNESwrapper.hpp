// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//
//                 Computational Structural Mechanics Lab
//                         University of Cambridge
//
// This software is copyrighted and all rights are retained by the University of
// Cambridge. All use, disclosure, and/or reproduction of any part not expressly
// authorised by the University of Cambridge is prohibited. (C) 2010
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

//! @file PetscSNESwrapper.hpp

#ifndef corlib_petscsneswrapper_h
#define corlib_petscsneswrapper_h
//------------------------------------------------------------------------------
#include "petscsnes.h"
#include <boost/numeric/ublas/matrix_sparse.hpp>

//------------------------------------------------------------------------------
namespace corlib{

    template<typename RESIDUAL> class PetscSNESwrapper;

    namespace ublas = boost::numeric::ublas;
}
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
template<typename RESIDUAL>
class corlib::PetscSNESwrapper
{
    typedef ublas::vector<double>                     UblasVec_;
    typedef std::vector< std::pair<unsigned,double> > ConstraintVec_;

public:
    //! Constructor with residual  pointer
    PetscSNESwrapper( RESIDUAL * residual, const unsigned size, 
                      const double constraintFactor );

private:
    //! Get the constraints
    void getConstraints(){ return residual_ -> fillConstraintVec( constraints_ ); }

    //! Distribute the solution
    void distributeSolution( UblasVec_ & ublasX );

    //! Insert constraint-values to increment
    void applyConstraintsToIncrement( UblasVec_ & ublasX );

    //! Compute new residual value F(X)
    void computeResidual( ) { return residual_ -> computeNodalResiduals( ); }

    //! Collect the residual
    void collectResidual( UblasVec_ & ublasF );

    //! Zero the values corresponding to constrained DOFs
    void applyConstraintsToResidual( UblasVec_ & ublasF );

public:
    //! Set an initial guess (here only zeros)
    void initialize() { VecSet( X_, 0. ); return; }

    //! The actual solve operation
    void solve() { SNESSolve( solver_, PETSC_NULL, X_ ); return; }

    //! Clear the allocated data
    PetscErrorCode clear();

    //! Function pointer wrapper
    static PetscErrorCode formFunction( SNES snes, Vec x, Vec f, void * ctx );

protected:
    RESIDUAL        *  residual_;         //!< Pointer to the residual object
    int            size_;                 //!< Number of degrees of freedom
    MPI_Comm       comm_;                 //!< MPI communicator
    SNES           solver_;               //!< Solver object
    Vec            F_;                    //!< Residual vector  F
    Vec            X_;                    //!< Solution vector X
    Mat            J_;                    //!< (Matrix-free) Jacobian matrix J=F'
    PetscErrorCode ierr_;                 //!< Error code
    double         constraintFactor_;     //!< multiplier for the constraints
    ConstraintVec_ constraints_;          //!< vector of <dofId,value>-pairs
};

//------------------------------------------------------------------------------
//! Constructor with residual  pointer
template<typename RESIDUAL>
corlib::PetscSNESwrapper<RESIDUAL>::PetscSNESwrapper( RESIDUAL * residual, 
                                                   const unsigned size,
                                                   const double constraintFactor )
    : residual_( residual ), size_( size ), 
      constraintFactor_( constraintFactor )
{
    // MPI communicator
    comm_ = PETSC_COMM_WORLD;
        
    // Create a solver object
    SNESCreate( comm_, &solver_ );
        
    // Create vectors X and F
    VecCreate( comm_, &X_ );
    VecSetSizes( X_, PETSC_DECIDE, size_ );
    VecSetFromOptions( X_ );
    VecDuplicate( X_, &F_ );

    // Insert the function for computing F(X)
    SNESSetFunction( solver_, F_, formFunction, this );

    // Create a matrix-free matrix J
    MatCreateSNESMF( solver_, &J_ );

    // Check for the -snes_mf flag given by caller
    PetscTruth flag;
    PetscOptionsHasName( PETSC_NULL, "-snes_mf", &flag );
    if ( !flag ) {
        std::cout << "Needs to be called with: -snes_mf  --> exit" << std::endl;
        exit( 1 );
    }
        
    // Set solver options
    SNESSetFromOptions( solver_ );

    // obtain constraints
    this -> getConstraints();
}

//------------------------------------------------------------------------------
/** Apply constraints to solution, then distribute it
 *  \param[in] ublasX ublas-vector containing the current increment
 */
template<typename RESIDUAL>
void corlib::PetscSNESwrapper<RESIDUAL>::distributeSolution( UblasVec_ & ublasX )
{ 
    this -> applyConstraintsToIncrement( ublasX );
    residual_ -> distributeSolution(  ublasX ); 
    return;
}

//-----------------------------------------------------------------------------
/** Apply constraints to the current solution: the values of constrained dofs
 *  are inserted into the solution vector. 
 *  \param[in,out] ublasX ublas-vector containing the current increment
*/
template<typename RESIDUAL>
void corlib::PetscSNESwrapper<RESIDUAL>::applyConstraintsToIncrement( UblasVec_ & ublasX )
{
    ConstraintVec_::iterator cstr = constraints_.begin();
    ConstraintVec_::iterator cend = constraints_.end();
    for (; cstr != cend; ++cstr ) {
        ublasX[ cstr->first ] = (cstr-> second) * constraintFactor_;
    }
    return;
}

//------------------------------------------------------------------------------
/** Collect the residual forces from the Residual object and apply constraints 
 *  to them.
 *  \param[out] ublasF ublas-Vector containing the residual forces
 */
template<typename RESIDUAL>
void corlib::PetscSNESwrapper<RESIDUAL>::collectResidual( UblasVec_ & ublasF ) 
{ 
    residual_ -> collectResidual(    ublasF ); 
    this -> applyConstraintsToResidual( ublasF );
    return;
}

//------------------------------------------------------------------------------
/** Set the entries corresponding to constrained DOFs to zero.
 *  \param[in,out] ublasF  ublas-vector containing the residual forces
 */
template<typename RESIDUAL>
void corlib::PetscSNESwrapper<RESIDUAL>::applyConstraintsToResidual( UblasVec_ & ublasF )
{
    ConstraintVec_::iterator cstr = constraints_.begin();
    ConstraintVec_::iterator cend = constraints_.end();
    for (; cstr != cend; ++cstr ) ublasF[ cstr->first ] = 0.;
    return;
}

//------------------------------------------------------------------------------
//! Clear all data which was allocated by Petsc
template<typename RESIDUAL>
PetscErrorCode corlib::PetscSNESwrapper<RESIDUAL>::clear()
{
    residual_ = NULL;
    ierr_ = SNESDestroy( solver_ ); CHKERRQ( ierr_ );
    ierr_ = VecDestroy( F_ );       CHKERRQ( ierr_ );
    ierr_ = VecDestroy( X_ );       CHKERRQ( ierr_ );
    ierr_ = MatDestroy( J_ );       CHKERRQ( ierr_ );
    return ierr_;
}

//------------------------------------------------------------------------------
/** Function whose pointer is given as a call-back to Petsc. The signature of
 *  this function is determined by Petsc and the SNES wrapper is hidden in the
 *  void-pointer.
 *  \param      snes The solver object for SNES
 *  \param[in]  x    Vector with current solution
 *  \param[out] f    Function f(x), i.e., the residual
 *  \param      ctx  PetscSNESwrapper hidden by the void-pointer
 */
template<typename RESIDUAL>
PetscErrorCode corlib::PetscSNESwrapper<RESIDUAL>::formFunction( SNES snes, Vec x, 
                                                              Vec f, void * ctx )
{
    // unmask the void-pointer
    PetscSNESwrapper * mySelf = (PetscSNESwrapper *) ctx;

    // get arrays and size from Petsc vectors
    PetscScalar *xx, *ff;
    int  dim;
    VecGetArray( x, &xx );
    VecGetArray( f, &ff );
    VecGetSize( x, &dim );

    // set up a ublas vector for X
    ublas::vector<double> ublasX( dim );
    std::copy( xx, xx+dim, ublasX.begin() );

    // scatter the increment to the nodes
    mySelf -> distributeSolution( ublasX );

    // compute F(X)
    mySelf -> computeResidual( );

    // set up a ublas vector for F
    ublas::vector<double> ublasF( dim );

    // collect the computed residual
    mySelf -> collectResidual( ublasF );
        
    // pass back the solution to Petsc object
    std::copy( ublasF.begin(), ublasF.end(), ff );

    // restore the Petsc vectors
    VecRestoreArray( x, &xx);
    VecRestoreArray( f, &ff);

    return 0;
}
//------------------------------------------------------------------------------
#endif
