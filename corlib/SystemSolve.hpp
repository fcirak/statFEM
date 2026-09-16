// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//
//                 Computational Structural Mechanics Lab
//                         University of Cambridge
//
// This software is copyrighted and all rights are retained by the University of
// Cambridge. All use, disclosure, and/or reproduction of any part not expressly
// authorised by the University of Cambridge is prohibited. (C) 2010
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

//! @file SystemSolve.hpp

#ifndef corlib_systemsolve_h
#define corlib_systemsolve_h

#include <ostream>
#include <vector>
#include <map>
#include <cmath>

#include <Eigen/Core>
#include <Eigen/Dense> //For dense LU operation

#include <corlib/linalg.hpp>
#include <corlib/verify.hpp>
#include <corlib/eigenX.hpp>

//------------------------------------------------------------------------------
namespace corlib{
    class SystemSolve;
}

//------------------------------------------------------------------------------
/** \brief System matrix handler and solver based on Eigen dense matrices
 */
class corlib::SystemSolve
{
public:
    //! type of solver to use
    enum solver {
        LU,        //!< uBLAS build-in LU decomposition
        CHOLESKY   //!< Cholesky decomposition
    };

private:
    typedef Eigen::MatrixXd                                 ElemStiff_;
    typedef Eigen::VectorXd                                 ElemForce_;
                                                            
    typedef std::vector< unsigned >                         DofVec_;
    typedef std::vector< double >                           ValueVec_;
    typedef std::vector< std::pair<unsigned,double> >       ConstraintVec_;
                                                            
    typedef Eigen::VectorXd                                 Vec_;
    typedef Eigen::MatrixXd                                 Mat_;

public:
    SystemSolve( const unsigned & numDof ) 
    {
        solverType_ = LU;
        K_.resize( numDof, numDof );
        K_.setZero();
        F_.resize( numDof );
        F_.setZero();
    }

    void insertToMatrix( const ElemStiff_ & eStiff, 
                         const DofVec_ & rowIndices,
                         const DofVec_ & colIndices );

    void insertToRhs( const ElemForce_ & eForce, const DofVec_ & dofIndices );

    void writeMatrix( std::ostream & out ) { corlib::writeCSVmatrix( out, K_ ); }
    void writeRhs(    std::ostream & out ) { corlib::writeVector( out, F_ ); }

    double normRhs( ) const { return F_.norm( ); }
    double normSol( ) const { return F_.norm( ); }

    double normMatrix( ) const { return K_.norm( ); }

    void applyConstraints( ConstraintVec_ & constraints, const double factor = 1., 
                           const double scalar = 1. );
    void finishAssembly( ) { } //to be compliant with Petsc sysmat

    //! Set system solver type
    //!
    //! \param[in]  sType  flag solver for solver type
    void setType( const enum solver sType ) { solverType_ = sType; return; }

    //! Solve the system
    void solveSystem( );

    void clearMatrix( ) { K_.resize( 0, 0 ); }
    void clearRhs( )  {   F_.resize( 0 ); }

    void giveSolution( const DofVec_ & dofIndices, ElemForce_ & result ) const;

private:
    enum solver solverType_ ;  //!< solver type
    Mat_ K_;      //!< global stiffness matrix
    Vec_ F_;      //!< global force vector
};

//------------------------------------------------------------------------------
//! insert element stiffness matrix into global system matrix
void corlib::SystemSolve::insertToMatrix( const ElemStiff_ & eStiff, 
                                       const DofVec_ & rowIndices,
                                       const DofVec_ & colIndices )
{
     for ( unsigned i = 0; i < rowIndices.size( ); i ++ ) {
        for ( unsigned j = 0; j < colIndices.size( ); j ++ ) {
            K_( rowIndices[ i ], colIndices[ j ] ) += eStiff( i, j );
        }
    }

    return;
}

//------------------------------------------------------------------------------
//! insert element stiffness matrix into global system matrix
void corlib::SystemSolve::insertToRhs( const ElemForce_ & eForce, 
                                    const DofVec_ & dofIndices )
{
    for ( unsigned i = 0; i < dofIndices.size( ); i ++ ) {
        F_( dofIndices[ i ] ) += eForce( i );
    }
    return;
}

//------------------------------------------------------------------------------
//! apply constraints by a sort of a penalty method
void corlib::SystemSolve::applyConstraints( ConstraintVec_ & constraints,
                                            const double factor,
                                            const double scalar )
{

    ConstraintVec_::const_iterator cIter = constraints.begin( );
    ConstraintVec_::const_iterator cEnd  = constraints.end( );

    for ( ; cIter != cEnd; ++ cIter ) {
        const unsigned index = cIter -> first;
        const double value   = factor * (cIter -> second);

        // main diagonal element
        const double diagEntry = K_( index, index );

        for ( unsigned colIndex = 0; colIndex < K_.cols(); ++ colIndex ) {
            // put column of prescribed dof to RHS
            F_( colIndex ) -= value * (K_( index, colIndex ) );
            K_( index, colIndex ) = 0.0;
            K_( colIndex, index ) = 0.0;
        }

        // new diagonal element
        const double newDiagEntry = diagEntry * scalar;
        // set new diagonal element and multiply the RHS entry with it
        K_( index, index ) = newDiagEntry;
        F_( index ) = value * newDiagEntry;

    }

    return;
}

//------------------------------------------------------------------------------
//! solve the system A x = b directly by means of a Cholesky decomposition
void corlib::SystemSolve::solveSystem( )
{
    switch ( solverType_ ) {
    case CHOLESKY : {
        Eigen::LLT< Mat_ > llt( K_ );
        Eigen::VectorXd x;
        x = llt.solve(F_);

        if ( ( K_ * x ).isApprox( F_ ) ) {
            F_ = x;
        }
        else {
            FTL_VERIFY_DESCRIPTIVE( true, "Cholesky decomposition failed\n" );
        }
    }
    break;
    case LU : {
        // create a permutation matrix for the LU-factorization
    	Eigen::FullPivLU< Mat_ > lu( K_ );
    	Eigen::VectorXd x;
    	x = lu.solve(F_);

    	if ( ( K_ * x ).isApprox( F_ ) ) {
    		F_ = x;
		}
		else {
			FTL_VERIFY_DESCRIPTIVE( true, "LU decomposition failed\n" );
		}
    }
    break;
    default : {
        FTL_VERIFY_DESCRIPTIVE( false, "Solver is not available\n" );
    }
    break;
    }

    return;
}


//------------------------------------------------------------------------------
//! return solution values 
void corlib::SystemSolve::giveSolution( const DofVec_ & dofIndices, 
                                        ElemForce_    & result ) const
{
    for ( unsigned i = 0; i < dofIndices.size( ); i ++ ) {
        result[ i ] = F_( dofIndices[i] );
    }
    return;
}

#endif
