// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//
//                 Computational Structural Mechanics Lab
//                         University of Cambridge
//
// This software is copyrighted and all rights are retained by the University of
// Cambridge. All use, disclosure, and/or reproduction of any part not expressly
// authorised by the University of Cambridge is prohibited. (C) 2010
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

//! @file SystemSolvePetsc.hpp

#ifndef corlib_systemsolvepetsc_h
#define corlib_systemsolvepetsc_h

//------------------------------------------------------------------------------
//! system includes
#include <ostream>
#include <vector>
#include <set>
#include <map>
#include <cmath>
//! petsc includes
#include <petscpc.h>
#include <petscksp.h>
//! boost includes
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>
//! corlib includes
#include <corlib/Triplet.hpp>

// verify PETSc version
#if ( PETSC_VERSION_MAJOR != 3 || ( PETSC_VERSION_MINOR != 1 && PETSC_VERSION_MINOR != 2 && PETSC_VERSION_MINOR != 4  ) )
 # error OpenFTL only supports PETSc versions 3.1, 3.2, 3.4. It looks like you have some other version.
# else 
  #if( PETSC_VERSION_MINOR == 4 )
    #define FTL_PETSC_3_4
  #elif( PETSC_VERSION_MINOR == 2 )
    #define FTL_PETSC_3_2
  #else 
    #define FTL_PETSC_3_1
  #endif
#endif

//------------------------------------------------------------------------------
namespace corlib{

    class SystemSolvePetsc;

    namespace ublas = boost::numeric::ublas;

    namespace detail_{

        //! Cast PetscScalar to NUMBERTYPE
        template< typename NUMBERTYPE >
        NUMBERTYPE castPetscScalar( const PetscScalar & r )
        {
            return static_cast<NUMBERTYPE>( r );
        }


        // Make a pointer to the data array out of a std::vector
        template<typename T>
        T *  makePetscPointer( std::vector<T> & array )
        {
            if ( array.size() ) return &(array[0]);
            return PETSC_NULL;
        }

    }
}

//------------------------------------------------------------------------------
/** \brief PETSc based linear system solver
 *
 *  \details  This class implements solution of a system of linear equations
 *  by PETSc solvers.
 *
 *  <h2>Selecting solution method and preconditioner</h2>
 *  Particular solution method and preconditioner
 *  have to be inserted using the function insertOptions(string).
 *  using the e.g. the flags -ksp_type and -pc_type. For instance 
 *  for an Additive Schwarz preconditioner for GMRES and overlap 10 elements:
 *     -ksp_type gmres -pc_type asm -pc_asm_overlap 10
 *  or for direct solution by MUMPS:
 *     -ksp_type preonly -pc_type lu -pc_factor_mat_solver_package mumps
 *
 *  <h2>Matrix storage</h2>
 *  The class uses preallocation of the PETSc matrix. In the pre-allocation 
 *  step the non-zero entries in the sparse PETSc MPIAIJ matrix are alloacted
 *  before the matrix is filled. The non-zero pattern is deduced from the
 *  #triplet_ intermediate matrix. The triplet is essentially a 3-tuple made of 
 *  (rowDof,colDof,value) --- see Triplet.hpp for details.
 *  (Until revision 2934, the version was without preallocation.)
 *  In any case the distribution pattern of the PETSc matrix, which is based
 *  on consecutive rows of the matrix, is also applied to the right-hand-side
 *  (RHS) vector. The RHS vector is handled like the matrix. It is temporarily
 *  stored in #rhsTwinlet_ and then copied to a RHS PETSc column-vector.
 *
 *  <h2>Distributed assembly</h2>
 *  Distributed assembly is controlled with the #overlappingDD_ flag.
 *  For #overlappingDD_ = true values are INSERTED into the Petsc matrix.
 *  This avoids communication at PETSc assembly stage of the matrix in case 
 *  of distributed meshes with an overlap consisting  of sufficient 
 *  number of layers of elements.
 *  For #overlappingDD_ = false values are ADDED into the Petsc matrix.
 *  This is to be used in case of distributed meshes with no overlap. 
 *
 *  <h2>Calling sequence</h2>
 *  The class is used as follows:
 *    0. insertOptions      - make solver choice (typically set once, used in 5.)
 *
 *    1. insertToMatrix     - called many times, inserts values into #triplet_
 *    2. insertToRhs        - called many times, inserts values into RHS #rhsTwinlet_
 *    3. finishAssembly     - assembles matrix and assembles RHS
 *    4. applyConstraints   - enforces boundary conditions to matrix and RHS
 *    5. solveSystem        - solves the system
 *    6. clearMatrix        - blanks matrix values (optionally keeps structure)
 *    7. clearRhs           - blanks RHS values (optionally keeps structure)
 */
class corlib::SystemSolvePetsc
{
private:
    typedef corlib::Triplet<PetscInt,PetscReal>             Triplet_;
    //! Type for (rowDof,value)-tuples
    typedef std::map<PetscInt,PetscReal>                    Twinlet_;

public:
    typedef ublas::matrix< double >                         SubMat;
    typedef ublas::vector< double >                         SubVec;
    typedef std::vector< unsigned >                         VecUInt;
    typedef std::vector< std::pair<unsigned,double> >       ConstraintVec;

public:
    //! Constructor with global DOF numbers partitioning is decided
    //! by Petsc
    //!
    //! \param[in]  numTotalDofs     Global number of DOFs
    SystemSolvePetsc( const unsigned numTotalDofs )
      : numTotalDofs_( static_cast<PetscInt>( numTotalDofs ) ),
        numLocalDofs_( -1 ),
        overlappingDD_ ( false ),
        comm_( PETSC_COMM_WORLD )
    {
        this -> initialise_();
    }

    //! Constructor already knowing the DOF partitioning
    //!
    //! \param[in]  numTotalDofs     Global number of DOFs
    //! \param[in]  numLocalDofs     Local number of DOFs of process
    //! \param[in]  filterLocalRows  Optional filter all row contributions
    //!                                of DOF not held on process
    //!                                (false is the default )
    //! \param[in]  comm             MPI communicator
    //!                                (PETSC_COMM_WORLD is the default )
    SystemSolvePetsc( const unsigned numTotalDofs,
                       const unsigned numLocalDofs,
                       bool filterLocalRows = false,
                       const MPI_Comm comm = PETSC_COMM_WORLD )
      : numTotalDofs_( static_cast<PetscInt>( numTotalDofs ) ),
        numLocalDofs_( static_cast<PetscInt>( numLocalDofs ) ),
        overlappingDD_ ( filterLocalRows ),
        comm_( comm )
    {
        this -> initialise_();
    }

protected:
    // produce initial state of object
    void initialise_()
    {
        FTL_VERIFY( numTotalDofs_ > 0 );
        FTL_VERIFY( ( numLocalDofs_ == -1 ) or ( numLocalDofs_ > 0 ) );

        // create matrix
        sysMat_ = NULL;
        this -> createMatrix_();
        // create RHS vector
        rhsVec_ = NULL;
        this -> createRhs_();
        // check
        this -> verifyOwnershipRanges_();

        // Petsc solver items
        solver_        = NULL;
        sol_           = NULL;
        solGathered_   = NULL;
        VSdistToLarge_ = NULL;
        myPetscOptions_.clear();
    }

    //! Return number of processors in communicator
    int getNumProcesses_()
    {
        int nProc = -1;
        MPI_Comm_size( comm_, &nProc );
        return nProc;
    }

    //! Return ID of process in communicator
    int getRank_()
    {
        int rank = -1;
        MPI_Comm_rank( comm_, &rank );
        return rank;
    }

    //! Create Petsc RHS vector
    void createRhs_()
    {
        if ( numLocalDofs_ >= 0 )
            ierr_ = VecCreateMPI( comm_, numLocalDofs_, numTotalDofs_, &rhsVec_ );
        else
            ierr_ = VecCreateMPI( comm_, PETSC_DECIDE,  numTotalDofs_, &rhsVec_ );
        CHKERRV( ierr_ );
    }

    //! Create Petsc matrix
    void createMatrix_()
    { 
        ierr_ = MatCreate( comm_, &sysMat_ ); CHKERRV( ierr_ );
        ierr_ = MatSetType( sysMat_, MATMPIAIJ ); CHKERRV( ierr_ );
        if ( numLocalDofs_ >= 0 )
            ierr_ = MatSetSizes( sysMat_, numLocalDofs_, numLocalDofs_,
                                 numTotalDofs_, numTotalDofs_ );
        else
            ierr_ = MatSetSizes( sysMat_, PETSC_DECIDE, PETSC_DECIDE,
                                 numTotalDofs_, numTotalDofs_ );
        ierr_ = MatSetUp( sysMat_ );
        CHKERRV( ierr_ );
    }

    //! Create Petsc solution
    void createSol_()
    {
        FTL_VERIFY_DESCRIPTIVE( rhsVec_ != NULL,    "Right hand side vector not ready." );
        ierr_ = VecDuplicate( rhsVec_, &sol_ ); CHKERRV( ierr_ );
    }

    //! Check if RHS vector and LHS matrix are distributed in the same way
    void verifyOwnershipRanges_()
    {
        // get index range [rhsLow,rhsUpp)
        PetscInt rhsLow = -1, rhsUpp = -1;
        ierr_ = VecGetOwnershipRange( rhsVec_, &rhsLow, &rhsUpp ); CHKERRV( ierr_ );
        FTL_VERIFY( rhsLow >= 0 );
        FTL_VERIFY( rhsUpp > rhsLow );

        // get index range [matLow,matUpp)
        PetscInt matLow = -1, matUpp = -1;
        ierr_ = MatGetOwnershipRange( sysMat_, &matLow, &matUpp ); CHKERRV( ierr_ );
        FTL_VERIFY( matLow >= 0 );
        FTL_VERIFY( matUpp > matLow );

        // check equality of ranges (row maps)
        FTL_VERIFY( matLow == rhsLow );
        FTL_VERIFY( matUpp == rhsUpp );
    }

    //! Erase matrix completely from memory
    void destroyMatrix_()
    {
        FTL_VERIFY( sysMat_ != NULL );
#ifdef FTL_PETSC_3_1
        ierr_ = MatDestroy( sysMat_ ); CHKERRV( ierr_ );
#else
        ierr_ = MatDestroy( &sysMat_ ); CHKERRV( ierr_ );
#endif
        sysMat_ = NULL;
    }
    //! Erase RHS completely from memory
    void destroyRhs_()
    {
        FTL_VERIFY( rhsVec_ != NULL );
#ifdef FTL_PETSC_3_1
        ierr_ = VecDestroy( rhsVec_ ); CHKERRV( ierr_ );
#else
        ierr_ = VecDestroy( &rhsVec_ ); CHKERRV( ierr_ );
#endif
        rhsVec_ = NULL;
    }

public:
    //! Destructor frees the Petsc structures
    ~SystemSolvePetsc()
    {
        if ( sysMat_ ) this -> destroyMatrix_();
        if ( rhsVec_ ) this -> destroyRhs_();
        if ( solver_ != NULL ) {
#ifdef FTL_PETSC_3_1
            ierr_ = KSPDestroy( solver_ ); CHKERRV( ierr_ );
#else
            ierr_ = KSPDestroy( &solver_ ); CHKERRV( ierr_ );
#endif
        }
        if ( sol_ != NULL ) {
#ifdef FTL_PETSC_3_1
            ierr_ = VecDestroy( sol_ ); CHKERRV( ierr_ );
#else
            ierr_ = VecDestroy( &sol_ ); CHKERRV( ierr_ );
#endif
        }
        if ( solGathered_ != NULL ) {
#ifdef FTL_PETSC_3_1
            ierr_ = VecDestroy( solGathered_ ); CHKERRV( ierr_ );
#else
            ierr_ = VecDestroy( &solGathered_ ); CHKERRV( ierr_ );
#endif
        }
        if ( VSdistToLarge_ != NULL ) {
#ifdef FTL_PETSC_3_1
            ierr_ = VecScatterDestroy( VSdistToLarge_ ); CHKERRV( ierr_ );
#else
            ierr_ = VecScatterDestroy( &VSdistToLarge_ ); CHKERRV( ierr_ );
#endif
        }

    }

    //--------------------------------------------------------------------------
    //! @name Insert operations
    //@{

    //! insert submatrix to system matrix - into triplet object
    void insertToMatrix( const SubMat & eStiff,
                         const VecUInt & rowIndices,
                         const VecUInt & colIndices );

    //! insert subvector to system rhs vector
    void insertToRhs( const SubVec & eForce, const VecUInt & dofIndices );
    //@}

    //--------------------------------------------------------------------------
    //! @name Finish assembly calls
    //@{

    //! Finalize assembly of matrix
    //! 
    //! Finish assembly of triplet, preallocate PETSc matrix,
    //! and copy #triplet_ into it.
    void finishMatAssembly();

    //! Finalize assembly of RHS vector
    //!
    //! Copy RHS vector from #rhsTwinlet_ to distributed PETSc #rhsVec_
    void finishRhsAssembly();

protected:
    //! Test wether matrix was assembled
    bool isMatrixAssembled_();

    //! Test wether RHS vector was assembled
    bool isRhsAssembled_();

public:
    //! Calls finishMatAssembly and finishRhsAssembly
    void finishAssembly();

    //@}

    //--------------------------------------------------------------------------
    //! @name Debug methods to write matrix and RHS someplace
    //@{
    void writeMatrix(   std::ostream & out ) { MatView( sysMat_, PETSC_VIEWER_STDOUT_WORLD ); }
    void writeRhs(      std::ostream & out ) { VecView( rhsVec_, PETSC_VIEWER_STDOUT_WORLD ); }
    void writeSolution( std::ostream & out ) { VecView( sol_, PETSC_VIEWER_STDOUT_WORLD ); }
    void writeMatrixMatlab( const std::string filename ) {
        PetscViewer myViewer;
        PetscViewerASCIIOpen( comm_, filename.c_str(), &myViewer );
        PetscViewerSetFormat( myViewer, PETSC_VIEWER_ASCII_MATLAB );
        MatView( sysMat_, myViewer );
#ifdef FTL_PETSC_3_1
        PetscViewerDestroy( myViewer );
#else
        PetscViewerDestroy( &myViewer );
#endif
    }
    void writeRhsMatlab( const std::string filename ) {
        PetscViewer myViewer;
        PetscViewerASCIIOpen( comm_, filename.c_str(), &myViewer );
        PetscViewerSetFormat( myViewer, PETSC_VIEWER_ASCII_MATLAB );
        VecView( rhsVec_, myViewer );
#ifdef FTL_PETSC_3_1
        PetscViewerDestroy( myViewer );
#else
        PetscViewerDestroy( &myViewer );
#endif
    }
    void writeSolutionMatlab( const std::string filename ) {
        PetscViewer myViewer;
        PetscViewerASCIIOpen( comm_ ,filename.c_str(), &myViewer );
        PetscViewerSetFormat( myViewer, PETSC_VIEWER_ASCII_MATLAB );
        VecView( sol_, myViewer );
#ifdef FTL_PETSC_3_1
        PetscViewerDestroy( myViewer );
#else
        PetscViewerDestroy( &myViewer );
#endif
    }
    //@}

    //--------------------------------------------------------------------------
    //! Compute norm of rhs vector
    double normRhs()
    {
        FTL_VERIFY( this -> isRhsAssembled_() );
        double result;
        ierr_ = VecNorm( rhsVec_, NORM_2, &result ); CHKERRQ( ierr_ );
        return result;
    }

    //! Compute 1-norm of rhs vector
    double norm1Rhs()
    {
        FTL_VERIFY( this -> isRhsAssembled_() );
        double result;
        ierr_ = VecNorm( rhsVec_, NORM_1, &result ); CHKERRQ( ierr_ );
        return result;
    }

    //--------------------------------------------------------------------------
    //! @name Constraint degrees of freedom (e.g. Dirichlet BC)
    //@{
    //! Get scalar suitable for diagonal while fixing BC
    double getDiagScalar() const { return diagScalar_; }

    //! Apply constraints
    void applyConstraints( ConstraintVec & constraints,
                           const double factor,
                           const double scalar );

    //! Set any DOF to zero
    void fixDOF( const unsigned index, const double scalar = 1. );

protected:
    //! Update the DOFs sets constaining the constraints
    //!
    //! \return True if constraint DOFs did not change
    bool updateConstraintDofs_()
    {
        std::vector<PetscInt> diffDof;
        std::set_symmetric_difference( newConstrDofs_.begin(), newConstrDofs_.end(),
                                       oldConstrDofs_.begin(), oldConstrDofs_.end(),
                                       std::back_inserter( diffDof ) );
        int sameConstraintsLocal = diffDof.empty() ? 1 : 0;
        
        int sameConstraintsTotal = -1;
        MPI_Allreduce( &sameConstraintsLocal, &sameConstraintsTotal, 1, MPI_INT, MPI_MIN, comm_ );
        
        // update constraint dofs sets
        std::swap( newConstrDofs_, oldConstrDofs_ );
        newConstrDofs_.clear();

        return ( sameConstraintsTotal == 1 );
    }
    //@}

    //--------------------------------------------------------------------------
    //! @name Solution
    //@{
public:
    //! Insert particular PETSc options to solver
    void insertOptions( const std::string additionalOptions )
    {
        myPetscOptions_.append( " " + additionalOptions + " " );
    }

    //! Solve the system
    void solveSystem( const bool reuseOperators = false,
                      const bool removeConstantNullSpace = false,
                      const MatStructure pcMatrixStructure = DIFFERENT_NONZERO_PATTERN );

    //! Compute norm of the solution vector
    double normSol()
    {
        FTL_VERIFY( sol_ );
        double result;
        ierr_ = VecNorm( sol_, NORM_2, &result ); CHKERRQ( ierr_ );
        return result;
    }

    //! Give number of iteration
    unsigned giveNumIterations()
    {
        FTL_VERIFY( solver_ );
        PetscInt numIterations;
        ierr_ = KSPGetIterationNumber( solver_, &numIterations ); CHKERRQ( ierr_ );
        unsigned result = static_cast<unsigned>( numIterations );
        return result;
    }

    //! Give reason of convergence
    int giveConvergedReason()
    {
        FTL_VERIFY( solver_ );
        KSPConvergedReason convergedReason;
        ierr_ = KSPGetConvergedReason( solver_, &convergedReason ); CHKERRQ( ierr_ );
        int result = static_cast<int>( convergedReason );
        return result;
    }

    //! Give preconditioner type
    std::string givePCType()
    {
        FTL_VERIFY( solver_ );
        PC precond;        
        ierr_ = KSPGetPC( solver_, &precond ); //CHKERRV( ierr_ );
#if (defined(FTL_PETSC_3_2) || defined(FTL_PETSC_3_1))
        const PCType pcType;
#else
        PCType pcType;
#endif
        ierr_ = PCGetType( precond, &pcType ); //CHKERRV( ierr_ );
        std::string result = static_cast<std::string>( pcType );
        return result;
    }

    //! Give iterative solver type
    std::string giveKSPType()
    {
        FTL_VERIFY( solver_ );
#if (defined(FTL_PETSC_3_2) || defined(FTL_PETSC_3_1))
        const KSPType kspType;
#else
        KSPType kspType;
#endif
        ierr_ = KSPGetType( solver_, &kspType ); //CHKERRV( ierr_ );
        std::string result = static_cast<std::string>( kspType );
        return result;
    }

    //@}

    //--------------------------------------------------------------------------
    //! Return solution
    void giveSolution( const VecUInt & dofIndices, SubVec & result ) const;

    //--------------------------------------------------------------------------
    //! Fill matrix with zeros
    //!
    //! \param[in] reCreate  This destroys the matrix and creates
    //!                      a new empty one. This should be faster than zeroing.
    void clearMatrix( const bool reCreate = true )
    {
        FTL_VERIFY( sysMat_ != NULL );
        if ( reCreate ) {
            this -> destroyMatrix_();
            this -> createMatrix_();
        }
        else {
            ierr_ = MatZeroEntries( sysMat_ ); CHKERRV( ierr_ );
        }
    }

    //! Fill RHS with zeros
    void clearRhs( const bool reCreate = true )
    {
        FTL_VERIFY( rhsVec_ != NULL );
        if ( reCreate ) {
            this -> destroyRhs_();
            this -> createRhs_();
        }
        else {
            ierr_ = VecSet( rhsVec_, 0. );  CHKERRV( ierr_ );
        }
    }
    //! Fill solution with zeros
    void clearSol()
    {
        FTL_VERIFY( sol_ != NULL );
        ierr_ = VecSet( sol_, 0. );  CHKERRV( ierr_ );
    }

protected:
    const PetscInt  numTotalDofs_;       //!< number of total dofs in the system
    const PetscInt  numLocalDofs_;       //!< number of local dofs, i.e. rows of the MPI matrix at the process
                                         //!< (Note: It is -1 if not set by user at construction)
    const bool  overlappingDD_;          //!< if true, PETSc will only INSERT values local to process

    const MPI_Comm  comm_;               //!< MPI communicator

    Triplet_   triplet_;                 //!< matrix in triplet format
    Mat        sysMat_;                  //!< PETSc system matrix
    double     diagScalar_;              //!< recommended value to be put on diagonal in fixing BC
    Twinlet_   rhsTwinlet_;              //!< RHS vector in form of map:globalDofId->value
    Vec        rhsVec_;                  //!< PETSC vector with RHS values
    
    std::set<PetscInt> oldConstrDofs_;   //!< previously constrained dofs
    std::set<PetscInt> newConstrDofs_;   //!< currently constrained dofs

    std::string myPetscOptions_;         //!< options specifying solver type from applications
    KSP        solver_;                  //!< solver context
    Vec        sol_;                     //!< distributed so lution vector
    Vec        solGathered_;             //!< large local solution vector
    VecScatter VSdistToLarge_;           //!< scatter for gathering local parts into large vector

    std::vector<PetscInt> subdomainDofIndices_; //!< global indices of the local solution vector (to be expected by nodes in local mesh) 
    typedef std::map<unsigned,unsigned> Global2LocalMap_; //! type for storage of global to local map
    Global2LocalMap_ global2LocalDofMap_;     //!< map from global dof indices to subdomain local indices

    PetscErrorCode ierr_;                //!< Petsc's error handler
};

//-----------------------------------------------------------------------------

#include "SystemSolvePetsc.ipp"

//-----------------------------------------------------------------------------

#endif
