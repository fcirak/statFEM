// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//
//                 Computational Structural Mechanics Lab
//                         University of Cambridge
//
// This software is copyrighted and all rights are retained by the University of
// Cambridge. All use, disclosure, and/or reproduction of any part not expressly
// authorised by the University of Cambridge is prohibited. (C) 2010
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

//! @file SystemSolvePetsc.ipp

//------------------------------------------------------------------------------
/** Insert element stiffness matrix into global system matrix
 */
void corlib::SystemSolvePetsc::insertToMatrix( const SubMat  & subMat,
                                               const VecUInt & rowIndices,
                                               const VecUInt & colIndices )
{
    // copy dense matrix into triplet
    for ( unsigned i = 0; i < rowIndices.size(); i++) {
        for ( unsigned j = 0; j < colIndices.size(); j++) {
            // insert matrix into triplet
            triplet_.insert( rowIndices[i], colIndices[j], subMat(i,j) );
        }
    }
    return;
}

//------------------------------------------------------------------------------
/** Insert subvector to the system's RHS vector
 */
void corlib::SystemSolvePetsc::insertToRhs( const SubVec  & subVec,
                                            const VecUInt & dofIndices )
{
    // simply loop over dofs
    const unsigned nDofs = dofIndices.size();
    for ( unsigned iDof = 0; iDof < nDofs; iDof++ ) {

        const PetscInt  index = static_cast<PetscInt>(  dofIndices[iDof] );
        const PetscReal value = static_cast<PetscReal>( subVec[iDof] );

        std::pair<Twinlet_::iterator,bool> rhsIter =
            rhsTwinlet_.insert( Twinlet_::value_type( index, value ) );
        if ( rhsIter.second == false )  // if index already in map, then
            rhsIter.first->second += value;  // update

    }
    return;
}

//------------------------------------------------------------------------------
/** Routine for Petsc's communication.
 *  To be called right after the finish of the MATRIX assembly and before
 *  the RHS assembly start.
 */
void corlib::SystemSolvePetsc::finishMatAssembly()
{
    // get index range [matLow,matHigh)
    PetscInt matLow, matUpp;
    ierr_ = MatGetOwnershipRange( sysMat_, &matLow, &matUpp ); CHKERRV( ierr_ );

    // let Triplet find data for preallocation
    unsigned locBlock = matUpp - matLow ;
    std::vector<PetscInt> nnzRowsIn(  locBlock, 0 );
    std::vector<PetscInt> nnzRowsOut( locBlock, 0 );

    // find number of nonzeros in [ matLow_, matHigh_ ) interval and outside on a row block corresponding
    // to the processor
    triplet_.getIntervalCounts( matLow,  matUpp,
                                nnzRowsIn, nnzRowsOut );

    // check data locality
    const PetscInt sum1 = std::accumulate( nnzRowsIn.begin(), nnzRowsIn.end(),
                                           0, std::plus<PetscInt>() );
    const PetscInt sum2 = std::accumulate( nnzRowsOut.begin(), nnzRowsOut.end(),
                                           0, std::plus<PetscInt>() );

    if ( sum2 > 0.5*sum1 ) {
        std::cout << "Warning: More than 50 percent of entries are nonlocal to process "
                  << this->getRank_() << " : "
                  << float(sum2)/float(sum1) * 100 << " percent" << std::endl;
        std::cout.flush();
    }

#ifdef SYSTEMSOLVEPETSC_PROFILE
    std::cout << "Number of local/non-local entries on rank " << this->getRank_() << " is "
              << float(sum1) << " / " << float(sum2) << ",  "
              << "nonlocal entries in percent: "
              << float(sum2) / float(sum1) * 100
              << std::endl;
    std::cout.flush();
#endif

    // preallocate the MPI matrix
    ierr_ = MatMPIAIJSetPreallocation( sysMat_, 0, &(nnzRowsIn[0]), 0, &(nnzRowsOut[0]) );
    CHKERRV( ierr_ );

    // how many entries to copy?
    unsigned numEntries = triplet_.nnz();

    // structure of the local solution will be deduced from the columns of the local matrix
    // copy entries from triplet to PETSc matrix
    std::set<PetscInt> subdomainDofIndicesSet;
    for ( unsigned i = 0; i < numEntries; i++ ) {

        // extract entry from triplet
        PetscInt row, col;
        PetscReal val;
        triplet_.getEntry( i, row, col, val );

        if ( not overlappingDD_ ) {
            // insert entry to existing system matrix for possible addition
            ierr_ = MatSetValue( sysMat_, row, col, val, ADD_VALUES );
            CHKERRV( ierr_ );
        }
        else if ( ( row >= matLow ) and ( row < matUpp ) ) {
            // if filtering is active, insert final entry to existing system matrix - no addition to be performed
            ierr_ = MatSetValue( sysMat_, row, col, val, INSERT_VALUES );
            CHKERRV( ierr_ );
        }

        subdomainDofIndicesSet.insert( col );
    }
    subdomainDofIndices_.resize( subdomainDofIndicesSet.size() );
    std::copy( subdomainDofIndicesSet.begin(), subdomainDofIndicesSet.end(), subdomainDofIndices_.begin() );

    // prepare auxiliary map for renumbering dofs in the global array
    global2LocalDofMap_.clear( );
    for ( unsigned ind = 0; ind < subdomainDofIndices_.size(); ++ind ) {
        global2LocalDofMap_.insert( std::make_pair( static_cast<unsigned>( subdomainDofIndices_[ind] ), ind ) );
    }

    // clear triplet after it has been copied to PETSc
    diagScalar_ = triplet_.getDiagScalar();
    triplet_.clear(); // this leads to triplet_.nnz()==0

    // do the assembly-begin and -end
#ifdef SYSTEMSOLVEPETSC_PROFILE
    const double time1 = MPI_Wtime();
#endif
    ierr_ = MatAssemblyBegin( sysMat_, MAT_FINAL_ASSEMBLY ); CHKERRV( ierr_ );
    ierr_ = MatAssemblyEnd(   sysMat_, MAT_FINAL_ASSEMBLY ); CHKERRV( ierr_ );
#ifdef SYSTEMSOLVEPETSC_PROFILE
    const double time2 = MPI_Wtime();
#endif

#ifdef SYSTEMSOLVEPETSC_PROFILE
    MatInfo info;
    MatGetInfo( sysMat_, MAT_GLOBAL_SUM, &info);
    const double mal   = info.mallocs;
    const double nz_a  = info.nz_allocated;
    const double nz_u  = info.nz_used;
    const double nz_un = info.nz_unneeded;

    if ( this->getRank_() == 0 ) {
        std::cout << "Time of assembly in Petsc [s]: " << time2 - time1
                  << std::endl << std::flush;
        std::cout << "Number of mallocs: "           << mal   << std::endl
                  << "Number of entries allocated: " << nz_a  << std::endl
                  << "Number of entries used: "      << nz_u  << std::endl
                  << "Number of entries unneeded: "  << nz_un << std::endl
                  << std::flush;
    }
#endif

    return;
}

//------------------------------------------------------------------------------
/** Routine for finishing RHS assembly.
 *  To be called right after the finish of the RHS assembly and before
 *  the constraint application.
 */
void corlib::SystemSolvePetsc::finishRhsAssembly()
{
    PetscInt rhsLow, rhsUpp;
    ierr_ = VecGetOwnershipRange( rhsVec_, &rhsLow, &rhsUpp ); CHKERRV( ierr_ );

    Twinlet_::iterator rhsIter = rhsTwinlet_.begin();
    for ( ; rhsIter != rhsTwinlet_.end(); ++rhsIter ) {

        const PetscInt index = rhsIter->first;
    
        if ( not overlappingDD_ ) {
            ierr_ = VecSetValues( rhsVec_, 1, &index, &(rhsIter->second), ADD_VALUES );
            CHKERRV( ierr_ );
        }
        else if ( ( index >= rhsLow ) and ( index <  rhsUpp ) ) {  // overlappingDD_=true
            ierr_ = VecSetValues( rhsVec_, 1, &index, &(rhsIter->second), INSERT_VALUES );
            CHKERRV( ierr_ );
        }

    }

    // clear RHS data in rhsTwinlet_
    rhsTwinlet_.clear();

    // assemble the vector
    ierr_ = VecAssemblyBegin( rhsVec_ ); CHKERRV( ierr_ );
    ierr_ = VecAssemblyEnd(   rhsVec_ ); CHKERRV( ierr_ );
    return;
}

//------------------------------------------------------------------------------
bool corlib::SystemSolvePetsc::isMatrixAssembled_()
{
    bool isAssembled = true;
    
    isAssembled = isAssembled and ( sysMat_ != NULL );
    FTL_VERIFY( sysMat_ != NULL );

    // check if matrix entries of #triplet_ have been flushed into #sysMat_
    isAssembled = isAssembled and ( triplet_.nnz() == 0 );
    FTL_VERIFY_DESCRIPTIVE( triplet_.nnz() == 0,
                            "You either have not assembled or have inserted"
                            " items after assembly. For your information, the"
                            " number of non-zeros in triplet_ is %d\n",
                            triplet_.nnz() );

    // check if matrix has been assembled
    MatInfo info;
    MatGetInfo( sysMat_, MAT_GLOBAL_SUM, &info );
    //const PetscLogDouble mallocs      = info.mallocs;
    //const PetscLogDouble nz_allocated = info.nz_allocated;
    //const PetscLogDouble nz_used      = info.nz_used;
    //const PetscLogDouble nz_unneeded  = info.nz_unneeded;
    //const PetscLogDouble assemblies   = info.assemblies;
    //isAssembled = isAssembled;// and ( static_cast<int>( nz_used ) > 0 );
    //FTL_VERIFY_DESCRIPTIVE( static_cast<int>( nz_used ) > 0,
    //                        "Petsc matrix appears empty.\n" );
    //isAssembled = isAssembled; //and ( static_cast<int>( assemblies ) > 0 );
    //FTL_VERIFY_DESCRIPTIVE( static_cast<int>( assemblies ) > 0,
    //                        "Petsc matrix has never been assembled.\n" );
    return isAssembled;
}

//------------------------------------------------------------------------------
bool corlib::SystemSolvePetsc::isRhsAssembled_()
{
    bool isAssembled = true;
    
    isAssembled = isAssembled and ( rhsVec_ != NULL );
    FTL_VERIFY( rhsVec_ != NULL );

    // check if matrix entries of #triplet_ have been flushed into #sysMat_
    isAssembled = isAssembled and ( rhsTwinlet_.empty() );
    FTL_VERIFY_DESCRIPTIVE( rhsTwinlet_.empty(),
                            "You either have not assembled or have inserted"
                            " items after assembly into the RHS vector."
                            " For your information, the number of non-zeros"
                            " in rhsTwinlet_ is %d\n",
                            rhsTwinlet_.size() );

    return isAssembled;
}

//------------------------------------------------------------------------------
/** Routine for finishing assemblies of both matrix and RHS.
 *
 *  This method copies the sparse matrix entries held in #triplet_ into
 *  the sparse Petsc matrix #sysMat_.
 */
void corlib::SystemSolvePetsc::finishAssembly()
{
    this -> finishMatAssembly();
    this -> finishRhsAssembly();
    return;
}



//------------------------------------------------------------------------------
/** Apply constraints to the system by replacing row 'i' with vector 'alpha*e_i',
 *  i.e., a zero vector with value 'alpha' at its i-th position, if the dof 'i'
 *  is prescribed. The corresponding value of the prescribed dof is set as the
 *  'i'-th value in the RHS-vector multiplied by a given factor.
 *  Note that Petsc requires the function 'MatZeroRows' to be called by all
 *  processes with all global indices even if not owned by some processes.
 */
void corlib::SystemSolvePetsc::applyConstraints( ConstraintVec & constraints,
                                                 const double factor,
                                                 const double scalar )
{
    FTL_VERIFY( this->isMatrixAssembled_() );
    FTL_VERIFY( this->isRhsAssembled_() );

    // number of constraints
    PetscInt numConstraints = static_cast<PetscInt>( constraints.size() );

    // Additional multiplier for numerical reasons (criterion to be established)
    const PetscScalar diagScalar = scalar;

    std::vector<PetscInt> globalDofs;
    std::vector<PetscScalar>  values;

    // Constraint iterators
    ConstraintVec::const_iterator cIter = constraints.begin();
    ConstraintVec::const_iterator cEnd  = constraints.end();
    // collect global dof indices and the correpsonding values
    for ( ; cIter != cEnd; ++cIter ) {
        PetscInt    row                = cIter -> first;
        PetscScalar prescribedSolValue = static_cast<PetscScalar>( cIter -> second);
        globalDofs.push_back( static_cast<PetscInt>( row ) );
#ifdef FTL_PETSC_3_1
        values.push_back( prescribedSolValue * diagScalar * factor );
#else
        values.push_back( prescribedSolValue * factor );
#endif
    }

#ifdef FTL_PETSC_3_1
    // set matrix rows to zero
    ierr_ = MatZeroRows( sysMat_, numConstraints, &(globalDofs[0]), diagScalar );
    CHKERRV( ierr_ );

    // set RHS entries to values (crashes if called with NULL pointers)
    if ( numConstraints ) {
        ierr_ = VecSetValues( rhsVec_, numConstraints, &(globalDofs[0]), &(values[0]), INSERT_VALUES );
        CHKERRV( ierr_ );
    }

    // perform communication in the rhs vector
    ierr_ = VecAssemblyBegin( rhsVec_ ); CHKERRV( ierr_ );
    ierr_ = VecAssemblyEnd(   rhsVec_ ); CHKERRV( ierr_ );
#else
        if ( sol_ == NULL ) {
        this -> createSol_( );
    }
    if ( numConstraints ) {
        ierr_ = VecSetValues( sol_, numConstraints, &(globalDofs[0]), &(values[0]), INSERT_VALUES );
        CHKERRV( ierr_ );
    }
    // perform communication in the rhs vector
    ierr_ = VecAssemblyBegin( sol_ ); CHKERRV( ierr_ );
    ierr_ = VecAssemblyEnd(   sol_ ); CHKERRV( ierr_ );

    // set matrix rows to zero
    //ierr_ = MatSetOption( sysMat_,MAT_NO_OFF_PROC_ZERO_ROWS,PETSC_TRUE );
    //CHKERRV( ierr_ );
    ierr_ = MatZeroRows( sysMat_, numConstraints, &(globalDofs[0]), diagScalar, sol_, rhsVec_ );
    CHKERRV( ierr_ );
#endif

    // add global DOFs which are constraint
    newConstrDofs_.insert( globalDofs.begin(), globalDofs.end() );

    return;
}

//------------------------------------------------------------------------------
/** Simply set the dof with index 'index' to zero
 *
 * \param[in]  index   Global index of DOF to fix
 * \param[in]  scalar  Value to put on diagonal (default 1.)
 */
void corlib::SystemSolvePetsc::fixDOF( const unsigned index, const double scalar )
{
    std::vector< std::pair<unsigned,double> > thisConstraint;
    thisConstraint.push_back( std::make_pair( index, 0. ) );
    this -> applyConstraints( thisConstraint, 1., scalar );
    thisConstraint.clear();
    return;
}

//------------------------------------------------------------------------------
/** Let Petsc solve the system
 *
 * \param[in]  reuseOperators             Reuse matrices
 * \param[in]  removeConstantNullSpace    Filter out null space (like rigid body modes etc)
 * \param[in]  pcMatrixStructure          Flag indicating information about
 *                                          the preconditioner matrix structure
 *                                          One of SAME_NONZERO_PATTERN,
 *                                                 DIFFERENT_NONZERO_PATTERN (default),
 *                                                 SAME_PRECONDITIONER
 *                                          For more information search for KSPSetOperators()
 */
void corlib::SystemSolvePetsc::solveSystem( const bool reuseOperators,
                                            const bool removeConstantNullSpace,
                                            const MatStructure pcMatrixStructure )
{
    // check if matrix has been assembled
    FTL_VERIFY( this->isMatrixAssembled_() );
    FTL_VERIFY( this->isRhsAssembled_() );

    // if this is the first call to this routine
    if ( solver_ == NULL ) {

        // solver object
        ierr_ = KSPCreate( comm_, &solver_ ); CHKERRV( ierr_ );

        // insert extra options if they were specified
        // Warning: PETSs inserts options globally at this point, so other instances of this SystemSolvePetsc get them as well
        //          to avoid unexpected behaviour, one should either specify desired options for any solver, or specify
        //          options for all instances of the solver
        FTL_VERIFY_DESCRIPTIVE( not myPetscOptions_.empty(),
                                "You must set Petsc solver options\n" );
        ierr_ = PetscOptionsInsertString( myPetscOptions_.c_str() ); CHKERRV( ierr_ );

        // set linear solver by options
        ierr_ = KSPSetFromOptions( solver_ ); CHKERRV( ierr_ );

        // avoid illegal settings
        FTL_VERIFY_DESCRIPTIVE( not reuseOperators,
                                "Nothing is loaded in solver, cannot reuse it. \n");

        // create solution vector
        if ( sol_ == NULL ) {
            this -> createSol_( );
        }

        // create a larger local solution vector for gathering - this vector contains all indices accessed by local elements
        PetscInt numDofsSub = subdomainDofIndices_.size();
        FTL_VERIFY_DESCRIPTIVE( numDofsSub > 0, "Zero size of vector of local indices. Perhaps missing call to loadSolIndices"); 
        FTL_VERIFY_DESCRIPTIVE( numDofsSub > numLocalDofs_, "Mismatch in size"); 
        ierr_ = VecCreateSeq( PETSC_COMM_SELF, numDofsSub, &solGathered_ ); CHKERRV( ierr_ );

        // prepare gathering scatter
        IS mySubIs;
#ifdef FTL_PETSC_3_1
        ierr_ = ISCreateGeneral( comm_, numDofsSub, &(subdomainDofIndices_[0]), &mySubIs); CHKERRV( ierr_ );
#else
        ierr_ = ISCreateGeneral( comm_, numDofsSub, &(subdomainDofIndices_[0]), PETSC_COPY_VALUES, &mySubIs); CHKERRV( ierr_ );
#endif
        ierr_ = VecScatterCreate( sol_, mySubIs, solGathered_, PETSC_NULL, &VSdistToLarge_ ); CHKERRV( ierr_ );
#ifdef FTL_PETSC_3_1
        ierr_ = ISDestroy( mySubIs );
#else
        ierr_ = ISDestroy( &mySubIs );
#endif
    }
    else {
        // clear solution
        this -> clearSol();
    }

    // check constraint DOFs
    const bool hasSameConstr = this -> updateConstraintDofs_();
    if ( reuseOperators ) FTL_VERIFY( hasSameConstr );

    // set operators - this is needed each time matrix changes to correctly update preconditioner
    // SAME_NONZERO_PATTERN may not always be the case - other options are DIFFERENT_NONZERO_PATTERN and SAME_PRECONDITIONER
    if ( not reuseOperators ) {
        ierr_ = KSPSetOperators( solver_, sysMat_, sysMat_, pcMatrixStructure ); CHKERRV( ierr_ );
    }

    // remove one-dimensional nullspace of constants if wanted
    if ( removeConstantNullSpace ) {
        MatNullSpace nullSpace;
        ierr_ = MatNullSpaceCreate( comm_, PETSC_TRUE, 0, PETSC_NULL, &nullSpace ); CHKERRV( ierr_ );
        ierr_ = KSPSetNullSpace( solver_, nullSpace ); CHKERRV( ierr_ );
#ifdef FTL_PETSC_3_1
        MatNullSpaceDestroy( nullSpace );
#else
        MatNullSpaceDestroy( &nullSpace );
#endif
    }

    // solve
    ierr_ = KSPSolve( solver_, rhsVec_, sol_ ); CHKERRV( ierr_ );

    // Compute residual res = A*u - b and store it into rhsVec_
    ierr_ = VecScale( rhsVec_, -1.0 ); CHKERRV( ierr_ );
    ierr_ = MatMultAdd( sysMat_, sol_, rhsVec_, rhsVec_ ); CHKERRV( ierr_ );

    // get global solution vector at each process
    ierr_ = VecScatterBegin( VSdistToLarge_, sol_, solGathered_, INSERT_VALUES, SCATTER_FORWARD ); CHKERRV( ierr_ );
    ierr_ = VecScatterEnd(   VSdistToLarge_, sol_, solGathered_, INSERT_VALUES, SCATTER_FORWARD ); CHKERRV( ierr_ );

    return;
}

//------------------------------------------------------------------------------
void corlib::SystemSolvePetsc::giveSolution( const VecUInt & dofIndices,
                                             SubVec & result ) const
{
    // PETSc's error code
    PetscErrorCode ierr;

    // Number of dof-indices
    PetscInt nDofs = static_cast<PetscInt>( dofIndices.size() );

    std::vector<PetscInt>    effectiveDofs;
    std::vector<PetscScalar> dofSolutions( nDofs );

    VecUInt::const_iterator dofIter = dofIndices.begin();
    VecUInt::const_iterator dofEnd  = dofIndices.end();
    SubVec::iterator valIter = result.begin();
    for ( ; dofIter != dofEnd; dofIter++ ) {

        // map it to local dof
        Global2LocalMap_::const_iterator pos = global2LocalDofMap_.find( *dofIter );
        FTL_VERIFY_DESCRIPTIVE( pos != global2LocalDofMap_.end(),
                                "Cannot remap index %d to local indices in solution distribution. \n ", *dofIter );
        unsigned indLoc = pos -> second;

        effectiveDofs.push_back( indLoc );
    }

    // simply get the dof solutions by using the LOCAL indices
    ierr = VecGetValues( solGathered_, nDofs, &(effectiveDofs[0]), &(dofSolutions[0]) );
    CHKERRV( ierr );

    boost::function<double(const PetscScalar &)> convert
        = boost::bind( &corlib::detail_::castPetscScalar<double>, _1 );
    std::transform( dofSolutions.begin(), dofSolutions.end(),
                    result.begin(), convert );


}

