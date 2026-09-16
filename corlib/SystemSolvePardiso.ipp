// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//
//                 Computational Structural Mechanics Lab
//                         University of Cambridge
//
// This software is copyrighted and all rights are retained by the University of
// Cambridge. All use, disclosure, and/or reproduction of any part not expressly
// authorised by the University of Cambridge is prohibited. (C) 2010
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

//! @file SystemSolvePardiso.ipp

#include<algorithm>
#include<vector>

// PARDISO prototype.
#define _WIN64

#if defined(_WIN32) || defined(_WIN64)
#define pardiso_ PARDISO
#else
#define PARDISO pardiso_
#endif

//------------------------------------------------------------------------------
void corlib::SystemSolvePardiso::pardisosv( const SolveMode_      solveMode,
                                            std::vector<int>    & colInd,
                                            std::vector<int>    & rowPtr,
                                            std::vector<double> & Anz)
{

    MKL_INT  n    = numDofs_;
    MKL_INT  nnz  = rowPtr[n]; // number of non-zeros

    MKL_INT mtype = 11;        // Real unsymmetric matrix
    MKL_INT nrhs  = 1;        // Number of right hand sides 

    MKL_INT maxfct = 1;      // Maximum number of numerical factorizations
    MKL_INT mnum   = 1;      // Which factorization to use
    MKL_INT msglvl = 0;      // Print statistical information in file 
    MKL_INT error  = 0;      // Initialize error flag 
    
    if ( solveMode == COMPLETE or solveMode == FACTORISE ) {

        /* -------------------------------------------------------------------- */
        /* .. Setup Pardiso control parameters. */
        /* -------------------------------------------------------------------- */
        std::fill(iparm_, iparm_+64, 0); // initialize with zeros

        iparm_[0] = 1;     // No solver default 
        iparm_[1] = 2;     // Fill-in reordering from METIS 
        iparm_[2] = mkl_get_max_threads(); // Numbers of processors, value of OMP_NUM_THREADS 
        iparm_[3] = 0;     // No iterative-direct algorithm
        iparm_[4] = 0;     // No user fill-in reducing permutation 
        iparm_[5] = 0;     // Write solution into x 
        iparm_[6] = 0;     // Not in use 
        iparm_[7] = 2;     // Max numbers of iterative refinement steps 
        iparm_[8] = 0;     // Not in use 
        iparm_[9] = 13;    // Perturb the pivot elements with 1E-13 
        iparm_[10] = 1;    // Use nonsymmetric permutation and scaling MPS 
        iparm_[11] = 0;    // Do not solve with transposed matrix
        iparm_[12] = 0;    // Maximum weighted matching algorithm is switched-off (default for symmetric).
        //   Try iparm_[12] = 1 in case of inappropriate accuracy 
        iparm_[13] = 0;    // Output: Number of perturbed pivots 
        iparm_[14] = 0;    // Not in use 
        iparm_[15] = 0;    // Not in use 
        iparm_[16] = 0;    // Not in use 
        iparm_[17] = -1;   // Output: Number of nonzeros in the factor LU 
        iparm_[18] = -1;   // Output: Mflops for LU factorization 
        iparm_[19] = 0;    // Output: Numbers of CG Iterations 

        // make row and column indices 1-based
        std::transform(rowPtr.begin(), rowPtr.end(), rowPtr.begin(), 
                       std::bind2nd(std::plus<double>(), 1));
        std::transform(colInd.begin(), colInd.end(), colInd.begin(), 
                       std::bind2nd(std::plus<double>(), 1));


        /* -------------------------------------------------------------------- */
        /* .. Initialize the internal solver memory pointer. This is only */
        /* necessary for the FIRST call of the PARDISO solver. */
        /* -------------------------------------------------------------------- */

        /* Internal solver memory pointer pt, 
           32-bit: int pt_[64]; 64-bit: long int pt_[64]
           or void *pt_[64] should be OK on both architectures */
        for ( unsigned i = 0; i < 64; i ++ ) pt_[i] = 0;

        /* -------------------------------------------------------------------- */
        /* .. Reordering and Symbolic Factorization. This step also allocates */
        /* all memory that is necessary for the factorization. */
        /* -------------------------------------------------------------------- */
        MKL_INT phase = 11;
    
        double ddum; // Double dummy 
        MKL_INT idum; // Integer dummy. 
    
        PARDISO (pt_, &maxfct, &mnum, &mtype, &phase,
                 &n, &(Anz[0]), &(rowPtr[0]), &(colInd[0]), &idum, &nrhs,
                 iparm_, &msglvl, &ddum, &ddum, &error);

        if (error != 0) {
            std::cout << std::endl
                      << "ERROR during analysis phase "  << error 
                      << std::endl;
            exit(1);
        }

        /* -------------------------------------------------------------------- */
        /* .. Numerical factorization. */
        /* -------------------------------------------------------------------- */
        phase = 22;
        PARDISO (pt_, &maxfct, &mnum, &mtype, &phase,
                 &n, &(Anz[0]), &(rowPtr[0]), &(colInd[0]), &idum, &nrhs,
                 iparm_, &msglvl, &ddum, &ddum, &error);

        if (error != 0) {
            std::cout << std::endl
                      << "ERROR during factorization phase "  << error 
                      << std::endl;
            exit(2);
        }

    }

    if ( solveMode == COMPLETE or solveMode == BACKSUBSTITUTE ) {
            
        /* -------------------------------------------------------------------- */
        /* .. Back substitution and iterative refinement. */
        /* -------------------------------------------------------------------- */
        MKL_INT phase = 33;
        iparm_[7] = 1; /* Max numbers of iterative refinement steps. */

        std::vector<double> rhs( F_.begin(), F_.end() );

        MKL_INT idum; // Integer dummy. 
   
        PARDISO (pt_, &maxfct, &mnum, &mtype, &phase,
                 &n, &(Anz[0]), &(rowPtr[0]), &(colInd[0]), &idum, &nrhs,
                 iparm_, &msglvl, &(rhs[0]), &(F_[0]), &error);

        if (error != 0) {
            std::cout << std::endl
                      << "ERROR during solution phase "  << error 
                      << std::endl;
            exit(3);
        }

    }

    if ( solveMode == COMPLETE or solveMode == RELEASE ) {

        /* -------------------------------------------------------------------- */
        /* .. Termination and release of memory. */
        /* -------------------------------------------------------------------- */
        MKL_INT phase = -1; /* Release internal memory. */

        double ddum; // Double dummy 
        MKL_INT idum; // Integer dummy. 

        PARDISO (pt_, &maxfct, &mnum, &mtype, &phase,
                 &n, &ddum, &(rowPtr[0]), &(colInd[0]), &idum, &nrhs,
                 iparm_, &msglvl, &ddum, &ddum, &error);

    }

    return;
}

