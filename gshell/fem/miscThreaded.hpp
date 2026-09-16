//
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//
//                 Computational Structural Mechanics Lab
//                         University of Cambridge
//
// This software is copyrighted and all rights are retained by the CSMLab. All
// use, disclosure, and/or reproduction of any part not expressly authorized by
// F Cirak is prohibited. (C) 2010.
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//

//! @author Burkhard Bornemann, Kosala Bandara
//! @date   2010

#ifndef gshell_fem_miscthreaded_h
#define gshell_fem_miscthreaded_h

//------------------------------------------------------------------------------
#include <vector>
#include <string>
#include <iostream>
#include <algorithm>

#include <omp.h>

//------------------------------------------------------------------------------
// declarations
namespace gshell {
    namespace fem {

        template< typename INPITER, typename FUNC >
        void for_each_threaded( INPITER begin, INPITER end, FUNC f );

    }
}

//==============================================================================
//! For-each loop OpenMP-threaded
//!
//! \param[in]   begin  Input iterator to the initial position of input sequence
//! \param[in]   end    Input iterator to the final position of input sequence
//! \param[in]   f      Unary function taking an element in the range as argument
template< typename INPITER, typename FUNC >
void gshell::fem::for_each_threaded( INPITER begin, INPITER end, FUNC f )
{
#pragma omp parallel
    {
        const int numThreads = omp_get_num_threads( );
        const int threadId = omp_get_thread_num( );
        const unsigned numElements = std::distance( begin, end );
        const unsigned elemsPerThread = numElements / numThreads;
        INPITER beginThread = begin;
        std::advance( beginThread, threadId*elemsPerThread );
        INPITER endThread = begin;
        std::advance( endThread, ( ( threadId == (numThreads-1) ) ? numElements : (threadId+1)*elemsPerThread ) );
        std::for_each( beginThread, endThread, f );
    }
    return;
}

#endif
