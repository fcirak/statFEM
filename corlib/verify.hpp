// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//
//                 Computational Structural Mechanics Lab
//                         University of Cambridge
//
// This software is copyrighted and all rights are retained by the University of
// Cambridge. All use, disclosure, and/or reproduction of any part not expressly
// authorised by the University of Cambridge is prohibited. (C) 2010
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

//! @file verify.hpp

//------------------------------------------------------------------------------
// copied and adapted from /usr/include/boost/assert.hpp
#ifndef corlib_verify_h
#define corlib_verify_h

#include <boost/current_function.hpp>
#include <boost/static_assert.hpp>
#include <iostream>
#include <cstdlib>
#include <cstdarg>
#include <cstdio>

//------------------------------------------------------------------------------
// declaration
namespace corlib
{
    void assertion_failed( char const * expr, char const * function, 
                           char const * file, long line );
    void assertion_failed_descriptive( char const * expr, char const * function, 
                                       char const * file, long line, ... );
} 

//------------------------------------------------------------------------------
// set macro
//
#define FTL_VERIFY(expr) ((expr)? ((void)0): ::corlib::assertion_failed(#expr, BOOST_CURRENT_FUNCTION, __FILE__, __LINE__ ))

// The final -1 is needed to check wether __VA_ARGS__ is existent.
#define FTL_VERIFY_DESCRIPTIVE(expr, ...) ((expr)? ((void)0): ::corlib::assertion_failed_descriptive(#expr, BOOST_CURRENT_FUNCTION, __FILE__, __LINE__, __VA_ARGS__, -1 ))


// macro for static assertion with message
#ifndef BOOST_NO_STATIC_ASSERT
#  define FTL_STATIC_ASSERT_MSG( B, Msg ) static_assert(B, Msg)
#else
#  define FTL_STATIC_ASSERT_MSG( B, Msg ) BOOST_STATIC_ASSERT( B )
#endif



//------------------------------------------------------------------------------
//! Output to std::cerr with error message
//!
//! \param[in]  expr        Boolean expression which is verified to be true.
//!                         Program is stopped if expression is false.
//! \param[in]  function    Function name in which verification is stated.
//! \param[in]  file        File in which verification has been called.
//! \param[in]  line        Line number on which verification is located.
void corlib::assertion_failed( char const * expr, char const * function, 
                               char const * file, long line )
{
    // print error message
    std::cerr << "(EE) Assertion of \"" << expr << "\" in function \"" << function << "\", "
              << std::endl
              << "(EE) line " << line  << " of file \"" << file << "\" failed! " << std::endl;
    abort( );  // dump core
    exit( -1 );
}

//------------------------------------------------------------------------------
//! Output to std::cerr with error message and a description
//!
//! \param[in]  expr        Boolean expression which is verified to be true.
//!                         Program is stopped if expression is false.
//! \param[in]  function    Function name in which verification is stated.
//! \param[in]  file        File in which verification has been called.
//! \param[in]  line        Line number on which verification is located.
//! \param[in]  ...         List of variable arguments which is finished
//!                         by -1. The variable arguments have to follow
//!                         the arguments of printf( const char *, ... ),
//!                         i.e. a format string followed by a variable
//!                         number of items.
void corlib::assertion_failed_descriptive( char const * expr, char const * function, 
                                           char const * file, long line, ... )
{
    // print user-specified message contained in '...' if existent
    va_list ap;  va_start( ap, line );
    const int i = va_arg( ap, int );
    va_end( ap );
    if ( i != -1 ) {
        va_list ap;  va_start( ap, line );
        const char * fmt = va_arg( ap, char * );
        fprintf( stderr, "(EE) " );
        vfprintf( stderr, fmt, ap );
        va_end( ap );
    }
    // print error message
    assertion_failed( expr, function, file, line );
}

#endif
