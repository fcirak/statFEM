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

#ifndef gshell_fem_linalg_h
#define gshell_fem_linalg_h

#include <cmath>

#include <corlib/verify.hpp>

namespace gshell {
    namespace fem {

        //! Linear algebra routines for small vectors and matrices
        namespace la {

            //==================================================================
            // Vector methods

            //------------------------------------------------------------------
            //! Set vector to zero
            //!
            //! \tparam          N   Number of vector components
            //! \param[in,out]   a   Vector to blanked
            template< unsigned N > inline
            void zero( double a[ N ] )
            {
                for ( unsigned i = 0; i < N; ++i )
                    a[ i ] = 0.0;
                return;
            }

            //------------------------------------------------------------------
            //! Multiply with scalar
            //!
            //! \tparam          N      Number of vector components
            //! \param[in,out]   a      Vector to be multiplied
            //! \param[in]       fact   Scalar to multiply with
            template< unsigned N > inline
            void multiply( double a[ N ],
                           const double & fact )
            {
                for ( unsigned i = 0; i < N; ++i )
                    a[ i ] *= fact;
                return;
            }

            //------------------------------------------------------------------
            //! Divide with scalar
            //!
            //! \tparam          N      Number of vector components
            //! \param[in,out]   a      Vector to be divided
            //! \param[in]       denom  Scalar to divide with
            template< unsigned N > inline
            void divide( double a[ N ],
                         const double & denom )
            {
                for ( unsigned i = 0; i < N; ++i )
                    a[ i ] /= denom;
                return;
            }

            //------------------------------------------------------------------
            //! Euclidian norm
            //!
            //! \tparam          N      Number of vector components
            //! \param[out]      s      Vector norm
            //! \param[in]       a      Vector whose norm is computed
            template< unsigned N > inline
            void norm2( double & s, const double a[ N ] )
            {
                s = 0.0;
                for ( unsigned i = 0; i < N; ++i )
                    s += a[ i ] * a[ i ];
                s = std::sqrt( s );
                return;
            }

            //------------------------------------------------------------------
            //! Inner/dot product
            //!
            //! \tparam          N      Number of vector components
            //! \param[out]      prod   Inner product
            //! \param[in]       a      First vector
            //! \param[in]       b      Second vector
            template< unsigned N > inline
            void innerProduct( double & prod,
                               const double a[ N ],
                               const double b[ N ] )
            {
                prod = 0.0;
                for ( unsigned i = 0; i < N; ++i )
                    prod += a[ i ] * b[ i ];
                return;
            }

            //------------------------------------------------------------------
            //! Add an inner/dot product
            //!
            //! \tparam          N      Number of vector components
            //! \param[in,out]   prod   Inner product is added onto this
            //! \param[in]       a      First vector
            //! \param[in]       b      Second vector
            template< unsigned N > inline
            void addInnerProduct( double & prod,
                                  const double a[ N ],
                                  const double b[ N ] )
            {
                for ( unsigned i = 0; i < N; ++i )
                    prod += a[ i ] * b[ i ];
                return;
            }

            //------------------------------------------------------------------
            //! Subtract an inner/dot product
            //!
            //! \tparam          N      Number of vector components
            //! \param[in,out]   prod   Inner product is subtracted from this
            //! \param[in]       a      First vector
            //! \param[in]       b      Second vector
            template< unsigned N > inline
            void subtractInnerProduct( double & prod,
                                       const double a[ N ],
                                       const double b[ N ] )
            {
                for ( unsigned i = 0; i < N; ++i )
                    prod -= a[ i ] * b[ i ];
                return;
            }

            //------------------------------------------------------------------
            //! Assign/copy vector from source vector
            //!
            //! \tparam          N      Number of vector components
            //! \param[out]      a      Assigned vector
            //! \param[in]       b      Source vector
            template< unsigned N > inline
            void assign( double a[ N ],
                         const double b[ N ] )
            {
                for ( unsigned i = 0; i < N; ++i )
                    a[ i ] = b[ i ];
                return;
            }

            //------------------------------------------------------------------
            //! Add a vector onto another
            //!
            //! \tparam          N      Number of vector components
            //! \param[in,out]   a      Vector onto which is added
            //! \param[in]       b      Source vector
            template< unsigned N > inline
            void add( double a[ N ],
                      const double b[ N ] )
            {
                for ( unsigned i = 0; i < N; ++i )
                    a[ i ] += b[ i ];
                return;
            }

            //------------------------------------------------------------------
            //! Update a vector by adding a scaled source vector
            //!
            //!\f[
            //!    \vec{a} := \vec{a} + fact * \vec{b}
            //!\f]
            //!
            //! \tparam          N      Number of vector components
            //! \param[in,out]   a      Vector onto which is added
            //! \param[in]       fact   Scalar
            //! \param[in]       b      Source vector
            template< unsigned N > inline
            void update( double a[ N ],
                         const double & fact,
                         const double b[ N ] )
            {
                for ( unsigned i = 0; i < N; ++i )
                    a[ i ] += fact * b[ i ];
                return;
            }

            //------------------------------------------------------------------
            //! Cross-product of two vectors
            //!
            //!\f[
            //!    \vec{c} := \vec{a} \times \vec{b}
            //!\f]
            //!
            //! \param[out]     c      Cross-product vector
            //! \param[in]      a      First source vector
            //! \param[in]      b      Second source vector
            inline
            void crossProduct( double c[ 3 ],
                               const double a[ 3 ],
                               const double b[ 3 ] )
            {
                c[ 0 ] = a[ 1 ] * b[ 2 ] - a[ 2 ] * b[ 1 ];
                c[ 1 ] = a[ 2 ] * b[ 0 ] - a[ 0 ] * b[ 2 ];
                c[ 2 ] = a[ 0 ] * b[ 1 ] - a[ 1 ] * b[ 0 ];
                return;
            }

            //------------------------------------------------------------------
            //! Add cross-product of two vectors
            //!
            //!\f[
            //!    \vec{c} := \vec{c} + \vec{a} \times \vec{b}
            //!\f]
            //!
            //! \param[in,out]  c      Cross-product vector is added onto
            //! \param[in]      a      First source vector
            //! \param[in]      b      Second source vector
            inline
            void addCrossProduct( double c[ 3 ],
                                  const double a[ 3 ],
                                  const double b[ 3 ] )
            {
                c[ 0 ] += a[ 1 ] * b[ 2 ] - a[ 2 ] * b[ 1 ];
                c[ 1 ] += a[ 2 ] * b[ 0 ] - a[ 0 ] * b[ 2 ];
                c[ 2 ] += a[ 0 ] * b[ 1 ] - a[ 1 ] * b[ 0 ];
                return;
            }

            //------------------------------------------------------------------
            //! Write vector to stream
            //!
            //! \tparam          N      Number of vector components
            //! \param[in,out]   os     Output stream
            //! \param[in]       a      Vector to write
            template< unsigned N >
            void write( std::ostream & os,
                        const double a[ N ] )
            {
                std::cout << "[" << N << "](";
                for ( unsigned i = 0; i < N; ++i ) {
                    os << a[ i ];
                    if ( i < N-1 ) os << ",";
                }
                os << ")";
                return;
            }


            //==================================================================
            // Matrix methods

            //------------------------------------------------------------------
            //! Determinant of square matrix
            //!
            //! \tparam          N      Number of rows/column
            //! \param[in]       m      Square matrix
            //! \return                 Determinant
            template< unsigned N >
            double determinant( const double m[ N ][ N ] );

            template< >
            double determinant( const double m[ 2 ][ 2 ] )
            {
                return 
                    m[ 0 ][ 0 ] * m[ 1 ][ 1 ] - m[ 0 ][ 1 ] * m[ 1 ][ 0 ];
            }

            template< >
            double determinant( const double m[ 3 ][ 3 ] )
            {
                return 
                    m[ 0 ][ 0 ] * m[ 1 ][ 1 ] * m[ 2 ][ 2 ] +
                    m[ 0 ][ 1 ] * m[ 1 ][ 2 ] * m[ 2 ][ 0 ] +
                    m[ 0 ][ 2 ] * m[ 1 ][ 0 ] * m[ 2 ][ 1 ] -
                    m[ 0 ][ 0 ] * m[ 1 ][ 2 ] * m[ 2 ][ 1 ] -
                    m[ 0 ][ 1 ] * m[ 1 ][ 0 ] * m[ 2 ][ 2 ] -
                    m[ 0 ][ 2 ] * m[ 1 ][ 1 ] * m[ 2 ][ 0 ];
            }

            //------------------------------------------------------------------
            //! Zero matrix
            //!
            //! \tparam          NROW   Number of matrix rows
            //! \tparam          NCOL   Number of matrix columns
            //! \param[out]      a      Matrix to blank
            template< unsigned NROW, unsigned NCOL > inline
            void zero( double a[ NROW ][ NCOL ] )
            {
                for ( unsigned i = 0; i < NROW; ++i )
                    for ( unsigned j = 0; j < NCOL; ++j )
                        a[ i ][ j ] = 0.0;
                return;
            }

            //------------------------------------------------------------------
            //! Identity matrix
            //!
            //! \tparam          N      Number of matrix rows/columns
            //! \param[out]      a      Idenity matrix
            template< unsigned N > inline
            void identity( double a[ N ][ N ] )
            {
                for ( unsigned i = 0; i < N; ++i )
                    for ( unsigned j = 0; j < N; ++j )
                        a[ i ][ j ] = 0.0;
                for ( unsigned i = 0; i < N; ++i )
                    a[ i ][ i ] = 1.;
                return;
            }

            //------------------------------------------------------------------
            //! Multiply matrix by scalar
            //!
            //! \tparam          NROW   Number of matrix rows
            //! \tparam          NCOL   Number of matrix columns
            //! \param[in,out]   a      Matrix
            //! \param[in]       fact   Scalar
            template< unsigned NROW, unsigned NCOL > inline
            void multiply( double a[ NROW ][ NCOL ],
                           const double & fact )
            {
                for ( unsigned i = 0; i < NROW; ++i )
                    for ( unsigned j = 0; j < NCOL; ++j )
                        a[ i ][ j ] *= fact;
                return;
            }

            //------------------------------------------------------------------
            //! Divide matrix by scalar
            //!
            //! \tparam          NROW   Number of matrix rows
            //! \tparam          NCOL   Number of matrix columns
            //! \param[in,out]   a      Matrix
            //! \param[in]       denom  Denominator
            template< unsigned NROW, unsigned NCOL > inline
            void divide( double a[ NROW ][ NCOL ],
                         const double & denom )
            {
                for ( unsigned i = 0; i < NROW; ++i )
                    for ( unsigned j = 0; j < NCOL; ++j )
                        a[ i ][ j ] /= denom;
                return;
            }

            //------------------------------------------------------------------
            //! Assign/copy matrix
            //!
            //! \tparam          NROW   Number of matrix rows
            //! \tparam          NCOL   Number of matrix columns
            //! \param[in,out]   a      Matrix to be assigned
            //! \param[in]       b      Source matrix
            template< unsigned NROW, unsigned NCOL > inline
            void assign( double a[ NROW ][ NCOL ],
                         const double b[ NROW ][ NCOL ] )
            {
                for ( unsigned i = 0; i < NROW; ++i )
                    for ( unsigned j = 0; j < NCOL; ++j )
                        a[ i ][ j ] = b[ i ][ j ];
                return;
            }

            //------------------------------------------------------------------
            //! Assign/copy matrix
            //!
            //! \tparam          NROW   Number of matrix rows
            //! \tparam          NCOL   Number of matrix columns
            //! \param[in,out]   a      Matrix to be assigned
            //! \param[in]       b      Source matrix
            template< unsigned NROW, unsigned NCOL > inline
            void assign( double * a,
                         const double b[ NROW ][ NCOL ] )
            {
                for ( unsigned i = 0; i < NROW; ++i )
                    for ( unsigned j = 0; j < NCOL; ++j )
                        a[ i*NROW+j ] = b[ i ][ j ];
                return;
            }

            //------------------------------------------------------------------
            //! Assign matrix by top-left submatrix block
            //!
            //! \tparam          NROW   Number of target matrix rows
            //! \tparam          NCOL   Number of target matrix columns
            //! \tparam          NROW1  Number of source matrix rows (NROW1>=NROW)
            //! \tparam          NCOL1  Number of source matrix columns (NCOL1>=NCOL)
            //! \param[in,out]   a      Matrix to be assigned
            //! \param[in]       b      Source matrix
            template< unsigned NROW, unsigned NCOL, unsigned NROW1, unsigned NCOL1 > inline
            void subAssign( double a[ NROW ][ NCOL ],
                            const double b[ NROW1 ][ NCOL1 ] )
            {
                assert( NROW1 >= NROW );
                assert( NCOL1 >= NCOL );
                for ( unsigned i = 0; i < NROW; ++i )
                    for ( unsigned j = 0; j < NCOL; ++j )
                        a[ i ][ j ] = b[ i ][ j ];
                return;
            }

            //------------------------------------------------------------------
            //! Add matrix to another
            //!
            //! \tparam          NROW   Number of matrix rows
            //! \tparam          NCOL   Number of matrix columns
            //! \param[in,out]   a      Matrix onto with is added
            //! \param[in]       b      Source matrix
            template< unsigned NROW, unsigned NCOL > inline
            void add( double a[ NROW ][ NCOL ],
                      const double b[ NROW ][ NCOL ] )
            {
                for ( unsigned i = 0; i < NROW; ++i )
                    for ( unsigned j = 0; j < NCOL; ++j )
                        a[ i ][ j ] += b[ i ][ j ];
                return;
            }

            //------------------------------------------------------------------
            //! Update matrix by adding scaled source matrix
            //!
            //!\f[
            //!    A := A + fact * B
            //!\f]
            //!
            //! \tparam          NROW   Number of matrix rows
            //! \tparam          NCOL   Number of matrix columns
            //! \param[in,out]   a      Matrix onto with is added
            //! \param[in]       fact   Scalar factor
            //! \param[in]       b      Source matrix
            template< unsigned NROW, unsigned NCOL > inline
            void update( double a[ NROW ][ NCOL ],
                         const double & fact,
                         const double b[ NROW ][ NCOL ] )
            {
                for ( unsigned i = 0; i < NROW; ++i )
                    for ( unsigned j = 0; j < NCOL; ++j )
                        a[ i ][ j ] += fact * b[ i ][ j ];
                return;
            }

            //------------------------------------------------------------------
            //! Invert square matrix
            //!
            //! \tparam          N      Number of rows/column
            //! \param[in, out]  m      Square matrix, inverse on exit
            //! \return                 Determinant
            template< unsigned N >
            double inverse( double m[ N ][ N ] );

            template< >
            double inverse( double m[ 2 ][ 2 ] )
            {
                double detM = gshell::fem::la::determinant( m );
                FTL_VERIFY( std::fabs( detM ) > 1.e-15 );

                const double aux = m[ 0 ][ 0 ];
                m[ 0 ][ 0 ] = m[ 1 ][ 1 ] / detM;
                m[ 0 ][ 1 ] /= -detM;
                m[ 1 ][ 0 ] /= -detM;
                m[ 1 ][ 1 ] = aux / detM;

                return detM;
            }

            template< >
            double inverse( double m[ 3 ][ 3 ] )
            {
                const double detM = gshell::fem::la::determinant( m );
                FTL_VERIFY_DESCRIPTIVE( std::fabs( detM ) > 1.e-15, "detM = %g\n", detM );

                // use a copy here
                double t[ 3 ][ 3 ];
                gshell::fem::la::assign<3,3>( t, m );
                // main-diagonal
                m[ 0 ][ 0 ] = ( t[ 1 ][ 1 ] * t[ 2 ][ 2 ] - t[ 1 ][ 2 ] * t[ 2 ][ 1 ] ) / detM;
                m[ 1 ][ 1 ] = ( t[ 0 ][ 0 ] * t[ 2 ][ 2 ] - t[ 0 ][ 2 ] * t[ 2 ][ 0 ] ) / detM;
                m[ 2 ][ 2 ] = ( t[ 0 ][ 0 ] * t[ 1 ][ 1 ] - t[ 0 ][ 1 ] * t[ 1 ][ 0 ] ) / detM;
                // off-diagonal
                m[ 0 ][ 1 ] = ( t[ 0 ][ 2 ] * t[ 2 ][ 1 ] - t[ 0 ][ 1 ] * t[ 2 ][ 2 ] ) / detM;
                m[ 0 ][ 2 ] = ( t[ 0 ][ 1 ] * t[ 1 ][ 2 ] - t[ 0 ][ 2 ] * t[ 1 ][ 1 ] ) / detM;
                m[ 1 ][ 0 ] = ( t[ 1 ][ 2 ] * t[ 2 ][ 0 ] - t[ 1 ][ 0 ] * t[ 2 ][ 2 ] ) / detM;
                m[ 1 ][ 2 ] = ( t[ 0 ][ 2 ] * t[ 1 ][ 0 ] - t[ 0 ][ 0 ] * t[ 1 ][ 2 ] ) / detM;
                m[ 2 ][ 0 ] = ( t[ 1 ][ 0 ] * t[ 2 ][ 1 ] - t[ 1 ][ 1 ] * t[ 2 ][ 0 ] ) / detM;
                m[ 2 ][ 1 ] = ( t[ 0 ][ 1 ] * t[ 2 ][ 0 ] - t[ 0 ][ 0 ] * t[ 2 ][ 1 ] ) / detM;
        
                return detM;
            }

            //------------------------------------------------------------------
            //! Write vector to stream -- format follows uBLAS IO
            //!
            //! \tparam          NROW   Number of matrix rows
            //! \tparam          NCOL   Number of matrix columns
            //! \param[in,out]   os     Output stream
            //! \param[in]       b      Matrix to write
            template< unsigned NROW, unsigned NCOL >
            void write( std::ostream & os,
                        const double b[ NROW ][ NCOL ] )
            {
                os << "[" << NROW << "," << NCOL << "](";
                for ( unsigned i = 0; i < NROW; ++i ) {
                    os << "(";
                    for ( unsigned j = 0; j < NCOL; ++j ) {
                        os << b[ i ][ j ];
                        if ( j < NCOL-1  ) os << ",";
                    }
                    os << ")";
                    if ( i < NROW-1 ) os << ",";
                }
                os << ")";
                return;
            }

        
            //==================================================================
            // Matrix-Vector methods

            //------------------------------------------------------------------
            //! Matrix-vector product
            //!
            //!\f[
            //!    a_i := B_{ij} \, c_j
            //!\f]
            //!
            //! \tparam         N      Number of vector components
            //! \param[out]     a      Resulting vector
            //! \param[in]      nInn   Inner number of indices
            //! \param[in]      b      Input (N)x(nInn)-matrix (row-major flat storage)
            //! \param[in]      c      Input (nInn)-vector
            template< unsigned N > inline
            void productN( double a[ N ],
                           const unsigned nInn,
                           const double * b,
                           const double * c )
            {
                for ( unsigned i = 0; i < N; ++i ) {
                    a[ i ] = 0.0;
                    for ( unsigned k = 0; k < nInn; ++k ) {
                        a[ i ] += b[ i*nInn+k ] * c[ k ];
                    }
                }
                return;
            }
        
            //------------------------------------------------------------------
            //! Matrix-vector product
            //!
            //!\f[
            //!    a_i := B_{ij} \, c_j
            //!\f]
            //!
            //! \tparam         N      Number of vector components
            //! \tparam         NINN   Number of matrix columns
            //! \param[out]     a      Resulting vector
            //! \param[in]      b      Input (N)x(NINN)-matrix
            //! \param[in]      c      Input (NINN)-vector
            template< unsigned N, unsigned NINN > inline
            void productN( double a[ N ],
                           const double b[ N ][ NINN ],
                           const double c[ NINN ] )
            {
                for ( unsigned i = 0; i < N; ++i ) {
                    a[ i ] = 0.0;
                    for ( unsigned k = 0; k < NINN; ++k ) {
                        a[ i ] += b[ i ][ k ] * c[ k ];
                    }
                }
                return;
            }

            //------------------------------------------------------------------
            //! Transposed matrix-vector product
            //!
            //!\f[
            //!    a_i := (B_{ji})^T \, c_j
            //!\f]
            //!
            //! \tparam         N      Number of vector components
            //! \tparam         NINN   Number of matrix rows
            //! \param[out]     a      Resulting vector
            //! \param[in]      b      Input matrix
            //! \param[in]      c      Input vector
            template< unsigned N, unsigned NINN > inline
            void productT( double a[ N ],
                           const double b[ NINN ][ N ],
                           const double c[ NINN ] )
            {
                for ( unsigned i = 0; i < N; ++i ) {
                    a[ i ] = 0.0;
                    for ( unsigned k = 0; k < NINN; ++k ) {
                        a[ i ] += b[ k ][ i ] * c[ k ];
                    }
                }
                return;
            }

            //==================================================================
            // Matrix-Matrix methods

            //------------------------------------------------------------------
            //! Matrix-matrix product
            //!
            //!\f[
            //!    A_{ij} := B_{ik} \, C_{kj}
            //!\f]
            //!
            //! \tparam          NROW   Number of resulting matrix rows
            //! \tparam          NCOL   Number of resulting matrix columns
            //! \param[out]      a      Resulting matrix
            //! \param[in]       nInn   Inner number of indices
            //! \param[in]       b      Left input (NROW)x(nInn)-matrix  (row-major flat storage)
            //! \param[in]       c      Right input (nInn)x(NCOL)-matrix  (row-major flat storage)
            template< unsigned NROW, unsigned NCOL > inline
            void productNN( double a[ NROW ][ NCOL ],
                            const unsigned nInn,
                            const double * b,
                            const double * c )
            {
                for ( unsigned i = 0; i < NROW; ++i ) {
                    for ( unsigned j = 0; j < NCOL; ++j ) {
                        double s = 0.0;
                        for ( unsigned k = 0; k < nInn; ++k ) {
                            s += b[ i*nInn+k ] * c[ k*NCOL+j ];
                        }
                        a[ i ][ j ] = s;
                    }
                }
                return;
            }

            //------------------------------------------------------------------
            //! Matrix-matrix product
            //!
            //!\f[
            //!    A_{ij} := B_{ik} \, C_{kj}
            //!\f]
            //!
            //! \tparam          NROW   Number of resulting matrix rows
            //! \tparam          NCOL   Number of resulting matrix columns
            //! \tparam          NINN   Number of inner columns/rows
            //! \param[out]      a      Resulting matrix
            //! \param[in]       b      Left input matrix
            //! \param[in]       c      Right input matrix
            template< unsigned NROW, unsigned NCOL, unsigned NINN > inline
            void productNN( double a[ NROW ][ NCOL ],
                            const double b[ NROW ][ NINN ],
                            const double c[ NINN ][ NCOL] )
            {
                for ( unsigned i = 0; i < NROW; ++i ) {
                    for ( unsigned j = 0; j < NCOL; ++j ) {
                        double s = 0.0;
                        for ( unsigned k = 0; k < NINN; ++k ) {
                            s += b[ i ][ k ] * c[ k ][ j ];
                        }
                        a[ i ][ j ] = s;
                    }
                }
                return;
            }

            //------------------------------------------------------------------
            //! Matrix-tranposed matrix product
            //!
            //!\f[
            //!    A_{ij} := B_{ik} \, (C_{jk})^T
            //!\f]
            //!
            //! \tparam          NROW   Number of resulting matrix rows
            //! \tparam          NCOL   Number of resulting matrix columns
            //! \param[out]      a      Resulting matrix
            //! \param[in]       b      Left input matrix
            //! \param[in]       c      Right input matrix
            template< unsigned NROW, unsigned NCOL > inline
            void productNT( double a[ NROW ][ NCOL ],
                            const unsigned nInn,
                            const double * b,
                            const double * c )
            {
                for ( unsigned i = 0; i < NROW; ++i ) {
                    for ( unsigned j = 0; j < NCOL; ++j ) {
                        double s = 0.0;
                        for ( unsigned k = 0; k < nInn; ++k ) {
                            s += b[ i*nInn+k ] * c[ j*nInn+k ];
                        }
                        a[ i ][ j ] = s;
                    }
                }
                return;
            }

            //------------------------------------------------------------------
            //! Matrix - transpose matrix product
            //!
            //!\f[
            //!    A_{ij} := B_{ik} \, (C_{jk})^T
            //!\f]
            //!
            //! \tparam          NROW   Number of resulting matrix rows
            //! \tparam          NCOL   Number of resulting matrix columns
            //! \tparam          NINN   Number of inner columns/rows
            //! \param[out]      a      Resulting matrix
            //! \param[in]       b      Left input matrix
            //! \param[in]       c      Right input matrix
            template< unsigned NROW, unsigned NCOL, unsigned NINN > inline
            void productNT( double a[ NROW ][ NCOL ],
                            const double b[ NROW ][ NINN ],
                            const double c[ NCOL ][ NINN ] )
            {
                for ( unsigned i = 0; i < NROW; ++i ) {
                    for ( unsigned j = 0; j < NCOL; ++j ) {
                        double s = 0.0;
                        for ( unsigned k = 0; k < NINN; ++k ) {
                            s += b[ i ][ k ] * c[ j ][ k ];
                        }
                        a[ i ][ j ] = s;
                    }
                }
                return;
            }

            //------------------------------------------------------------------
            //! Transpose Matrix - matrix product
            //!
            //!\f[
            //!    A_{ij} := (B_{ki})^T \, C_{kj}
            //!\f]
            //!
            //! \tparam          NROW   Number of resulting matrix rows
            //! \tparam          NCOL   Number of resulting matrix columns
            //! \tparam          NINN   Number of inner columns/rows
            //! \param[out]      a      Resulting matrix
            //! \param[in]       b      Left input matrix
            //! \param[in]       c      Right input matrix
            template< unsigned NROW, unsigned NCOL, unsigned NINN > inline
            void productTN( double a[ NROW ][ NCOL ],
                            const double b[ NINN ][ NROW ],
                            const double c[ NINN ][ NCOL ] )
            {
                for ( unsigned i = 0; i < NROW; ++i ) {
                    for ( unsigned j = 0; j < NCOL; ++j ) {
                        double s = 0.0;
                        for ( unsigned k = 0; k < NINN; ++k ) {
                            s += b[ k ][ i ] * c[ k ][ j ];
                        }
                        a[ i ][ j ] = s;
                    }
                }
                return;
            }

            //------------------------------------------------------------------
            //! Transpose matrix - transpose matrix - product
            //!
            //!\f[
            //!    A_{ij} := (B_{ki})^T \, (C_{jk})^T
            //!\f]
            //!
            //! \tparam          NROW   Number of resulting matrix rows
            //! \tparam          NCOL   Number of resulting matrix columns
            //! \tparam          NINN   Number of inner columns/rows
            //! \param[out]      a      Resulting matrix
            //! \param[in]       b      Left input matrix
            //! \param[in]       c      Right input matrix
            template< unsigned NROW, unsigned NCOL, unsigned NINN > inline
            void productTT( double a[ NROW ][ NCOL ],
                            const double b[ NINN ][ NROW ],
                            const double c[ NCOL ][ NINN ] )
            {
                for ( unsigned i = 0; i < NROW; ++i ) {
                    for ( unsigned j = 0; j < NCOL; ++j ) {
                        double s = 0.0;
                        for ( unsigned k = 0; k < NINN; ++k ) {
                            s += b[ k ][ i ] * c[ j ][ k ];
                        }
                        a[ i ][ j ] = s;
                    }
                }
                return;
            }

            //------------------------------------------------------------------
            //! Inner matrix product / double contraction
            //!
            //! \tparam          NROW   Number of matrix rows
            //! \tparam          NCOL   Number of matrix columns
            //! \param[out]      prod   Resulting inner product
            //! \param[in]       a      First input matrix
            //! \param[in]       b      Second input matrix
            template< unsigned NROW, unsigned NCOL > inline
            void innerProduct( double & prod,
                               const double a[ NROW ][ NCOL ],
                               const double b[ NROW ][ NCOL ] )
            {
                prod = 0.0;
                for ( unsigned i = 0; i < NROW; ++i )
                    for ( unsigned j = 0; j < NCOL; ++j )
                        prod += a[ i ][ j ] * b[ i ][ j ];
                return;
            }

            //------------------------------------------------------------------
            //! Add an inner matrix product / double contraction
            //!
            //! \tparam          NROW   Number of matrix rows
            //! \tparam          NCOL   Number of matrix columns
            //! \param[in,out]   prod   Add resulting inner product
            //! \param[in]       a      First input matrix
            //! \param[in]       b      Second input matrix
            template< unsigned NROW, unsigned NCOL > inline
            void addInnerProduct( double & prod,
                                  const double a[ NROW ][ NCOL ],
                                  const double b[ NROW ][ NCOL ] )
            {
                for ( unsigned i = 0; i < NROW; ++i )
                    for ( unsigned j = 0; j < NCOL; ++j )
                        prod += a[ i ][ j ] * b[ i ][ j ];
                return;
            }



            //==================================================================
            // Tensor methods

            //------------------------------------------------------------------
            //! Zero 4-th order tensor
            //!
            //! \tparam          N0     First index dimension
            //! \tparam          N1     Second index dimension
            //! \tparam          N2     Third index dimension
            //! \tparam          N3     Fourth index dimension
            //! \param[out]      a      Tensor to be blanked
            template< unsigned N0, unsigned N1, unsigned N2, unsigned N3 >
            void zero( double a[ N0 ][ N1 ][ N2 ][ N3 ] )
            {
                for ( unsigned i = 0; i < N0; ++i )
                    for ( unsigned j = 0; j < N1; ++j )
                        for ( unsigned k = 0; k < N2; ++k )
                            for ( unsigned l = 0; l < N3; ++l )
                                a[ i ][ j ][ k ][ l ] = 0.0;
                return;
            }

        }
    }
}

#endif
