// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//
//                 Computational Structural Mechanics Lab
//                         University of Cambridge
//
// This software is copyrighted and all rights are retained by the University of
// Cambridge. All use, disclosure, and/or reproduction of any part not expressly
// authorised by the University of Cambridge is prohibited. (C) 2010
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

//! @file linalg.hpp

#ifndef corlib_linalg_h
#define corlib_linalg_h
//------------------------------------------------------------------------------
#include <ostream>
#include <iomanip>
#include <cmath>
#include <algorithm>

#include <boost/numeric/ublas/triangular.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/multi_array.hpp>

#include <Eigen/Core>
#include <Eigen/Geometry> // For cross product
#include <Eigen/SparseCholesky>

#include <corlib/verify.hpp>
#include <corlib/eigenX.hpp>

namespace corlib {
    
    template< unsigned DIM >
    double trace( const eigenX::MatrixSd< DIM > & m );
} 


namespace corlib {

    //------------------------------------------------------------------------------
    //! @name Determinants of small matrices (up to 3x3)
    //@{

    //! Determinant of a 1x1 matrix -> the matrix entry itself
    //!
    //! \param[in]  m   matrix 1x1
    //! \return         determinant, i.e. the matrix itself
    double determinant( const eigenX::MatrixSd< 1 > & m )
    {
        return m( 0, 0 );
    }

    //! Determinant of a 2x2 matrix
    //!
    //! \param[in]  m   matrix 2x2
    //! \return         determinant
    double determinant( const eigenX::MatrixSd< 2 > & m )
    {
        return m( 0, 0 ) * m( 1, 1 ) - m( 0, 1 ) * m( 1, 0 );
    }

    //! Determinant of a 3x3 matrix
    //!
    //! \param[in]  m   matrix 3x3
    //! \return         determinant
    double determinant( const eigenX::MatrixSd< 3 > & m )
    {
        return 
            m( 0, 0 ) * m( 1, 1 ) * m( 2, 2 ) +
            m( 0, 1 ) * m( 1, 2 ) * m( 2, 0 ) +
            m( 0, 2 ) * m( 1, 0 ) * m( 2, 1 ) -
            m( 0, 0 ) * m( 1, 2 ) * m( 2, 1 ) -
            m( 0, 1 ) * m( 1, 0 ) * m( 2, 2 ) -
            m( 0, 2 ) * m( 1, 1 ) * m( 2, 0 );
    }

    //! Determinant of a 2x2 matrix provided as two vectors
    //!
    //! \param[in]  m0  first vector
    //! \param[in]  m1  second vector
    //! \return         determinant
    double determinant( const eigenX::VectorSd< 2 > & m0,
                        const eigenX::VectorSd< 2 > & m1 )
    {
        return m0( 0 ) * m1( 1 ) - m0( 1 ) * m1( 0 );
    }

    //@}

    //------------------------------------------------------------------------------
    //! @name Inverse matrices, stored in-place
    //@{

    //! Inverse of a 1x1 matrix
    //!
    //! \param[in,out]  m   matrix 1x1, its inverse on exit
    //! \return             determinant
    double inverse( eigenX::MatrixSd< 1 > & m )
    {
        const double detM =  m( 0, 0 );
        const double tolerance = 1.e-15;
        FTL_VERIFY( std::fabs( detM ) > tolerance );
        m( 0, 0 ) = 1./ detM;
        return detM;
    }

    //! Inverse of a 2x2 matrix
    //!
    //! \param[in,out]  m   matrix 2x2, its inverse on exit
    //! \return             determinant
    double inverse( eigenX::MatrixSd< 2 > & m )
    {
        const double tolerance = 1.e-15;
        const double detM = determinant( m );
        FTL_VERIFY( std::fabs( detM ) > tolerance );
        const double aux = m( 0, 0 );
        m( 0, 0 ) = m( 1, 1 ) / detM;
        m( 0, 1 ) /= -detM;
        m( 1, 0 ) /= -detM;
        m( 1, 1 ) = aux / detM;
        
        return detM;
    }

    //! Inverse of a 3x3 matrix
    //!
    //! \param[in,out]  m   matrix 3x3, its inverse on exit
    //! \return             determinant
    double inverse( eigenX::MatrixSd< 3 > & m )
    {
        const double tolerance = 1.e-15;
        const double detM = determinant( m );
        FTL_VERIFY_DESCRIPTIVE( std::fabs( detM ) > tolerance, " detM = %g \n ", detM );

        // use a copy here
        const eigenX::MatrixSd< 3 > t( m );
        // main-diagonal
        m( 0, 0 ) = ( t( 1, 1 ) * t( 2, 2 ) - t( 1, 2 ) * t( 2, 1 ) ) / detM;
        m( 1, 1 ) = ( t( 0, 0 ) * t( 2, 2 ) - t( 0, 2 ) * t( 2, 0 ) ) / detM;
        m( 2, 2 ) = ( t( 0, 0 ) * t( 1, 1 ) - t( 0, 1 ) * t( 1, 0 ) ) / detM;
        // off-diagonal
        m( 0, 1 ) = ( t( 0, 2 ) * t( 2, 1 ) - t( 0, 1 ) * t( 2, 2 ) ) / detM;
        m( 0, 2 ) = ( t( 0, 1 ) * t( 1, 2 ) - t( 0, 2 ) * t( 1, 1 ) ) / detM;
        m( 1, 0 ) = ( t( 1, 2 ) * t( 2, 0 ) - t( 1, 0 ) * t( 2, 2 ) ) / detM;
        m( 1, 2 ) = ( t( 0, 2 ) * t( 1, 0 ) - t( 0, 0 ) * t( 1, 2 ) ) / detM;
        m( 2, 0 ) = ( t( 1, 0 ) * t( 2, 1 ) - t( 1, 1 ) * t( 2, 0 ) ) / detM;
        m( 2, 1 ) = ( t( 0, 1 ) * t( 2, 0 ) - t( 0, 0 ) * t( 2, 1 ) ) / detM;
        
        return detM;
    }
    //@}



    //--------------------------------------------------------------------------
    //! @name Computation of the contra-variant basis and Jacobian
    //@{

    /** Computation where the dimensions of the manifold is lower
     *  than the dimension of the embedding space.     
     *  \tparam DIM  Dimension of the embedding space
     *  \tparam LDIM Dimension of the embedded manifold     */
    template<unsigned DIM, unsigned LDIM>
    struct ContraVariantBasis
    {
        /** Compute the contra-variant basis and Jacobian via the metric tensor
         *  \param[in]  coVariant     Co-variant basis: 
         *                            \f$ g_\alpha [i] = d x_i / d \xi_\alpha \f$
         *  \param[out] contraVariant Contra-variant basis:
         *                            \f$ g^\alpha [i] = d \xi_\alpha / d x_i \f$
         *  \return     Jacobian
         */
        static double compute( const eigenX::MatrixSd<DIM,LDIM> & coVariant,
                               eigenX::MatrixSd<DIM,LDIM> & contraVariant )
        {
            //! - Co-variant metric coefficients:
            //!   \f$ G_{\alpha\beta} = <g_\alpha, g_\beta> \f$
            eigenX::MatrixSd<LDIM,LDIM> metricTens
                = coVariant.transpose( ) * coVariant;

            //! - Inverse of the basis (inplace): 
            //!   \f$ G^{\alpha\beta} = (G_{\alpha\beta})^{-1}\f$
            //!
            //! - Determinant (metric) is computed at the same time:
            //!   \f$ det J = \sqrt{ | (G_{\alpha\beta}) |} \f$
            const double detJ = std::sqrt( corlib::inverse( metricTens ) );
        
            //! - Contra-variant basis: 
            //!   \f$ g^\alpha = G^{\alpha\beta} g_\beta\f$
            contraVariant = coVariant * metricTens;

            return detJ;
        }
    };

    /** Computation with equal dimensions of manifold and embedding space
     *  \tparam DIM  Spatial dimension     */
    template<unsigned DIM>
    struct ContraVariantBasis<DIM,DIM>
    {
        /** Compute the co-variant basis by simple inversion
         *  \param[in]  coVariant     Co-variant basis: 
         *                            \f$ g_\alpha [i] = d x_i / d \xi_\alpha \f$
         *  \param[out] contraVariant Contra-variant basis:
         *                            \f$ g^\alpha [i] = d \xi_\alpha / d x_i \f$
         *  \return     Jacobian
         */
        static double compute( const eigenX::MatrixSd< DIM > & coVariant,
                               eigenX::MatrixSd< DIM > & contraVariant )
        {
            //! - Co-variante basis gives a square matrix, simply transpose and invert:
            //!   \f$  (g^\alpha [i] ) = (g_\alpha [i])^{-T} \f$
            //!
            //! - the determinant is computed at the same time:
            //!   \f$ det J = | (g_\alpha [ i ]) | \f$
            contraVariant = coVariant.transpose( );
            const double detJ = corlib::inverse( contraVariant );

            return detJ;
        }
    };
    //@}

    //------------------------------------------------------------------------------
    //! Compute the trace of a bounded matrix
    template< unsigned DIM >
    double trace( const eigenX::MatrixSd< DIM > & m )
    {
        return m.trace( );
    }

    //--------------------------------------------------------------------------
    //! Sum matrix entries
    template< typename MATRIX >
    double sumOfEntries( const MATRIX & mat )
    {
        double result = 0;
        for ( unsigned i = 0; i < mat.rows(); i ++ ) {
            for ( unsigned j = 0; j < mat.cols(); j ++ ) {
                result += mat(i,j);
            }
        }
        return result;
    }

    //--------------------------------------------------------------------------
    //! Inner product (double contraction) of two matrices
    //!
    //! \tparam     MATRIX    Matrix type
    //!
    //! \param[in]  a         First input matrix \f$[a_{ij}]\f$
    //! \param[in]  b         Second input matrix  \f$[b_{ij}]\f$
    //! \return               Inner product \f$\sum_i \sum_j ( a_{ij} b_{ij} )\f$
    template< typename MATRIX >
    double inner_prod( const MATRIX & a,
                       const MATRIX & b )
    {
        const unsigned rows = a.rows();
        const unsigned cols = a.cols();
        double c = 0.0;
        for ( unsigned i=0; i<rows; ++i )
            for ( unsigned j=0; j<cols; ++j )
                c += a( i, j ) * b( i, j );
        return c;
    }


    //------------------------------------------------------------------------------
    //! cross-product of two 3-dimensional vectors: c = a x b
    eigenX::VectorSd< 3 > cross_prod( const eigenX::VectorSd< 3 > & a,
                                      const eigenX::VectorSd< 3 > & b )
    {
        return a.cross( b );
    }

    inline
    eigenX::VectorSd< 3 > cross_prod( const eigenX::MatrixSd<3,2> & tngt )
    {
        eigenX::VectorSd< 3 > result;
        result(0) = tngt( 1, 0 ) * tngt( 2, 1 ) - tngt( 2, 0 ) * tngt( 1, 1 );
        result(1) = tngt( 2, 0 ) * tngt( 0, 1 ) - tngt( 0, 0 ) * tngt( 2, 1 );
        result(2) = tngt( 0, 0 ) * tngt( 1, 1 ) - tngt( 1, 0 ) * tngt( 0, 1 );
        return result;
    }

    inline
    eigenX::VectorSd< 2 > cross_prod( const eigenX::MatrixSd<2,1> & tngt )
    {
        eigenX::VectorSd< 2 > result;
        result(0) =  tngt(1,0);
        result(1) = -tngt(0,0);
        return result;
    }

    /*
    inline
    eigenX::VectorSd< 2 > cross_prod( const eigenX::VectorSd< 2 > & tngt )
    {
        eigenX::VectorSd< 2 > result;
        result(0) =  tngt(1);
        result(1) = -tngt(0);
        return result;
    }
	*/

    //--------------------------------------------------------------------------
    //! Cross-product matrix of a vector
    //!
    //! The cross-product matrix is defined as
    //!\f[
    //!    \vec{a} \times \vec{b} = A(\vec{a}) \cdot \vec{b}
    //!    \quad\mbox{with}\quad
    //!    A(\vec{a}) = [\vec{a}]_\times
    //!               = \left[\begin{array}{ccc}
    //!                       0 & -a_z & a_y \\ a_z & 0 & -a_x \\ -a_y & a_x & 0
    //!                 \end{array}\right]
    //!\f]
    //! The cross-product matrix is skew-symmetric.
    //!
    //! \param[in]  a   Vector \f$\vec{a}\f$
    //! \return         Cross-product matrix \f$[\vec{a}]_\times\f$
    eigenX::MatrixSd< 3 > crossProdMatrix( const eigenX::VectorSd< 3 > & a )
    {
        eigenX::MatrixSd< 3 > m;
        m( 0, 0 ) =      0.;    m( 0, 1 ) = -a( 2 );    m( 0, 2 ) =  a( 1 );
        m( 1, 0 ) =  a( 2 );    m( 1, 1 ) =      0.;    m( 1, 2 ) = -a( 0 );
        m( 2, 0 ) = -a( 1 );    m( 2, 1 ) =  a( 0 );    m( 2, 2 ) =      0.;
        return m;
    }

    //--------------------------------------------------------------------------
    //! Rotation matrix of rotation around axis with an angle (Rodrigues' formula)
    //!
    //! The rotation matrix around a axis vector \f$\vec{n}=\vec{a}/|\vec{a}|\f$
    //! with angle \f$\theta\f$ is given as
    //!\f[
    //!    R(\vec{n},\theta) = I + \sin(\theta) \, [\vec{n}]_\times
    //!                      + (1-\cos\theta) \, (\vec{n}\otimes\vec{n} - I)
    //!\f]
    //! with identity matrix \f$I\f$ and cross-product matrix \f$[\vec{n}]_\times\f$.
    //!
    //! For details, see for instance,
    //! http://en.wikipedia.org/wiki/Rodrigues%27_rotation_formula
    //!
    //! \param[in]   axis    Rotation axis vector \f$\vec{a}\f$
    //! \param[in]   angle   Rotation angle \f$\theta\f$
    //! \return              Rotation matrix \f$R(\vec{n},\theta)\f$
    eigenX::MatrixSd< 3 > rotationMatrix( const eigenX::VectorSd< 3 > & axis,
                                          const double & angle )
    {
        const eigenX::VectorSd< 3 > n { axis * ( 1. / axis.norm( ) ) };
        const eigenX::MatrixSd< 3 > rot {
            ( 1. - std::cos( angle ) ) * ( n * n.transpose( ) ) +
            std::cos( angle ) * eigenX::MatrixSd< 3 >::Identity( ) +
            std::sin( angle ) * corlib::crossProdMatrix( n ) };
        return rot;
    }

    //--------------------------------------------------------------------------
    //! Time-derivative of rotation matrix around axis with an angle
    //!
    //! The time-derivative of the rotation matrix around
    //! an axis unit-vector \f$\vec{n}(t)\f$ with angle \f$\theta(t)\f$ is given as
    //! \f[
    //!    \dot{R} = \frac{d R(\vec{n},\theta)}{d t}
    //!        = \dot{\theta} \Big(
    //!             \cos(\theta) \, [\vec{n}]_\times
    //!             \sin(\theta) \, (\vec{n}\otimes\vec{n} - I)
    //!        \Big)
    //!        + \sin(\theta) \, [\dot{\vec{n}}]_\times
    //!        + (1-\cos\theta) \, \big(\dot{\vec{n}}\otimes\vec{n} + \vec{n}\otimes\dot{\vec{n}}\big)
    //! \f]
    //! with identity matrix \f$I\f$ and cross-product matrix \f$[\vec{n}]_\times\f$.
    //!
    //! \param[in]   axis        Rotation axis unit-vector \f$\vec{n}\f$
    //! \param[in]   angle       Rotation angle \f$\theta\f$
    //! \param[in]   axisDeriv   Time-derivative of unit-axis vector \f$\dot{\vec{n}}\f$
    //! \param[in]   angleDeriv  Time-derivative rotation angle \f$\dot{\theta}\f$
    //! \return                  Time-derivative of rotation matrix \f$\dot{R}(\vec{n},\theta)\f$
    eigenX::MatrixSd< 3 > rotationMatrixDeriv( const eigenX::VectorSd< 3 > & axis,
                                               const double & angle,
                                               const eigenX::VectorSd< 3 > & axisDeriv,
                                               const double & angleDeriv )
    {
        const eigenX::VectorSd< 3 > n { axis };
        const eigenX::MatrixSd< 3 > rotDeriv {
            -std::sin( angle ) * eigenX::MatrixSd< 3 >::Identity( ) * angleDeriv +
            std::cos( angle ) * corlib::crossProdMatrix( n )  * angleDeriv +
            std::sin( angle ) * ( n * n.transpose( ) ) * angleDeriv +
            std::sin( angle ) * corlib::crossProdMatrix( axisDeriv ) +
            ( 1. - std::cos( angle ) ) * ( n * axisDeriv.transpose( ) ) +
            ( 1. - std::cos( angle ) ) * ( axisDeriv * n.transpose( ) ) };
        return rotDeriv;
    }

    //--------------------------------------------------------------------------
    //! Metric of small rectangular matrices (up to 3x3)
    //!
    //! These methods call the determinant function in case of
    //! <i>square</i> matrices. In the case of <i>non-square</i> matrices the
    //! metric is determined by a vector norm or cross product.
    double metric( const eigenX::MatrixSd< 1 > & m )
    {
        return determinant( m );
    }

    double metric( const eigenX::MatrixSd< 2 > & m )
    {
        return determinant( m );
    }

    double metric( const eigenX::MatrixSd< 3 > & m )
    {
        return determinant( m );
    }

    double metric( const eigenX::MatrixSd< 2, 1 > & m )
    {
        return m.col( 0 ).norm( );
    }

    double metric( const Eigen::Matrix< double, 3, 1 > & m )
    {
        return m.col( 0 ).norm( );
    }
    
    double metric( const Eigen::Matrix< double, 3, 2 > & m )
    {
        const eigenX::VectorSd< 3 > n =
            m.col( 0 ).cross( m.col( 1 ) );
        return n.norm( );
    }

    //------------------------------------------------------------------------------
    //! write matrix in CSV format to stream
    template< typename MAT >
    inline std :: ostream & writeCSVmatrix( std :: ostream & out, 
                                            const MAT & mat )
    {
        out << std::setprecision( 6 );
        for ( unsigned m = 0; m < mat.rows( ); m ++ ) {
            for ( unsigned n = 0; n < mat.cols( ); n ++ ) {
                out << mat( m, n ) << "  ";
            }
            out << std::endl;
        }
        return out;
    }

    //------------------------------------------------------------------------------
    //! write vector 
    template< typename VEC >
    inline std :: ostream & writeVector( std :: ostream & out, 
                                         const VEC & vec,
                                         const std::string spacer = "\n" )
    {
        out << std :: setprecision( 6 );
        for ( unsigned m = 0; m < vec.size( ); m ++ ) {
            out << vec[ m ] << spacer;
        }
        return out;
    }

    //------------------------------------------------------------------------------
    //! Cholesky decomposition (copied from: http://www.guwi17.de/ublas/examples/ )
    /** \brief decompose the symmetric positive definit matrix A into product L L^T.
     *
     * \param A square symmetric positive definite input matrix (only the lower 
     *        triangle is accessed)
     * \param L lower triangular output matrix 
     * \return nonzero if decompositon fails (the value is 1 + the number of 
     *         the failing row)
     */
    template< typename SPARSEMAT, typename LOWERMAT >
    inline size_t cholesky_decompose(const SPARSEMAT & A, LOWERMAT & L)
    {

    	L = A.llt( ).matrixL( );
        return 0;      
    }


    //------------------------------------------------------------------------------
    /** \brief solve system L L^T x = b inplace using a Cholesky decomposition
     *  \details copied from: http://www.guwi17.de/ublas/examples/ 
     *  \param L a triangular matrix
     *  \param x input: right hand side b; output: solution x
     */
    template< typename LOWERMAT, typename SPARSEVEC >
    inline void cholesky_solve(const LOWERMAT & L, SPARSEVEC & x )
    {
    	SPARSEVEC b( x );
    	x.setZero( );
    	x = L.llt().solve( b );
    }

    //==========================================================================
    // General fixed-size Matrix-Matrix product implementation

    //! \cond SKIPDOX
    namespace detail_{

        //----------------------------------------------------------------------
        // structure to switch between a+=b and a=b
        template<bool YES> struct AddTo;

        template<> 
        struct AddTo<false> // only assign b to a
        { 
            static void apply( double & a, const double b ) { a  = b; }
        };

        template<> 
        struct AddTo<true> // add b to a
        {
            static void apply( double & a, const double b ) { a += b; }
        };

        //----------------------------------------------------------------------
        // structure to help for transposed matrix access
        template<unsigned M, unsigned N, bool YES> struct TransIndex;

        // no transposition
        template<unsigned M, unsigned N>
        struct TransIndex<M,N,false>
        {
            static const unsigned numFirst  = M;
            static const unsigned numSecond = N;

            static unsigned getFirst(  const unsigned row, 
                                       const unsigned col ) { return row; }
            static unsigned getSecond( const unsigned row, 
                                       const unsigned col ) { return col; }
        };

        // transposition
        template<unsigned M, unsigned N>
        struct TransIndex<M,N,true>
        {
            static const unsigned numFirst  = N; 
            static const unsigned numSecond = M;
            
            static unsigned getFirst(  const unsigned row, 
                                       const unsigned col ) { return col; }
            static unsigned getSecond( const unsigned row, 
                                       const unsigned col ) { return row; }
            
        };
    }
    //! \endcond

    //------------------------------------------------------------------------------
    /** \brief Matrix-matrix multiplication with optional transposition
     *  \details The compute function of this object multiplies two matrices A and B
     *  and stores the result in D. Both operands can be transposed and the storage
     *  in D can be assignment or addition, i.e.
     *  \f[
     *         D =  op(A) * op(B)  \quad and \quad D+= op(A) * op(B)
     *  \f]
     *  with op = trans or op = identity are implemented.
     *  \tparam  MA      Number of rows    of A
     *  \tparam  NA      Number of columns of A
     *  \tparam  MB      Number of rows    of B
     *  \tparam  NB      Number of columns of B
     *  \tparam  TRANSA  True if A should be transposed
     *  \tparam  TRANSB  True if B should be transposed
     *  \tparam  ADDTO   True if the result should be added to D
     */
    template<unsigned MA, unsigned NA, unsigned MB, unsigned NB,
             bool TRANSA = false, bool TRANSB = false, 
             bool ADDTO = false>
    struct Gemm
    {
        //! Helper for transposition
        typedef detail_::TransIndex<MA,NA,TRANSA> TransA;
        typedef detail_::TransIndex<MB,NB,TRANSB> TransB;

        //! Sanity check 
        FTL_STATIC_ASSERT_MSG( ( TransA::numSecond == TransB::numFirst ), 
                               "Matrix sizes do not fit" );

        //! Size of inner products
        static const unsigned N = TransA::numSecond;

        //! Function to multiply 'a' and 'b' and store the result in 'd'
        static void compute( const eigenX::MatrixSd<MA,NA> & a,
                             const eigenX::MatrixSd<MB,NB> & b,
                             eigenX::MatrixSd<TransA::numFirst,
                                              TransB::numSecond> & d )
        {
            for ( unsigned r = 0; r < TransA::numFirst; r ++ ) {
                for ( unsigned c = 0; c < TransB::numSecond; c ++ ) {
                    double innerProd = 0.;
                    // do inner product
                    for ( unsigned k = 0; k < N; k ++ ) {
                        // Compute  sum_k( a(r,k) * c(k,c) )
                        // where the indices are swapped according to 
                        // desired transposition
                        innerProd += 
                            a( TransA::getFirst(  r, k), TransA::getSecond( r, k) ) *
                            b( TransB::getSecond( c, k), TransB::getFirst(  c, k) );

                    }
                                   
                    // storage of result
                    detail_::AddTo<ADDTO>::apply( d(r,c), innerProd ); 
                }
            }
        }
    };

    //--------------------------------------------------------------------------
    //! Convenience function for C = A * B
    template<typename M1, typename M2, typename M3>
    void aMatProdB( const M1 & a, const M2 & b, M3 & c ) 
    {
        Gemm<M1::max_size1, M1::max_size2, 
             M2::max_size1, M2::max_size2>::compute( a, b, c);
    }


    //--------------------------------------------------------------------------
    //! Convenience function for C = A * B'
    template<typename M1, typename M2, typename M3>
    void aMatProdBT( const M1 & a, const M2 & b, M3 & c ) 
    {
        Gemm<M1::max_size1, M1::max_size2, 
             M2::max_size1, M2::max_size2, false, true>::compute( a, b, c);
    }

    //--------------------------------------------------------------------------
    //! Convenience function for C = A' * B
    template<typename M1, typename M2, typename M3>
    void aTMatProdB( const M1 & a, const M2 & b, M3 & c ) 
    {
        Gemm<M1::max_size1, M1::max_size2, 
             M2::max_size1, M2::max_size2, true>::compute( a, b, c);
    }

}//namespace corlib




#endif



