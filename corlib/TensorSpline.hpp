// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//
//                 Computational Structural Mechanics Lab
//                         University of Cambridge
//
// This software is copyrighted and all rights are retained by the University of
// Cambridge. All use, disclosure, and/or reproduction of any part not expressly
// authorised by the University of Cambridge is prohibited. (C) 2010
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

//! @file TensorSpline.hpp
//! @todo   make the calls recursive

#ifndef corlib_tensorspline_h
#define corlib_tensorspline_h
//------------------------------------------------------------------------------
//! Eigen includes
#include <Eigen/Core>
//! corlib includes
#include <corlib/eigenX.hpp>
#include <corlib/BSpline.hpp>
#include <corlib/Shape.hpp>
#include <corlib/ShapefunTraits.hpp>
//------------------------------------------------------------------------------
namespace corlib{
    template<unsigned DEGREE, unsigned DIM> class TensorSpline;
    namespace eigenX = corlib::eigenX;
}

//------------------------------------------------------------------------------
namespace corlib{
        
    //--------------------------------------------------------------------------
    /** \brief Evaluate one-dimensional splines as tensor-products (DIM=1)
     *  \details 
     *  <pre>
     *
     *   0----------1 --> x 
     *
     *  </pre>
     */
    template<unsigned DEGREE>
    class TensorSpline<DEGREE,1>
    {
    public:
        static const unsigned      degree       = DEGREE;
        static const unsigned      numFunctions = degree + 1;
        static const unsigned      dim          = 1;
        static const corlib::shape myShape      = corlib::LINE;

        static const corlib::SfunClass sfc      = TENSORSPLINE;

        typedef eigenX::VectorSd<dim>                VecDim;
        typedef eigenX::VectorSd<numFunctions>       VecNFun;
        typedef eigenX::MatrixSd<dim, numFunctions>  MatDimNFun;
        typedef Eigen::Matrix<VecNFun, dim, dim>     MatVecNFunDimDim;

        typedef corlib::BSpline<degree>  BSpline;

    public:
        //! Evaluate the underlying 1D-spline
        void evaluate( const VecDim & xi, VecNFun & phi ) const
        {
            spline_.evaluate( xi(0), phi );
            return;
        }

        //! Evaluate the underlying 1D-spline's gradient
        void evaluateGradient( const VecDim & xi, MatDimNFun & dphi ) const
        {
            typename BSpline::VecNspl dphi1 = dphi.row( 0 );
            spline_.evaluateGradient( xi(0), dphi1 );
            dphi.row( 0 ) = dphi1;
            return;
        }

        //! Evaluate the underlying 1D-spline's hessian
        void evaluateHessian( const VecDim & xi, MatVecNFunDimDim & ddphi ) const
        {                
            typename BSpline::VecNspl ddphi1;
            spline_.evaluateHessian( xi(0), ddphi1 );
               
            for ( unsigned j = 0; j < numFunctions; j ++ ) {
                ddphi(0,0)(j) = ddphi1(j);
            }
            return;
        }

    private:
        BSpline spline_; //!< underlying 1D-spline
    };

    //--------------------------------------------------------------------------
    /** \brief Evaluate one-dimensional splines as tensor-products (DIM=2)
     *  \details
     *  <pre>
     *
     *            2-----------3          
     *            |           |          
     *            |           |          
     *            |           |          
     *            |           |          
     *            |           |          
     *            0-----------1          
     *  </pre>                          
     *
     */
    template<unsigned DEGREE>
    class TensorSpline<DEGREE,2>
    {
    public:
        static const unsigned     degree         = DEGREE;
        static const unsigned     numFunctions1D = degree+1;
        static const unsigned     numFunctions   = numFunctions1D * numFunctions1D;
        static const unsigned     dim            = 2;
        static const corlib::shape myShape       = corlib::PIXEL;

        static const corlib::SfunClass sfc       = TENSORSPLINE;

        typedef eigenX::VectorSd<dim>                VecDim;
        typedef eigenX::VectorSd<numFunctions>       VecNFun;
        typedef eigenX::MatrixSd<dim, numFunctions>  MatDimNFun;
        typedef Eigen::Matrix<VecNFun, dim, dim>     MatVecNFunDimDim;
    private:
        typedef corlib::BSpline<degree>  BSpline;

    public:
        //! Evaluate underlying 1D-spline and assemble it to
        //! a 2-variate spline by using the tensor-product.
        //!
        //! Ordering of spline basis functions with n=degree
        //! in vector phi
        //!\verbatim
        //!     ^ xi_1
        //!     |
        //!     |
        //!
        //! n*(n+1)+0...n*(n+1)+1...n*(n+1)+2...  ..n*(n+1)+n
        //!     :           :           :              :
        //!
        //!     :           :           :              :
        //! 1*(n+1)+0...1*(n+1)+1...1*(n+1)+2...  ..1*(n+1)+n
        //!     :           :           :              :
        //!     0...........1...........2.......  .....n       ---> xi_0
        //!\endverbatim
        void evaluate( const VecDim & xi, VecNFun & phi ) const
        {
            // arrays for one-dimensional spline evaluation
            typename BSpline::VecNspl phi1, phi2;
            // evluate one-dimensional spline at coodinate directions
            spline_.evaluate( xi(0), phi1 ); 
            spline_.evaluate( xi(1), phi2 ); 
            // insert tensor product into array
            unsigned offSet = 0;
            for ( unsigned j = 0; j < numFunctions1D; j ++ ) {
                for ( unsigned i = 0; i < numFunctions1D; i ++ ) {
                    phi( offSet++ ) = phi1(i) * phi2(j);
                }
            }
            return;
        }

        //----------------------------------------------------------------------
        /** Evaluate underlying 1D-spline's derivatives and perform tensor 
         *  product. The output array is of dimensions dim x numFunctions.
         *  In this array, the column number refers to the number of the 
         *  considered tensor-product spline as documented above in #evaluate.
         *  The row refers to the direction of derivative. Hence, 
         *  \f[
         *       dphi[i,n] = \frac{\partial B_n(\xi_1,\xi_2)}{\partial \xi_i}
         *  \f]
         *  \param[in]  xi   Local coordinate of evaluation
         *  \param[out] dphi Array of derivatives
         */
        void evaluateGradient( const VecDim & xi, MatDimNFun & dphi ) const
        {
            // arrays for one-dimensional spline evaluation
            typename BSpline::VecNspl phi1, phi2, dphi1, dphi2;
            // evaluate one-dimensional spline at coodinate directions
            spline_.evaluate( xi(0), phi1 ); 
            spline_.evaluate( xi(1), phi2 ); 
            // evaluate gradient of one-dimensional spline at coodinate directions
            spline_.evaluateGradient( xi(0), dphi1 ); 
            spline_.evaluateGradient( xi(1), dphi2 ); 
            // insert tensor product into result array
            unsigned offSet = 0;
            for ( unsigned j = 0; j < numFunctions1D; j ++ ) {
                for ( unsigned i = 0; i < numFunctions1D; i ++ ) {
                    dphi( 0, offSet ) = dphi1( i ) *  phi2( j );
                    dphi( 1, offSet ) =  phi1( i ) * dphi2( j );
                    offSet ++;
                }
            }
            return;
        }

        //----------------------------------------------------------------------
        /** Evaluate the underlying 1D-spline's hessian
         *  Similar to the above, the underlying one-dimensional b-spline's
         *  derivatives are evaluated at given local coordinate and arranged in
         *  the output array. This array has dimensions dim x dim x numFunctions.
         *  Its entries are 
         *  \f[
         *    ddphi[i,j,n] = 
         *     \frac{\partial^2 B_n(\xi_1,\xi_2)}{\partial \xi_i \partial \xi_j}
         *  \f]
         *  \param[in]  xi    Local coordinate of evaluation
         *  \param[out] ddphi Array of derivatives
         */
        void evaluateHessian (const VecDim & xi, MatVecNFunDimDim & ddphi ) const
        {
            // arrays for one-dimensional spline evaluation
            typename BSpline::VecNspl phi1, phi2, dphi1, dphi2, ddphi1, ddphi2;
            // evaluate one-dimensional spline at coodinate directions
            spline_.evaluate( xi(0), phi1 ); 
            spline_.evaluate( xi(1), phi2 ); 
            // evaluate gradient of one-dimensional spline at coodinate directions
            spline_.evaluateGradient( xi(0), dphi1 ); 
            spline_.evaluateGradient( xi(1), dphi2 ); 
            // evaluate hessian of one-dimensional spline at coodinate directions
            spline_.evaluateHessian( xi(0), ddphi1 ); 
            spline_.evaluateHessian( xi(1), ddphi2 ); 

            unsigned offSet = 0;
            for ( unsigned j = 0; j < numFunctions1D; j ++ ) {
                for ( unsigned i = 0; i < numFunctions1D; i ++ ) {
                    ddphi(0, 0)(offSet) = ddphi1(i)*phi2(j); 
                    ddphi(0, 1)(offSet) = ddphi(1, 0)(offSet) = dphi1(i)*dphi2(j);
                    ddphi(1, 1)(offSet) = phi1(i)*ddphi2(j); 
                    offSet++;
                }
            }
            return;
        }

    private:
        BSpline spline_; //!< underlying 1D-spline
    };

    //--------------------------------------------------------------------------
    /** \brief Evaluate one-dimensional splines as tensor-products (DIM=3)
     *  \details
     *  <pre>
     *              A z
     *              |
     *              |
     *              4----------6
     *              |\         |\
     *              | \        | \
     *              |  \       |  \
     *              |   5----------7
     *              |   |      |   |
     *              0---|------2 --|----> y
     *               \  |       \  |
     *                \ |        \ |
     *                 \|         \|
     *                  1----------3
     *                   \
     *                    \
     *                    `'x
     *  </pre>
     */
    template<unsigned DEGREE>
    class TensorSpline<DEGREE,3>
    {
    public:
        static const unsigned     degree         = DEGREE;
        static const unsigned     numFunctions1D = degree+1;
        static const unsigned     numFunctions   = numFunctions1D * numFunctions1D * numFunctions1D;
        static const unsigned     dim            = 3;
        static const corlib::shape myShape       = corlib::VOXEL;

        static const corlib::SfunClass sfc       = TENSORSPLINE;

        typedef eigenX::VectorSd<dim>               VecDim;
        typedef eigenX::VectorSd<numFunctions>      VecNFun;
        typedef eigenX::MatrixSd<dim, numFunctions> MatDimNFun;
        typedef Eigen::Matrix<VecNFun, dim, dim>    MatVecNFunDimDim;

    private:
        typedef corlib::BSpline<degree>  BSpline;

    public:
        //! Evaluate underlying 1D spline and generate tensor-product result
        //! See corlib::TensorSpline<degree,2>#evaluate
        void evaluate( const VecDim & xi, VecNFun & phi ) const
        {
            // arrays for one-dimensional spline evaluation
            typename BSpline::VecNspl phi1, phi2, phi3;
            // evluate one-dimensional spline at coodinate directions
            spline_.evaluate( xi(0), phi1 ); 
            spline_.evaluate( xi(1), phi2 ); 
            spline_.evaluate( xi(2), phi3 ); 
            // insert tensor product into array
            unsigned offSet = 0;
            for ( unsigned k = 0; k < numFunctions1D; k ++ ) {
                for ( unsigned j = 0; j < numFunctions1D; j ++ ) {
                    for ( unsigned i = 0; i < numFunctions1D; i ++ ) {
                        phi( offSet++ ) = phi1(i) * phi2(j) * phi3(k);
                    }
                }
            }
            return;
        }

        //! Evaluate underlying 1D spline and derivatives and combine to tensor-product
        //! See corlib::TensorSpline<degree,2>#evaluateGradient
        void evaluateGradient( const VecDim & xi, MatDimNFun & dphi ) const
        {
            // arrays for one-dimensional spline evaluation
            typename BSpline::VecNspl phi1, phi2, phi3, dphi1, dphi2, dphi3;
            // evaluate one-dimensional spline at coodinate directions
            spline_.evaluate( xi(0), phi1 ); 
            spline_.evaluate( xi(1), phi2 ); 
            spline_.evaluate( xi(2), phi3 ); 
            // evaluate gradient of one-dimensional spline at coodinate directions
            spline_.evaluateGradient( xi(0), dphi1 ); 
            spline_.evaluateGradient( xi(1), dphi2 ); 
            spline_.evaluateGradient( xi(2), dphi3 ); 
            // insert tensor product into result array
            unsigned offSet = 0;
            for ( unsigned k = 0; k < numFunctions1D; k ++ ) {
                for ( unsigned j = 0; j < numFunctions1D; j ++ ) {
                    for ( unsigned i = 0; i < numFunctions1D; i ++ ) {
                        dphi( 0, offSet ) = dphi1( i ) *  phi2( j ) *  phi3( k );
                        dphi( 1, offSet ) =  phi1( i ) * dphi2( j ) *  phi3( k );
                        dphi( 2, offSet ) =  phi1( i ) *  phi2( j ) * dphi3( k );
                        offSet ++;
                    }
                }
            }
            return;
        }

        //! Evaluate the underlying 1D-spline's hessian
        //! See corlib::TensorSpline<degree,2>#evaluateHessian
        void evaluateHessian( const VecDim & xi, MatVecNFunDimDim & ddphi ) const
        {
            // arrays for one-dimensional spline evaluation
            typename BSpline::VecNspl   phi1,   phi2,   phi3;
            typename BSpline::VecNspl  dphi1,  dphi2,  dphi3;
            typename BSpline::VecNspl ddphi1, ddphi2, ddphi3;
            // evaluate one-dimensional spline at coodinate directions
            spline_.evaluate( xi(0), phi1 ); 
            spline_.evaluate( xi(1), phi2 ); 
            spline_.evaluate( xi(2), phi3 ); 
            // evaluate gradient of one-dimensional spline at coodinate directions
            spline_.evaluateGradient( xi(0), dphi1 ); 
            spline_.evaluateGradient( xi(1), dphi2 ); 
            spline_.evaluateGradient( xi(2), dphi3 ); 
            // evaluate hessian of one-dimensional spline at coodinate directions
            spline_.evaluateHessian( xi(0), ddphi1 ); 
            spline_.evaluateHessian( xi(1), ddphi2 ); 
            spline_.evaluateHessian( xi(2), ddphi3 ); 

            unsigned offSet = 0;
            for ( unsigned k = 0; k < numFunctions1D; k ++ ) {
                for ( unsigned j = 0; j < numFunctions1D; j ++ ) {
                    for ( unsigned i = 0; i < numFunctions1D; i ++ ) {
                        ddphi(0, 0)(offSet) = ddphi1(i) *   phi2(j) *   phi3(k);
                        ddphi(0, 1)(offSet) =  dphi1(i) *  dphi2(j) *   phi3(k);
                        ddphi(0, 2)(offSet) =  dphi1(i) *   phi2(j) *  dphi3(k);
                        ddphi(1, 0)(offSet) = ddphi(0, 1)(offSet);  // sym
                        ddphi(1, 1)(offSet) =   phi1(i) * ddphi2(j) *   phi3(k);
                        ddphi(1, 2)(offSet) =   phi1(i) *  dphi2(j) *  dphi3(k);
                        ddphi(2, 0)(offSet) = ddphi(0, 2)(offSet);  // sym
                        ddphi(2, 1)(offSet) = ddphi(1, 2)(offSet);  // sym
                        ddphi(2, 2)(offSet) =   phi1(i) *   phi2(j) * ddphi3(k);
                        offSet++;
                    }
                }
            }
            return;
        }

    private:
        BSpline spline_; //! Underlying B-spline
    };
    //--------------------------------------------------------------------------

}

//------------------------------------------------------------------------------
#endif
