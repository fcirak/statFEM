// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//
//                 Computational Structural Mechanics Lab
//                         University of Cambridge
//
// This software is copyrighted and all rights are retained by the University of
// Cambridge. All use, disclosure, and/or reproduction of any part not expressly
// authorised by the University of Cambridge is prohibited. (C) 2010
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

//! @file SimplexQuadrature.hpp

#ifndef corlib_simplexquadrature_h
#define corlib_simplexquadrature_h
//------------------------------------------------------------------------------
#include <utility>
#include <array>

#include <corlib/GaussLegendre.hpp>
#include <corlib/TensorQuadrature.hpp>
#include <corlib/misc.hpp>
#include <corlib/eigenX.hpp>

namespace corlib{
    template<unsigned DIM, unsigned QPS> class SimplexQuadrature;

    namespace detail_{
        //! Helper for the construction of the quadrature object
        template<unsigned DIM, unsigned QPS> 
        struct SimplexQuadratureConstructor;
    }

    template<unsigned DIM, unsigned DEGREE> struct SimplexQuadratureOrder;
}


//------------------------------------------------------------------------------
/** \brief %Quadrature rules for simplices (triangles and tetrahedra)
 *  \details Object generates and stores the quadrature points and weights for
 *  the numerical integration over an N-simplex (N=2 or N=3).
 *  <h4>References</h4>
 *  -  D.A. Dunavant,
 *  High degree efficient symmetrical Gaussian quadrature rules for the
 *  triangle, IJNME, 1985 (21)
 *  - P. Keast
 *  Moderate-degree tetrahedral quadrature formulas, CMAME, 1986 (55)
 *  \tparam DIM Spatial dimension (2 or 3)
 *  \tparam QPS Number of quadrature points
 */
template<unsigned DIM, unsigned QPS>
class corlib::SimplexQuadrature
{
public:
    static const unsigned dim         = DIM;
    static const unsigned numPoints   = QPS;
    static const corlib::shape theShape = corlib::Simplex<dim>::theShape;

    typedef eigenX::VectorSd<dim>                                                  VecDim;
    typedef std::pair<double, VecDim>                                    PairDoubleVecDim;
    typedef typename std::array<PairDoubleVecDim, numPoints>::const_iterator     QuadIter;

    //! Begin of array iterator
    QuadIter begin() const  { return weightsAndPoints_.begin(); }
    //! End of array iterator
    QuadIter end()   const  { return weightsAndPoints_.end();   }

    //! Constructor which generates the array of (weight,point)-pairs
    SimplexQuadrature() 
    {
        detail_::SimplexQuadratureConstructor<dim,numPoints>()( weightsAndPoints_ );
    }

private:
    //! Array of weight and point pairs for the tensor product rule
    std::array< PairDoubleVecDim, numPoints > weightsAndPoints_;
};

//------------------------------------------------------------------------------
namespace corlib{
    namespace detail_{
        //! \cond SKIPDOX  // documentation of this ??

        template<unsigned QPS>
        struct SimplexQuadratureConstructor<0,QPS>
        {
            void operator()( std::array<std::pair<double, eigenX::VectorSd<0> >, QPS> &
                             weightsAndPoints )
            {
                return;
            }
        };


        //---------------------------------------------------------------------
        //! One-dimensional quadrature copies from the TensorQuadrature
        template<unsigned QPS>
        struct SimplexQuadratureConstructor<1,QPS>
        {
            void operator()( std::array<std::pair<double, eigenX::VectorSd<1> >, QPS> &
                             weightsAndPoints )
            {
                corlib::TensorQuadrature<corlib::GaussLegendre<QPS>,1> tensorQuadrature;
                std::copy( tensorQuadrature.begin(), tensorQuadrature.end(), weightsAndPoints.begin() );
                return;
            }
        };

        //----------------------------------------------------------------------
        //! @name Quadrature points and weights for triangle 
        //! \f$(\xi,\eta)\in[0,1]\times[0,1-\xi]\f$
        //@{

#if 0
        //! Determine quadrature points and weights on triangle
        //! \f$\{ (\xi_1,\xi_2) | 0\leq\xi_1,\xi_2\leq1 and \xi_1+\xi_2\leq1 \}\f$
        //! using the Gauss-Legendre weights on square
        //! \f$\{ (s_1,s_2) | 0 \leq s_1,s_2 \leq 1 \}\f$
        //! using the map
        //!\f[
        //!     \xi_1 = s_1  \quad and \quad \xi_2 = (1 - s_1) s_2
        //!\f]
        //! with the Jacobian determinant
        //!\f[
        //!     \left|\frac{\partial (\xi_1,\xi_2)}{\partial (s_1,s_2)}\right| = 1 - s_1
        //!\f]
        //! Therefore, every Gauss-Legendre point on the square
        //! is easily transformed to the triangle with the above relations.
        //!
        //! <h4>References</h4>
        //! -  J.C Lachat and J.O. Watson,
        //! Effective numerical treatment of boundary integral equations,
        //! IJNME, 1976(10)
        //! -  M.G. Duffy,
        //! Quadrature over a pyramid or cube of integrands with a singularity at a vertex,
        //! SIAM JNA, 1982 (21)
        //! - H.T. Rathod, K.V. Nagaraja, B. Venkatesudu and N.L. Ramesh,
        //! Gauss Legendre quadrature over a triangle, J. Indian Inst. Sci., 2004, vol 84, p 183-188.
        //! http://journal.library.iisc.ernet.in/vol200405/paper6/rathod.pdf
        template<unsigned QPS>
        struct SimplexQuadratureConstructor<2,QPS>
        {
            // compute the integral part of the DIM-th root of QPS
            static const unsigned numPoints1D = corlib::NthRoot<QPS,2>::value;
            // total number of points is the previous result to the DIM-th power
            static const unsigned numPoints   = corlib::mToTheN<numPoints1D,2>::value;

            void operator()( std::array< std::pair<double, eigenX::VectorSd<2> >, numPoints > &
                             weightsAndPoints )
            {
                corlib::detail_::TensorQuadratureConstructor< corlib::GaussLegendre<numPoints1D>, 2 >()( weightsAndPoints );
                for ( unsigned i = 0; i < numPoints; ++i ) {
                    const double fact = 1. - weightsAndPoints[ i ].second[ 0 ];
                    weightsAndPoints[ i ].second[ 1 ] *= fact;
                    weightsAndPoints[ i ].first *= fact;
                }
            }
        };
#endif

        //! 1st order rule integrating linear polynomials exactly
        template<>
        struct SimplexQuadratureConstructor<2,1>
        {
            void operator()( std::array<std::pair<double,eigenX::VectorSd<2> >, 1> &
                             weightsAndPoints )
            {
                eigenX::VectorSd<2> point;
                point(0) = 0.333333333333333; point(1) = 0.333333333333333;
                weightsAndPoints[0] = std::make_pair( 0.5, point );
            }
        };

        //! 2nd order rule integrating quadratic (\f$\xi^a\eta^b\f$ with \f$a+b\leq2\f$) polynomials exactly
        template<>
        struct SimplexQuadratureConstructor<2,3>
        {
            void operator()( std::array<std::pair<double,eigenX::VectorSd<2> >, 3> &
                             weightsAndPoints )
            {
                eigenX::VectorSd<2> point;
                point(0) = 0.666666666666667; point(1) = 0.166666666666667;
                weightsAndPoints[0] = std::make_pair( 0.166666666666666, point );

                point(0) = 0.166666666666667; point(1) = 0.666666666666667;
                weightsAndPoints[1] = std::make_pair( 0.166666666666666, point );

                point(0) = 0.166666666666667; point(1) = 0.166666666666667;
                weightsAndPoints[2] = std::make_pair( 0.166666666666666, point );
            }
        };

        //! 3rd order rule integrating cubic (\f$\xi^a\eta^b\f$ with \f$a+b\leq3\f$) polynomials exactly
        template<>
        struct SimplexQuadratureConstructor<2,4>
        {
            void operator()( std::array<std::pair<double,eigenX::VectorSd<2> >, 4> &
                             weightsAndPoints )
            {
                eigenX::VectorSd<2> point;

                point(0) = 0.333333333333333; point(1) = 0.333333333333333;
                weightsAndPoints[0] = std::make_pair( -0.28125, point );

                point(0) = 0.6; point(1) = 0.2;
                weightsAndPoints[1] = std::make_pair( 0.260416666666667, point );

                point(0) = 0.2; point(1) = 0.6;
                weightsAndPoints[2] = std::make_pair( 0.260416666666667, point );

                point(0) = 0.2; point(1) = 0.2;
                weightsAndPoints[3] = std::make_pair( 0.260416666666667, point );
            }
        };

        //! 4th order rule integrating quartic (\f$\xi^a\eta^b\f$ with \f$a+b\leq4\f$) polynomials exactly
        template<>
        struct SimplexQuadratureConstructor<2,6>
        {
            void operator()( std::array<std::pair<double,eigenX::VectorSd<2> >, 6> &
                             weightsAndPoints )
            {
                eigenX::VectorSd<2> point;

                point(0) = 0.10810301816807;  point(1) = 0.445948490915965;
                weightsAndPoints[0] = std::make_pair( 0.111690794839005, point );

                point(0) = 0.816847572980459; point(1) = 0.091576213509771;
                weightsAndPoints[1] = std::make_pair( 0.054975871827661, point );

                point(0) = 0.445948490915965; point(1) = 0.10810301816807;
                weightsAndPoints[2] = std::make_pair( 0.111690794839005, point );

                point(0) = 0.445948490915965; point(1) = 0.445948490915965;
                weightsAndPoints[3] = std::make_pair( 0.111690794839005, point );

                point(0) = 0.091576213509771; point(1) = 0.816847572980459;
                weightsAndPoints[4] = std::make_pair( 0.054975871827661, point );

                point(0) = 0.091576213509771; point(1) = 0.091576213509771;
                weightsAndPoints[5] = std::make_pair( 0.054975871827661, point );
            }
        };

        //! 5th order rule integrating quintic (\f$\xi^a\eta^b\f$ with \f$a+b\leq5\f$) polynomials exactly
        template<>
        struct SimplexQuadratureConstructor<2,7>
        {
            void operator()( std::array<std::pair<double,eigenX::VectorSd<2> >, 7> &
                             weightsAndPoints )
            {
                eigenX::VectorSd<2> point;

                point(0) = 0.333333333333333; point(1) = 0.333333333333333;
                weightsAndPoints[0] = std::make_pair( 0.1125, point );

                point(0) = 0.05971587178977;  point(1) = 0.470142064105115;
                weightsAndPoints[1] = std::make_pair( 0.066197076394253, point );

                point(0) = 0.797426985353087; point(1) = 0.101286507323456;
                weightsAndPoints[2] = std::make_pair( 0.0629695902724135, point );

                point(0) = 0.470142064105115; point(1) = 0.05971587178977;
                weightsAndPoints[3] = std::make_pair( 0.066197076394253, point );

                point(0) = 0.470142064105115; point(1) = 0.470142064105115;
                weightsAndPoints[4] = std::make_pair( 0.066197076394253, point );

                point(0) = 0.101286507323456; point(1) = 0.797426985353087;
                weightsAndPoints[5] = std::make_pair( 0.0629695902724135, point );

                point(0) = 0.101286507323456; point(1) = 0.101286507323456;
                weightsAndPoints[6] = std::make_pair( 0.0629695902724135, point );
            }
        };

        //! 7th order rule integrating septic (\f$\xi^a\eta^b\f$ with \f$a+b\leq7\f$) polynomials exactly
        template<>
        struct SimplexQuadratureConstructor<2,13>
        {
            void operator()( std::array<std::pair<double,eigenX::VectorSd<2> >, 13> &
                             weightsAndPoints )
            {
                eigenX::VectorSd<2> point;

                const double a0 = 0.333333333333333;
                const double w0 = -0.074785022233835;

                point(0) = a0; point(1) = a0;
                weightsAndPoints[0] = std::make_pair( w0, point );

                const double a1 = 0.260345966079038;
                const double b1 = 0.479308067841923;
                const double w1 = 0.087807628716602;

                point(0) = a1;  point(1) = a1;
                weightsAndPoints[1] = std::make_pair( w1, point );

                point(0) = a1; point(1) = b1;
                weightsAndPoints[2] = std::make_pair( w1, point );

                point(0) = b1; point(1) = a1;
                weightsAndPoints[3] = std::make_pair( w1, point );

                const double a2 = 0.065130102902216;
                const double b2 = 0.869739794195568;
                const double w2 = 0.02667361780442;

                point(0) = a2;  point(1) = a2;
                weightsAndPoints[4] = std::make_pair( w2, point );

                point(0) = a2; point(1) = b2;
                weightsAndPoints[5] = std::make_pair( w2, point );

                point(0) = b2; point(1) = a2;
                weightsAndPoints[6] = std::make_pair( w2, point );

                const double a3 = 0.048690315425316;
                const double b3 = 0.312865496004875;
                const double c3 = 0.638444188569809;
                const double w3 = 0.038556880445128;

                point(0) = a3; point(1) = b3;
                weightsAndPoints[7] = std::make_pair( w3, point );

                point(0) = b3; point(1) = a3;
                weightsAndPoints[8] = std::make_pair( w3, point );

                point(0) = a3; point(1) = c3;
                weightsAndPoints[9] = std::make_pair( w3, point );

                point(0) = b3; point(1) = c3;
                weightsAndPoints[10] = std::make_pair( w3, point );

                point(0) = c3; point(1) = a3;
                weightsAndPoints[11] = std::make_pair( w3, point );

                point(0) = c3; point(1) = b3;
                weightsAndPoints[12] = std::make_pair( w3, point );
            }
        };

        //! 9th order rule integrating 9-th order polynomials exactly
        template<>
        struct SimplexQuadratureConstructor<2,19>
        {
            void operator()( std::array<std::pair<double,eigenX::VectorSd<2> >, 19> &
                             weightsAndPoints )
            {
                eigenX::VectorSd<2> point;

                const double xi[19] = {
                    .3333333333333329 , .0206349616025250 , .1258208170141270 , .6235929287619349 ,
                    .9105409732110950 , .0368384120547360 , .4896825191987380 , .4896825191987380 ,
                    .4370895914929371 , .4370895914929371 , .1882035356190330 , .1882035356190330 ,
                    .0447295133944530 , .0447295133944530 , .0368384120547360 , .2219629891607660 ,
                    .7411985987844981 , .2219629891607660 , .7411985987844981  };

                const double eta[19] = {
                    .3333333333333329 , .4896825191987380 , .4370895914929371 , .1882035356190330 ,
                    .0447295133944530 , .2219629891607660 , .0206349616025250 , .4896825191987380 ,
                    .1258208170141270 , .4370895914929371 , .6235929287619350 , .1882035356190330 ,
                    .9105409732110950 , .0447295133944530 , .7411985987844981 , .7411985987844981 ,
                    .0368384120547360 , .0368384120547360 , .2219629891607660  };

                const double omega[19] = {
                    .0971357962827990 , .0313347002271390 , .0778275410047740 , .0796477389272100 ,
                    .0255776756586980 , .0432835393772890 , .0313347002271390 , .0313347002271390 ,
                    .0778275410047740 , .0778275410047740 , .0796477389272100 , .0796477389272100 ,
                    .0255776756586980 , .0255776756586980 , .0432835393772890 , .0432835393772890 ,
                    .0432835393772890 , .0432835393772890 , .0432835393772890  };

                for ( unsigned i = 0; i < 19; i ++ ) {
                    point(0) = xi[i]; point(1) = eta[i];
                    weightsAndPoints[i] = std::make_pair( 0.5 * omega[i], point );
                }

            }
        };

        //@}

        //----------------------------------------------------------------------
        //! @name Quadrature points and weights for tetrahedron
        //! \f$(\xi,\eta,\zeta)\in[0,1]\times[0,1-\xi]\times[1-\xi-\eta]\f$
        //@{

        //! 1st order rule integrating linear polynomials exactly
        template<>
        struct SimplexQuadratureConstructor<3,1>
        {
            void operator()( std::array<std::pair<double,eigenX::VectorSd<3> >, 1> &
                             weightsAndPoints )
            {
                eigenX::VectorSd<3> point;
                // Midpoint and volume
                point( 0 ) = point( 1 ) = point( 2 ) = 0.25;
                weightsAndPoints[ 0 ] = std::make_pair( 0.1666666666666666, point );
            }
        };

        //! 2nd order rule integrating quadratic (\f$\xi^a\eta^b\zeta^c\f$ with \f$a+b+c\leq2\f$) polynomials exactly
        template<>
        struct SimplexQuadratureConstructor<3,4>
        {
            void operator()( std::array<std::pair<double,eigenX::VectorSd<3> >, 4> &
                             weightsAndPoints )
            {
                const double a = 0.13819660112501051518; // (5-sqrt(5))/20
                const double b = 0.58541019662496845446; // 1 - 3*a
                const double w = 0.04166666666666666667; // (1/4)*(1/6)
                eigenX::VectorSd<3> p;
                p(0) = a; p(1) = a; p(2) = a; weightsAndPoints[ 0 ] = std::make_pair( w, p );
                p(0) = b; p(1) = a; p(2) = a; weightsAndPoints[ 1 ] = std::make_pair( w, p );
                p(0) = a; p(1) = b; p(2) = a; weightsAndPoints[ 2 ] = std::make_pair( w, p );
                p(0) = a; p(1) = a; p(2) = b; weightsAndPoints[ 3 ] = std::make_pair( w, p );
            }
        };

        //! 3rd order rule integrating cubic (\f$\xi^a\eta^b\zeta^c\f$ with \f$a+b+c\leq3\f$) polynomials exactly
        template<>
        struct SimplexQuadratureConstructor<3,5>
        {
            void operator()( std::array<std::pair<double,eigenX::VectorSd<3> >, 5> &
                             weightsAndPoints )
            {
                eigenX::VectorSd<3> p;
                const double x = 0.1666666666666667; // 1/6
                // weights: -4/30, 3/40
                p(0) = 0.25; p(1) = 0.25; p(2) = 0.25; weightsAndPoints[0] = std::make_pair( -0.13333333333333333333, p );
                p(0) = x;    p(1) = x;    p(2) = x;    weightsAndPoints[1] = std::make_pair( 0.075, p );
                p(0) = 0.5;  p(1) = x;    p(2) = x;    weightsAndPoints[2] = std::make_pair( 0.075, p );
                p(0) = x;    p(1) = 0.5;  p(2) = x;    weightsAndPoints[3] = std::make_pair( 0.075, p );
                p(0) = x;    p(1) = x;    p(2) = 0.5;  weightsAndPoints[4] = std::make_pair( 0.075, p );
            }
        };

        //! 4th order rule integrating quartic (\f$\xi^a\eta^b\zeta^c\f$ with \f$a+b+c\leq4\f$) polynomials exactly
        //!
        //! See http://www.cs.rpi.edu/~flaherje/pdf/fea6.pdf
        template<>
        struct SimplexQuadratureConstructor<3,11>
        {
            void operator()( std::array<std::pair<double,eigenX::VectorSd<3> >, 11> &
                             weightsAndPoints )
            {
                eigenX::VectorSd<3> p;
                p(0) = 0.25; p(1) = 0.25; p(2) = 0.25; weightsAndPoints[0] = std::make_pair( -0.01315555555555555556, p );

                double x = 0.071428571428571, y = 0.785714285714286, w = 0.0076222222222222;
                p(0) = x;    p(1) = x;    p(2) = x;    weightsAndPoints[1] = std::make_pair( w, p );
                p(0) = y;    p(1) = x;    p(2) = x;    weightsAndPoints[2] = std::make_pair( w, p );
                p(0) = x;    p(1) = y;    p(2) = x;    weightsAndPoints[3] = std::make_pair( w, p );
                p(0) = x;    p(1) = x;    p(2) = y;    weightsAndPoints[4] = std::make_pair( w, p );

                x = 0.100596423833201, y = 0.399403576166799, w = 0.024888888888889;
                p(0) = y;    p(1) = y;    p(2) = x;    weightsAndPoints[5] = std::make_pair( w, p );
                p(0) = y;    p(1) = x;    p(2) = x;    weightsAndPoints[6] = std::make_pair( w, p );
                p(0) = x;    p(1) = y;    p(2) = x;    weightsAndPoints[7] = std::make_pair( w, p );
                p(0) = x;    p(1) = x;    p(2) = y;    weightsAndPoints[8] = std::make_pair( w, p );
                p(0) = y;    p(1) = x;    p(2) = y;    weightsAndPoints[9] = std::make_pair( w, p );
                p(0) = x;    p(1) = y;    p(2) = y;    weightsAndPoints[10] = std::make_pair( w, p );
            }
        };

        //! 5th order rule integrating quintic (\f$\xi^a\eta^b\zeta^c\f$ with \f$a+b+c\leq5\f$) polynomials exactly
        //!
        //! See http://www.cs.rpi.edu/~flaherje/pdf/fea6.pdf
        template<>
        struct SimplexQuadratureConstructor<3,15>
        {
            void operator()( std::array<std::pair<double,eigenX::VectorSd<3> >, 15> &
                             weightsAndPoints )
            {
                eigenX::VectorSd<3> p;
                p(0) = 0.25; p(1) = 0.25; p(2) = 0.25; weightsAndPoints[0]  = std::make_pair( 0.030283678097089, p );

                double x = 0.333333333333333, y = 0.0, w = 0.006026785714286;
                p(0) = x;    p(1) = x;    p(2) = x;    weightsAndPoints[1]  = std::make_pair( w, p );
                p(0) = y;    p(1) = x;    p(2) = x;    weightsAndPoints[2]  = std::make_pair( w, p );
                p(0) = x;    p(1) = y;    p(2) = x;    weightsAndPoints[3]  = std::make_pair( w, p );
                p(0) = x;    p(1) = x;    p(2) = y;    weightsAndPoints[4]  = std::make_pair( w, p );

                x = 0.090909090909091, y = 0.727272727272727, w = 0.011645249086029;
                p(0) = x;    p(1) = x;    p(2) = x;    weightsAndPoints[5]  = std::make_pair( w, p );
                p(0) = y;    p(1) = x;    p(2) = x;    weightsAndPoints[6]  = std::make_pair( w, p );
                p(0) = x;    p(1) = y;    p(2) = x;    weightsAndPoints[7]  = std::make_pair( w, p );
                p(0) = x;    p(1) = x;    p(2) = y;    weightsAndPoints[8]  = std::make_pair( w, p );

                x = 0.066550153573664, y = 0.433449846426336, w = 0.010949141561386;
                p(0) = y;    p(1) = y;    p(2) = x;    weightsAndPoints[9]  = std::make_pair( w, p );
                p(0) = y;    p(1) = x;    p(2) = x;    weightsAndPoints[10] = std::make_pair( w, p );
                p(0) = x;    p(1) = y;    p(2) = x;    weightsAndPoints[11] = std::make_pair( w, p );
                p(0) = x;    p(1) = x;    p(2) = y;    weightsAndPoints[12] = std::make_pair( w, p );
                p(0) = y;    p(1) = x;    p(2) = y;    weightsAndPoints[13] = std::make_pair( w, p );
                p(0) = x;    p(1) = y;    p(2) = y;    weightsAndPoints[14] = std::make_pair( w, p );
            }
        };

        //@}

        //! \endcond

    } // end namespace detail_
} // end namespace corlib

//------------------------------------------------------------------------------
namespace corlib {
    //! One-dimensional quadrature DEGREE = 2 * N - 1;
    template<unsigned DEGREE> 
    struct SimplexQuadratureOrder<1,DEGREE>
    {
        static const unsigned numPoints = static_cast<unsigned>( (DEGREE + 2 )/2 );
    };


    //--------------------------------------------------------------------------
    // Triangle rules
    //! \cond SKIPDOX
    template<> 
    struct SimplexQuadratureOrder<2,1> { static const unsigned numPoints = 1; };

    template<> 
    struct SimplexQuadratureOrder<2,2> { static const unsigned numPoints = 3; };

    template<> 
    struct SimplexQuadratureOrder<2,3> { static const unsigned numPoints = 4; };

    template<> 
    struct SimplexQuadratureOrder<2,4> { static const unsigned numPoints = 6; };

    template<> 
    struct SimplexQuadratureOrder<2,5> { static const unsigned numPoints = 7; };

    //--------------------------------------------------------------------------
    // Tetrehedron rules
    template<> 
    struct SimplexQuadratureOrder<3,1> { static const unsigned numPoints = 1; };

    template<> 
    struct SimplexQuadratureOrder<3,2> { static const unsigned numPoints = 4; };

    template<> 
    struct SimplexQuadratureOrder<3,3> { static const unsigned numPoints = 5; };

    template<> 
    struct SimplexQuadratureOrder<3,4> { static const unsigned numPoints = 11; };

    template<> 
    struct SimplexQuadratureOrder<3,5> { static const unsigned numPoints = 15; };
    //! \endcond

}


#endif
