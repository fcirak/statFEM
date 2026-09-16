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

#ifndef subdiv_surf_boxspline_h
#define subdiv_surf_boxspline_h

//------------------------------------------------------------------------------
#include <Eigen/Core>
#include <corlib/eigenX.hpp>

//------------------------------------------------------------------------------
namespace subdiv{
    namespace surf{

        class BoxSpline;

        namespace eigenX = corlib::eigenX;

    }
}

//------------------------------------------------------------------------------
/// Quartic box spline on triangular grid with 5 generating
/// direction vectors \f$(1,1),(1,0),(0,1),(0,1),(1,1)\f$.
/// The class computes shape functions of the 12 box-splines which have
/// non-zero support over an triangular element.
///
/// Ordering of functions around element with vertices (4,7,8)
///\verbatim
///    9--12        
///    |`  |`
///    | ` | `      w=xi_2 ^
///    |  `|  `            |
///    5---8--11           1 8
///    |`  |`  |`          | |`
///    | ` |*` | `         | |*`
///    |  `|**`|  `        | |**`
///    2---4---7--10       0 4---7
///     `  |`  |`  |
///      ` | ` | ` |         0---1--> xi_1=v
///       `|  `|  `|
///        1---3---6
///\endverbatim
class subdiv::surf::BoxSpline
{
public:
    static const unsigned numSplines   = 12;
    static const unsigned localDim     = 2;
    
    typedef eigenX::VectorSd<localDim>                VecLDim;
    typedef eigenX::VectorSd<numSplines>              VecNspl;
    typedef eigenX::MatrixSd<localDim,numSplines>     MatLDimNspl;
    typedef Eigen::Matrix<VecNspl, 2, 2>              MatVecNspl2x2;
public:
    /** Evaluation of BoxSpline FE shape functions of order DEGREE
     *  \param[in]  xi  Local coordinates from [0,1]
     *  \param[out] phi Result \f$\phi_i(\xi)\f$      */
    void evaluate( const VecLDim & xi, VecNspl & phi ) const;

    /** Evaluate the first derivative of the Box Spline FE functions
     *  \param[in]  xi   Local coordinates from [0,1]
     *  \param[out] dphi Result \f$ \phi^\prime_i(\xi) \f$     */
    void evaluateGradient( const VecLDim & xi, MatLDimNspl & dphi ) const;

    /** Evaluate the second derivative of the Box Spline FE functions
     *  \param[in]  xi    Local coordinates from [0,1]
     *  \param[out] ddphi Result \f$ \phi^{\prime\prime}_i(\xi) \f$ */
    void evaluateHessian( const VecLDim & xi, MatVecNspl2x2 & ddphi ) const;
};

//------------------------------------------------------------------------------
void subdiv::surf::BoxSpline::evaluate( const VecLDim & xi, VecNspl & phi ) const
{ 
    const double v = xi( 0 );
    const double w = xi( 1 );
    const double u = 1.0 - v - w;

    const double v2 = v * v;
    const double w2 = w * w;
    const double u2 = u * u;

    const double v3 = v2 * v;
    const double w3 = w2 * w;
    const double u3 = u2 * u;

    const double v4 = v3 * v;
    const double w4 = w3 * w;
    const double u4 = u3 * u;

 
    phi(0)=(u4 + 2.*u3*v)/12.;
    
    phi(1)=(u4 + 2.*u3*w)/12.;
 
    phi(2)=(u4 + 2.*u3*w + 6.*u3*v + 6.*u2*v*w + 
            12.*u2*v2 + 6.*u*v2*w + 6.*u*v3 + 2.*v3*w + 
            v4)/12.;
    
    phi(3)=(6.*u4 + 24.*u3*w + 24.*u2*w2 + 8.*u*w3 + 
            w4 + 24.*u3*v + 60.*u2*v*w + 36.*u*v*w2 + 
            6.*v*w3 + 24.*u2*v2 + 36.*u*v2*w + 12.*v2*w2 + 
            8.*u*v3 + 6.*v3*w + v4)/12.;
    
    phi(4)=(u4 + 6.*u3*w + 12.*u2*w2 + 6.*u*w3 + 
            w4 + 2.*u3*v + 6.*u2*v*w + 6.*u*v*w2 + 
            2.*v*w3)/12.;
    
    phi(5)=(2.*u*v3 + v4)/12.;
    
    phi(6)=(u4 + 6.*u3*w + 12.*u2*w2 + 
            6.*u*w3 + w4 + 8.*u3*v + 36.*u2*v*w + 
            36.*u*v*w2 + 8.*v*w3 + 24.*u2*v2 + 60.*u*v2*w + 
            24.*v2*w2 + 24.*u*v3 + 24.*v3*w + 
            6.*v4)/12.;
    
    phi(7)=(u4 + 8.*u3*w + 24.*u2*w2 + 24.*u*w3 + 
            6.*w4 + 6.*u3*v + 36.*u2*v*w + 60.*u*v*w2 + 
            24.*v*w3 + 12.*u2*v2 + 36.*u*v2*w + 
            24.*v2*w2 + 6.*u*v3 + 8.*v3*w + v4)/12.;
    
    phi(8)=(2.*u*w3 + w4)/12.;
  
    phi(9)=(2.*v3*w + v4)/12.;
    
    
    phi(10)=(2.*u*w3 + w4 + 6.*u*v*w2 + 6.*v*w3 + 
             6.*u*v2*w + 12.*v2*w2 + 2.*u*v3 + 
             6.*v3*w + v4)/12.;
    
    phi(11)=(w4 + 2.*v*w3)/12.;

    return;
}

//------------------------------------------------------------------------------
void subdiv::surf::BoxSpline::evaluateGradient( 
    const VecLDim & xi,
    MatLDimNspl & dphi
    ) const
{ 
    const double v = xi( 0 );
    const double w = xi( 1 );
    const double u = 1.0 - v - w;

    const double v2 = v * v;
    const double w2 = w * w;
    const double u2 = u * u;

    const double v3 = v2 * v;
    const double w3 = w2 * w;
    const double u3 = u2 * u;


    dphi(0,0) = (-6.0*v*u2 - 2.0*u3)/12.0;
    dphi(1,0) = (-6.0*v*u2 - 4.0*u3)/12.0;
 
    dphi(0,1) = (-4.0*u3-6.0*u2*w)/12.0;
    dphi(1,1) = (-2.0*u3-6.0*u2*w)/12.0;
    
    dphi(0,2) = (-2.0*v3-6.0*v2*u
                 + 6.0*v*u2+2.0*u3)/12.0;
    dphi(1,2) = (-4.0*v3-18.0*v2*u
                 - 12.0*v*u2-2.0*u3
                 - 6.0*v2*w-12.0*v*u*w
                 - 6.0*u2*w)/12.0;
    
    dphi(0,3) = (-4.0*v3-24.0*v2*u
                 - 24.0*v*u2-18.0*v2*w 
                 - 48.0*v*u*w-12.0*u2*w
                 - 12.0*v*w2 - 12.0*u*w2
                 - 2.0*w3)/12.0;
    
    dphi(1,3) = (-2.0*v3-12.0*v2*u
                 - 12.0*v*u2-12.0*v2*w
                 - 48.0*v*u*w-24.0*u2*w
                 - 18.0*v*w2-24.0*u*w2
                 - 4.0*w3)/12.0;
    
    dphi(0,4) = (-6.0*v*u2-2.0*u3
                 - 12.0*v*u*w-12.0*u2*w
                 - 6.0*v*w2-18.0*u*w2
                 - 4.0*w3)/12.0;
    
    dphi(1,4) = (2.0*u3+6.0*u2*w
                 - 6.0*u*w2-2.0*w3)/12.0;
    
    dphi(0,5) = (2.0*v3+6.0*v2*u)/12.0;
    dphi(1,5) = -v3/6.0;
    
    dphi(0,6) = (24.0*v2*u+24.0*v*u2
                 + 4.0*u3+12.0*v2*w
                 + 48.0*v*u*w+18.0*u2*w
                 + 12.0*v*w2+12.0*u*w2
                 + 2.0*w3)/12.0;
    
    dphi(1,6) = (12.0*v2*u+12.0*v*u2
                 + 2.0*u3-12.0*v2*w
                 + 6.0*u2*w-12.0*v*w2
                 - 6.0*u*w2-2.0*w3)/12.0;
    
    dphi(0,7) = (-2.0*v3-6.0*v2*u
                 + 6.0*v*u2+2.0*u3
                 - 12.0*v2*w+12.0*u2*w
                 - 12.0*v*w2+12.0*(1.0-v-w)*w2)/12.0;
    
    dphi(1,7) = (2.0*v3+12.0*v2*u
                 + 18.0*v*u2+4.0*u3
                 + 12.0*v2*w+48.0*v*u*w
                 + 24.0*u2*w+12.0*v*w2
                 + 24.0*u*w2)/12.0;
    
    dphi(0,8) = -w3/6.0;
    dphi(1,8) = (6.0*u*w2+2.0*w3)/12.0;
    
    dphi(0,9) = (4.0*v3+6.0*v2*w)/12.0;
    dphi(1,9) = v3/6.0;
    
    dphi(0,10)= (2.0*v3+6.0*v2*u
                 + 12.0*v2*w+12.0*v*u*w
                 + 18.0*v*w2+6.0*u*w2
                 + 4.0*w3)/12.0;
    
    dphi(1,10)= (4.0*v3+6.0*v2*u
                 + 18.0*v2*w+12.0*v*u*w
                 + 12.0*v*w2+6.0*u*w2
                 + 2.0*w3)/12.0;
    
    dphi(0,11) = w3/6.0;
    dphi(1,11) = (6.0*v*w2+4.0*w3)/12.0;
   
    return;
}

//------------------------------------------------------------------------------
void subdiv::surf::BoxSpline::evaluateHessian( 
    const VecLDim & xi,
    MatVecNspl2x2 & ddphi
    ) const
{
    const double v = xi( 0 );
    const double w = xi( 1 );
    const double u = 1.0 - v - w;

    const double v2 = v * v;
    const double w2 = w * w;
    const double u2 = u * u;


    ddphi(0,0)(0) = v*u;
    ddphi(1,1)(0) = v*u+u2;
    ddphi(0,1)(0) = (12.0*v*u+6.0*u2)/12.0;
 
    ddphi(0,0)(1) = u2+u*w;
    ddphi(1,1)(1) = u*w;
    ddphi(0,1)(1) = (6.0*u2+12.0*u*w)/12.0;
 
    ddphi(0,0)(2) = -2.0*v*u;
    ddphi(1,1)(2) = v2+v*u+v*w+u*w;
    ddphi(0,1)(2) = (6.0*v2-12.0*v*u
                     -6.0*u2)/12.0;
 
    ddphi(0,0)(3) = v2-2.0*u2+ v*w-2.0*u*w;
    ddphi(1,1)(3) = -2.0*v*u-2.0*u2+ v*w+w2;
    ddphi(0,1)(3) = (6.0*v2-12.0*u2+ 24.0*v*w+6.0*w2)/12.0;
 
    ddphi(0,0)(4) = v*u+v*w+u*w+ w2;
    ddphi(1,1)(4) = -2.0*u*w;
    ddphi(0,1)(4) = (-6.0*u2-12.0*u*w + 6.0*w2)/12.0;
 
    ddphi(0,0)(5) = v*u;
    ddphi(1,1)(5) = 0.0;
    ddphi(0,1)(5) = -v2/2.0;
 
    ddphi(0,0)(6) = (-24.0*v2+12.0*u2-24.0*v*w + 12.0*u*w)/12.0;
    ddphi(1,1)(6) = (-24.0*v2-24.0*v*u-24.0*v*w - 24.0*u*w)/12.0;
    ddphi(0,1)(6) = (-12.0*v2+6.0*u2-24.0*v*w - 12.0*u*w-6.0*w2)/12.0;
 
    ddphi(0,0)(7) = -2.0*v*u-2.0*v*w-2.0*u*w - 2.0*w2;
    ddphi(1,1)(7) = v*u+u2-2.0*v*w - 2.0*w2;
    ddphi(0,1)(7) = (-6.0*v2-12.0*v*u+6.0*u2- 24.0*v*w-12.0*w2)/12.0;
 
    ddphi(0,0)(8) = 0.0;
    ddphi(1,1)(8) = u*w;
    ddphi(0,1)(8) = -w2/2.0; 
 
    ddphi(0,0)(9) = (12.0*v2+12.0*v*w)/12.0;
    ddphi(1,1)(9) = 0.0;
    ddphi(0,1)(9) = v2/2.0;
 
    ddphi(0,0)(10)= (12.0*v*u+12.0*v*w+12.0*u*w + 12.0*w2)/12.0;
    ddphi(1,1)(10)= v2+v*u+v*w+u*w;
    ddphi(0,1)(10)= (6.0*v2+12.0*v*u+24.0*v*w+ 12.0*u*w+6.0*w2)/12.0;
 
    ddphi(0,0)(11)= 0.0;
    ddphi(1,1)(11)= v*w+w2;
    ddphi(0,1)(11)= w2/2.0;

    for ( unsigned i = 0; i < numSplines; i++ )
        ddphi(1,0)(i)= ddphi(0,1)(i);

    return;
}

#endif
