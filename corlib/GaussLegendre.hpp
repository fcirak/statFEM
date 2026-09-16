// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//
//                 Computational Structural Mechanics Lab
//                         University of Cambridge
//
// This software is copyrighted and all rights are retained by the University of
// Cambridge. All use, disclosure, and/or reproduction of any part not expressly
// authorised by the University of Cambridge is prohibited. (C) 2010
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

//! @file   GaussLegendre.hpp

#ifndef corlib_gausslegendre_h
#define corlib_gausslegendre_h
//------------------------------------------------------------------------------
#include <utility>
#include <array>

//------------------------------------------------------------------------------
namespace corlib{
    template<unsigned NPTS> class GaussLegendre;
}

//------------------------------------------------------------------------------
/** \brief Gauss-Legendre quadrature rules with NPTS points on (0,1) 
 *  \tparam NPTS  Number of points
 */
template<unsigned NPTS> 
class corlib::GaussLegendre
{
public:
    static const unsigned numPoints = NPTS;

    //! Specialized constructor for different NPTS
    GaussLegendre();

    //! Random access to any entry in the array of (weight,point) pairs
    std::pair<double,double> operator[]( const unsigned index ) const
    {
        return weightsAndPoints_[index];
    }

private:
    //! Pair of weights and points representing the NPTS-point quadrature rule
    std::array< std::pair<double,double>, numPoints > weightsAndPoints_;
};


//------------------------------------------------------------------------------
/** The implementation of Gauss-Legendre quadrature rules for
 *        the reference interval (0,1).
 */
//------------------------------------------------------------------------------
//! \cond SKIPDOX
namespace corlib{
        
    //! 1-point rule
    template<> GaussLegendre<1>::GaussLegendre()
    {
        weightsAndPoints_[0] =  std::make_pair( 1, 0.5);
    }


    //! 2-point rule
    template<> GaussLegendre<2>::GaussLegendre()
    {
        weightsAndPoints_[0] =  std::make_pair( 0.5, 0.788675134594813);
        weightsAndPoints_[1] =  std::make_pair( 0.5, 0.211324865405187);
    }


    //! 3-point rule
    template<> GaussLegendre<3>::GaussLegendre()
    {
        weightsAndPoints_[0] =  std::make_pair( 0.277777777777777, 0.887298334620741);
        weightsAndPoints_[1] =  std::make_pair( 0.444444444444444, 0.5);
        weightsAndPoints_[2] =  std::make_pair( 0.277777777777777, 0.112701665379259);
    }


    //! 4-point rule
    template<> GaussLegendre<4>::GaussLegendre()
    {
        weightsAndPoints_[0] =  std::make_pair( 0.173927422568727, 0.930568155797026);
        weightsAndPoints_[1] =  std::make_pair( 0.326072577431273, 0.669990521792428);
        weightsAndPoints_[2] =  std::make_pair( 0.326072577431273, 0.330009478207572);
        weightsAndPoints_[3] =  std::make_pair( 0.173927422568727, 0.069431844202974);
    }


    //! 5-point rule
    template<> GaussLegendre<5>::GaussLegendre()
    {
        weightsAndPoints_[0] =  std::make_pair( 0.118463442528095, 0.953089922969332);
        weightsAndPoints_[1] =  std::make_pair( 0.239314335249683, 0.769234655052841);
        weightsAndPoints_[2] =  std::make_pair( 0.284444444444444, 0.5);
        weightsAndPoints_[3] =  std::make_pair( 0.239314335249683, 0.230765344947159);
        weightsAndPoints_[4] =  std::make_pair( 0.118463442528095, 0.046910077030668);
    }


    //! 6-point rule
    template<> GaussLegendre<6>::GaussLegendre()
    {
        weightsAndPoints_[0] =  std::make_pair( 0.085662246189585, 0.966234757101576);
        weightsAndPoints_[1] =  std::make_pair( 0.180380786524069, 0.830604693233132);
        weightsAndPoints_[2] =  std::make_pair( 0.233956967286345, 0.619309593041598);
        weightsAndPoints_[3] =  std::make_pair( 0.233956967286345, 0.380690406958402);
        weightsAndPoints_[4] =  std::make_pair( 0.180380786524069, 0.169395306766868);
        weightsAndPoints_[5] =  std::make_pair( 0.085662246189585, 0.033765242898424);
    }


    //! 7-point rule
    template<> GaussLegendre<7>::GaussLegendre()
    {
        weightsAndPoints_[0] =  std::make_pair( 0.064742483084435, 0.974553956171379);
        weightsAndPoints_[1] =  std::make_pair( 0.139852695744638, 0.870765592799697);
        weightsAndPoints_[2] =  std::make_pair( 0.190915025252559, 0.702922575688699);
        weightsAndPoints_[3] =  std::make_pair( 0.208979591836735, 0.5);
        weightsAndPoints_[4] =  std::make_pair( 0.190915025252559, 0.297077424311301);
        weightsAndPoints_[5] =  std::make_pair( 0.139852695744638, 0.129234407200303);
        weightsAndPoints_[6] =  std::make_pair( 0.064742483084435, 0.025446043828621);
    }


    //! 8-point rule
    template<> GaussLegendre<8>::GaussLegendre()
    {
        weightsAndPoints_[0] =  std::make_pair( 0.0506142681451885, 0.980144928248768);
        weightsAndPoints_[1] =  std::make_pair( 0.111190517226687, 0.898333238706813);
        weightsAndPoints_[2] =  std::make_pair( 0.156853322938943, 0.762766204958164);
        weightsAndPoints_[3] =  std::make_pair( 0.181341891689181, 0.591717321247825);
        weightsAndPoints_[4] =  std::make_pair( 0.181341891689181, 0.408282678752175);
        weightsAndPoints_[5] =  std::make_pair( 0.156853322938943, 0.237233795041836);
        weightsAndPoints_[6] =  std::make_pair( 0.111190517226687, 0.101666761293187);
        weightsAndPoints_[7] =  std::make_pair( 0.0506142681451885, 0.019855071751232);
    }


    //! 9-point rule
    template<> GaussLegendre<9>::GaussLegendre()
    {
        weightsAndPoints_[0] =  std::make_pair( 0.0406371941807875, 0.984080119753813);
        weightsAndPoints_[1] =  std::make_pair( 0.0903240803474285, 0.918015553663318);
        weightsAndPoints_[2] =  std::make_pair( 0.130305348201468, 0.806685716350295);
        weightsAndPoints_[3] =  std::make_pair( 0.156173538520001, 0.662126711701905);
        weightsAndPoints_[4] =  std::make_pair( 0.16511967750063, 0.5);
        weightsAndPoints_[5] =  std::make_pair( 0.156173538520001, 0.337873288298095);
        weightsAndPoints_[6] =  std::make_pair( 0.130305348201468, 0.193314283649705);
        weightsAndPoints_[7] =  std::make_pair( 0.0903240803474285, 0.081984446336682);
        weightsAndPoints_[8] =  std::make_pair( 0.0406371941807875, 0.015919880246187);
    }


    //! 10-point rule
    template<> GaussLegendre<10>::GaussLegendre()
    {
        weightsAndPoints_[0] =  std::make_pair( 0.033335672154344, 0.986953264258586);
        weightsAndPoints_[1] =  std::make_pair( 0.0747256745752905, 0.932531683344493);
        weightsAndPoints_[2] =  std::make_pair( 0.109543181257991, 0.839704784149512);
        weightsAndPoints_[3] =  std::make_pair( 0.134633359654998, 0.716697697064623);
        weightsAndPoints_[4] =  std::make_pair( 0.147762112357376, 0.574437169490815);
        weightsAndPoints_[5] =  std::make_pair( 0.147762112357376, 0.425562830509185);
        weightsAndPoints_[6] =  std::make_pair( 0.134633359654998, 0.283302302935377);
        weightsAndPoints_[7] =  std::make_pair( 0.109543181257991, 0.160295215850488);
        weightsAndPoints_[8] =  std::make_pair( 0.0747256745752905, 0.0674683166555075);
        weightsAndPoints_[9] =  std::make_pair( 0.033335672154344, 0.013046735741414);
    }

}
//! \endcond
//------------------------------------------------------------------------------
#endif
