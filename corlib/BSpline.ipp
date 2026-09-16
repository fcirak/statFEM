// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//
//                 Computational Structural Mechanics Lab
//                         University of Cambridge
//
// This software is copyrighted and all rights are retained by the University of
// Cambridge. All use, disclosure, and/or reproduction of any part not expressly
// authorised by the University of Cambridge is prohibited. (C) 2010
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

//! @file BSpline.ipp

#include <corlib/verify.hpp>

namespace corlib{

    //--------------------------------------------------------------------------
    //! \brief Characteristic function, i.e., constant spline
    template<> 
    void BSpline<0>::evaluate( const double & xi, VecNspl & phi ) const
    {
        phi(0) = 1.;
        return;
    }

    // template<>
    // void BSpline<0>::evaluateGradient( const double & xi, VecNspl & dphi ) 
    // is not implemented --> trigger compile-time error

    //--------------------------------------------------------------------------
    /** \brief Linear splines on (0,1)
     *  \details <pre>
     *
     *   1 **------------+------------+-------------+------------+------------##
     *     + ****        +            + 1-x ******  +            +        #### +
     *     |    ****                      x ######                     ####    |
     *     |        ****                                           ####        |
     * 0.8 ++          ****                                     ####          ++
     *     |       phi_1   ****                             ####  phi_2        |
     *     |                  ****                       ####                  |
     *     |                      ***                 ###                      |
     * 0.6 ++                        ****         ####                        ++
     *     |                             ***   ###                             |
     *     |                                ***                                |
     *     |                             ###   ***                             |
     * 0.4 ++                        ####         ****                        ++
     *     |                      ###                 ***                      |
     *     |                  ####                       ****                  |
     *     |               ####                             ****               |
     * 0.2 ++          ####                                     ****          ++
     *     |        ####                                           ****        |
     *     |    ####                                                   ****    |
     *     + ####        +            +             +            +        **** +
     *   0 ##------------+------------+-------------+------------+------------**
     *     0            0.2          0.4           0.6          0.8            1
     *  </pre>
     */
    template<> 
    void BSpline<1>::evaluate( const double & xi, VecNspl & phi ) const
    {
        phi(0) = 1. - xi;
        phi(1) = xi;
        return;
    }

    //! Gradient of linear splines
    template<>
    void BSpline<1>::evaluateGradient( const double & xi, VecNspl & dphi ) const
    {
        dphi(0) = -1.;
        dphi(1) =  1.;
        return;
    }
        
    //--------------------------------------------------------------------------
    /** \brief Quadratic splines on (0,1)
     *  \details <pre>
     *
     * 0.8 ++------------+------------+-------------+------------+------------++
     *     +             +            ###############            +             +
     *     |                  #########   phi_2     #########                  |
     * 0.7 ++            ######                             ######            ++
     *     |         ####                                         ####         |
     * 0.6 ++    ####                                                 ####    ++
     *     |  ###                                                         ###  |
     * 0.5 ####                                                             ###*
     *     |$$$                                                             ***|
     *     |  $$$                        1./2.*x*x ******                 ***  |
     * 0.4 ++    $$$  phi_1           -x*x+x+1./2. ######      phi_3   ***    ++
     *     |        $$$$         1./2.*(1-x)*(1-x) $$$$$$          ****        |
     * 0.3 ++          $$$$                                     ****          ++
     *     |               $$$$                             ****               |
     *     |                  $$$$$                     *****                  |
     * 0.2 ++                     $$$$$             *****                     ++
     *     |                          $$$$$$   ******                          |
     * 0.1 ++                             ******$                             ++
     *     |                       *******       $$$$$$$                       |
     *     +             **********   +             +   $$$$$$$$$$             +
     *   0 **************+------------+-------------+------------+$$$$$$$$$$$$$$
     *     0            0.2          0.4           0.6          0.8            1
     * </pre>
     */
    template<>
    void BSpline<2>::evaluate( const double & xi, VecNspl & phi ) const
    {
        phi(0) = 0.5 * (1.-xi) * (1.-xi);
        phi(1) = 0.5 + xi - xi * xi;
        phi(2) = 0.5 * xi * xi;
        return;
    }

    //! Gradient of quadratic splines
    template<>
    void BSpline<2>::evaluateGradient( const double & xi, VecNspl & dphi ) const
    {
        dphi(0) = - (1.-xi);
        dphi(1) = 1. - 2. * xi;
        dphi(2) = xi;
        return;
    }

    //! Hessian of quadratic splines
    template<>
    void BSpline<2>::evaluateHessian( const double & xi, VecNspl & ddphi ) const
    {
        ddphi(0) = 1.;
        ddphi(1) = -2.;
        ddphi(2) = 1.;
        return;
    }

    //--------------------------------------------------------------------------
    /** \brief Cubic splines on (0,1)
     *  \details <pre>
     *
     *    0.7 ++------------+------------+-------------+------------+------------++
     *        ###########   +            +             +            +   $$$$$$$$$$$
     *        |         #######                                   $$$$$$$         |
     *    0.6 ++        phi_2  ######                       $$$$$$ phi_3         ++
     *        |                      ####               $$$$                      |
     *        |                          ####       $$$$                          |
     *    0.5 ++                             ####$$$                             ++
     *        |                             $$$$ ####                             |
     *        |                         $$$$$       #####                         |
     *    0.4 ++                     $$$                 ###                     ++
     *        |                  $$$$                       ####                  |
     *    0.3 ++              $$$$                             ####              ++
     *        |           $$$$                                     ####           |
     *        |       $$$$                                             ####       |
     *    0.2 ++ $$$$$              (1.-x)*(1.-x)*(1.-x)/6. ******         ##### ++
     *        **$$                  (3.*x*x*x-6.*x*x+4.)/6. ######             ####
     *        |******         (-3.*x*x*x+3.*x*x+3.*x+1.)/6. $$$$$$          %%%%%%|
     *    0.1 ++     *******                       x*x*x/6. %%%%%%   %%%%%%%     ++
     *        |      phi_1  ********                         %%%%%%%%  phi_4      |
     *        +             +       ***************%%%%%%%%%%       +             +
     *      0 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%---*********************************
     *        0            0.2          0.4           0.6          0.8            1
     *  </pre>
     */
    template<>
    void BSpline<3>::evaluate( const double & xi, VecNspl & phi ) const
    {
        phi(0) = (1. - xi) * (1. - xi) * (1. - xi) / 6.;
        phi(1) = (3. * xi * xi * xi - 6. * xi * xi + 4.)/ 6.;
        phi(2) = (-3. * xi * xi * xi + 3. * xi * xi + 3. * xi + 1.)/ 6.;
        phi(3) = xi * xi * xi / 6.;
        return;
    }

    //!  Gradient of cubic splines
    template<>
    void BSpline<3>::evaluateGradient( const double & xi, VecNspl & dphi ) const
    {
        dphi(0) = - (1. - xi) * (1. - xi) / 2.;
        dphi(1) = (9. * xi * xi - 12. * xi) / 6.;
        dphi(2) = (-9. * xi * xi + 6. * xi + 3.) / 6.;
        dphi(3) = xi * xi / 2.;
        return;
    }

    //!  Hessian of cubic splines
    template<>
    void BSpline<3>::evaluateHessian( const double & xi, VecNspl & ddphi ) const
    {
        ddphi(0) = 1. - xi;
        ddphi(1) = 3. * xi - 2.;
        ddphi(2) = 1. - 3. * xi;
        ddphi(3) = xi;
        return;
    }

    //--------------------------------------------------------------------------
    /** Evaluate higher order splines by means of recursion
     *  \f[
     *  \phi^n_k = \frac{1}{n} [(\xi+n-k) \phi^{n-1}_{k-1} + (k+1-\xi) \phi^{n-1}_k]
     *  \f]
     */
    template<unsigned DEGREE>
    void BSpline<DEGREE>::evaluate( const double & xi, VecNspl & phi ) const
    {
        BSpline<degree-1> lowerOrderSpline;
        typename BSpline<degree-1>::VecNspl lowerPhi;
        lowerOrderSpline.evaluate( xi, lowerPhi );
        const double factor = 1./static_cast<double>(degree);

        phi[0] = factor * (1. - xi) * lowerPhi[0];
        for ( unsigned k = 1; k < numSplines-1; k ++ ) {
            phi[k] = factor * ( (xi+degree-k) * lowerPhi[k-1] + (k+1.-xi) * lowerPhi[k] );
        }
        phi[degree] = factor * (xi * lowerPhi[degree-1] );

        return;
    }
        
    //--------------------------------------------------------------------------
    /** Evaluate higher order splines by means of recursion
     *  \f[
     *  \frac{\partial \phi^n_k}{\partial \xi} = \phi^{n-1}_{k-1} - \phi^{n-1}_k
     *  \f]
     */
    template<unsigned DEGREE>
    void BSpline<DEGREE>::evaluateGradient( const double & xi, VecNspl & dphi ) const
    {
        static_assert( (degree > 0), "Cannot compute gradient for lower order" );
        BSpline<degree-1> lowerOrderSpline;
        typename BSpline<degree-1>::VecNspl lowerPhi;
        lowerOrderSpline.evaluate( xi, lowerPhi );

        dphi[0] = -lowerPhi[0];
        for ( unsigned k = 1; k < numSplines-1; k ++ ) {
            dphi[k] = lowerPhi[k-1] - lowerPhi[k];
        }
        dphi[degree] = lowerPhi[degree-1];

        return;
    }

    //--------------------------------------------------------------------------
    /** Evaluate higher order splines by means of recursion
     *  \f[
     *       \frac{\partial^2 \phi^n_k}{\partial \xi^2} = 
     *       \phi^{n-1}_{k-2} - 2 \phi^{k-1}_{k-1} + \phi^{k-1}_k
     *  \f]
     */
    template<unsigned DEGREE>
    void BSpline<DEGREE>::evaluateHessian( const double & xi, VecNspl & ddphi ) const
    {
        static_assert( (degree > 1), "Cannot compute Hessian for lower order" );
        BSpline<degree-2> lowerOrderSpline;
        typename BSpline<degree-2>::VecNspl lowerPhi;
        lowerOrderSpline.evaluate( xi, lowerPhi );

        ddphi[0] = lowerPhi[0];
        ddphi[1] = - 2. * lowerPhi[0] + lowerPhi[1];
        for ( unsigned k = 2; k < numSplines-2; k ++ ) {
            ddphi[k] = lowerPhi[k-2] - 2. * lowerPhi[k-1] + lowerPhi[k];
        }
        ddphi[degree-1] = lowerPhi[degree-3] - 2. * lowerPhi[degree-2];
        ddphi[degree]   = lowerPhi[degree-2];

        return;
    }

}


