// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//
//                 Computational Structural Mechanics Lab
//                         University of Cambridge
//
// This software is copyrighted and all rights are retained by the University of
// Cambridge. All use, disclosure, and/or reproduction of any part not expressly
// authorised by the University of Cambridge is prohibited. (C) 2010
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

//! @file BackwardEulerAssembler.hpp

#ifndef corlib_backwardeulerassembler_h
#define corlib_backwardeulerassembler_h
//------------------------------------------------------------------------------
//! System includes
#include <functional>
//! Corlib includes
#include <corlib/NewmarkAssembler.hpp>

//------------------------------------------------------------------------------
namespace corlib {

    template<typename SYSSOLVER, typename ELEMENT, 
             typename MASS, typename STIFF>
    class BackwardEulerMatrixAssembler;

    template<typename SYSSOLVER, typename ELEMENT, typename MASS,
             typename STIFF = detail_::DoNothing<ELEMENT> >
    class BackwardEulerInertiaForceAssembler;

    template<typename NODE>
    class BackwardEulerSolutionDistributor;
}


//------------------------------------------------------------------------------
/**  \brief Assembly of system matrix for an Euler Backward time integrator 
 *   \details
 *   Using Euler backward for time integration, the system becomes
 *   \f[
 *       \left[ \frac{1}{\Delta t^2} M + \frac{1}{\Delta t} C + K \right]
 *         (x^{(k+1)}_{n+1} - x^{(k)}) = RHS
 *   \f]
 *  \tparam SYSSOLVER Type of system solver to receive the final system matrix
 *  \tparam ELEMENT   Type of finite element which delivers the matrices
 *  \tparam MASS      A functor which gives the element mass matrix
 *  \tparam STIFF     A functor which gives the element's linearised stiffness
 */
template<typename SYSSOLVER, typename ELEMENT, typename MASS, typename STIFF>
class corlib::BackwardEulerMatrixAssembler
    : public corlib::NewmarkMatrixAssembler<SYSSOLVER, ELEMENT, MASS, STIFF>
{
public:
    //! Constructor which allows Backward Euler
    //!
    //! \param[in]   dt       time step size 
    //! \param[in]   mass     functor for the mass matrix computation
    //! \param[in]   stiff    functor for the stiffness matrix computation
    //! \param[in]   d1,d2    Rayleigh damping: C = d1 * M + d2 * K
    BackwardEulerMatrixAssembler( const double dt, 
                                  MASS  & mass, 
                                  STIFF & stiff, 
                                  const double d1 = 0.,
                                  const double d2 = 0. )
        : corlib::NewmarkMatrixAssembler<SYSSOLVER,ELEMENT, 
                                         MASS,STIFF>( 1./(dt * dt) + d1/dt, 
                                                      1.           + d2/dt,
                                                      mass, stiff )
    { }
};

//------------------------------------------------------------------------------
/** \brief Inertia forces for the BackwardEuler time stepping algorithm
 *
 *  \details The complete RHS in case of a Newmark time integration scheme reads
 *  \f[
 *     RHS = f^{ext}_{n+1} - f^{int}(x_{n+1}^{(k)})
 *      - M \left[ \frac{1}{\Delta t^2} (x^{(k)}_{n+1} - x_n) +
 *                 \frac{-1}{\Delta t} v_n + 0 a_n \right]
 *      - C \left[ \frac{1}{\Delta t} (x^{(k)}_{n+1} - x_n) +
 *                 0 v_n + 0 a_n \right]
 *  \f]
 *  It is assumed the the external and internal forces, \f$f^{ext}\f$ and
 *  \f$ f^{int} \f$ are computed separately. 
 *  \tparam SYSSOLVER  Final storage of the computed result
 *  \tparam ELEMENT Type of element to give the data
 *  \tparam MASS      A functor which gives the element mass matrix
 *  \tparam STIFF     A functor which gives the element's linearised stiffness
 */
template<typename SYSSOLVER, typename ELEMENT, 
         typename MASS,   typename STIFF>
class corlib::BackwardEulerInertiaForceAssembler
    : public corlib::NewmarkInertiaForceAssembler<SYSSOLVER, ELEMENT, MASS, STIFF>
{
public:
    //! Constructor
    //!
    //! \param[in]   dt      time step size
    //! \param[in]   mass    functor for the mass matrix computation
    //! \param[in]   stiff   (Optional) stiffness functor (needed for damping)
    //! \param[in]   d1,d2   Rayleigh damping: C = d1 * M + d2 * K
    BackwardEulerInertiaForceAssembler( const double dt, MASS & mass, 
                                        STIFF stiff = STIFF(), 
                                        const double d1 = 0., 
                                        const double d2 = 0. )
        : corlib::NewmarkInertiaForceAssembler<SYSSOLVER,ELEMENT,MASS,
                                               STIFF>( 1./dt,         0.0, 0.0,
                                                       1./(dt*dt), -1./dt, 0.0,               
                                                       mass, stiff,
                                                       d1, d2 )
    { }
};

//------------------------------------------------------------------------------
/** \brief Update the nodal velocities and accelerations a la Euler backward
 *  \details
 *  The update rules of the Euler backward time integration scheme are
 *  \f[
 *      v_{n+1} = \frac{1}{\Delta t} (x_{n+1} - x_n) + 0 v_n + 0 a_n 
 *  \f]
 *  and
 *  \f[
 *      a_{n+1} = \frac{1}{\Delta t^2} (x_{n+1} - x_n) +
 *                \frac{-1}{\Delta t} v_n + 0 a_n
 *  \f]
 *  \tparam NODE    Type of node
 */
template<typename NODE> 
class corlib::BackwardEulerSolutionDistributor
    : public corlib::NewmarkSolutionDistributor<NODE>
{
public:
    //! Constructor
    //!
    //! \param[in]   dt       time step size
    BackwardEulerSolutionDistributor( const double dt ) 
    : corlib::NewmarkSolutionDistributor<NODE>( 1./dt,         0.0, 0.0,                
                                                1./(dt*dt), -1./dt, 0.0 )
    { }
};

//------------------------------------------------------------------------------
#endif
