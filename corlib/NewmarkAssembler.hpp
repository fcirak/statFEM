// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//
//                 Computational Structural Mechanics Lab
//                         University of Cambridge
//
// This software is copyrighted and all rights are retained by the University of
// Cambridge. All use, disclosure, and/or reproduction of any part not expressly
// authorised by the University of Cambridge is prohibited. (C) 2010
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

//! @file NewmarkAssembler.hpp

#ifndef corlib_newmarkassembler_h
#define corlib_newmarkassembler_h
//------------------------------------------------------------------------------
//! standard library includes
#include <algorithm>
//! boost includes
#include <boost/function.hpp>
//! Eigen includes
#include <Eigen/Core>
//! corlib includes
#include <corlib/eigenX.hpp>
#include <corlib/verify.hpp>

//------------------------------------------------------------------------------
namespace corlib {

    namespace eigenX = corlib::eigenX;

    namespace detail_{

        //----------------------------------------------------------------------
        //! \cond SKIPDOX
        //! Dummy for default template parameters
        template<typename ELEMENT>
        struct DoNothing
            : public boost::function<void( const ELEMENT*, 
                                           Eigen::MatrixXd& )>
        {
            void operator()( const ELEMENT* ep, 
                             Eigen::MatrixXd& mat ) const
            {
                return;
            }
        };
        //! \endcond
    }

    //--------------------------------------------------------------------------
    //! Forward declarations
    template<typename SYSSOLVER,typename ELEMENT, 
             typename MASS, typename STIFF>
    class NewmarkMatrixAssembler;

    template< typename SYSSOLVER, typename ELEMENT, typename MASS, 
              typename STIFF = detail_::DoNothing<ELEMENT> >
    class NewmarkInertiaForceAssembler;

    template<typename NODE>
    class NewmarkSolutionDistributor;
}

//------------------------------------------------------------------------------
//! \brief Assemble the dynamic stiffness for the Newmark method
//! \details Using a Newmark time integration routine, a linear system arises
//! of the form
//! \f[
//!      \left[ \frac{1}{\beta \Delta t^2} M + \frac{\gamma}{\beta \Delta t} C
//!             + K \right] ( x_{n+1}^{(k+1)} - x_{n+1}^{(k)} ) = RHS
//! \f]
//! Under the assumption of Rayleigh damping, i.e., \f$ C = d_1 M + d_2 K \f$, 
//! the system can be simplified.
//! \tparam SYSSOLVER Type of system solver to receive the final system matrix
//! \tparam ELEMENT   Type of finite element which delivers the matrices
//! \tparam MASS      A functor which gives the element mass matrix
//! \tparam STIFF     A functor which gives the element's linearised stiffness
template<typename SYSSOLVER, typename ELEMENT, 
         typename MASS, typename STIFF>
class corlib::NewmarkMatrixAssembler
    : public boost::function<void (const ELEMENT*, SYSSOLVER*) >
{
public:
    //--------------------------------------------------------------------------
    //! Constructor
    //!
    //! \param[in]   beta     in (0,1/2]
    //! \param[in]   gamma    in (0,1]
    //! \param[in]   dt       time step size 
    //! \param[in]   mass     Functor for the mass matrix computation
    //! \param[in]   stiff    Functor for the stiffness matrix computation
    //! \param[in]   d1, d2   Rayleigh damping: C = d1 * M + d2 * K
    NewmarkMatrixAssembler( const double beta,
                            const double gamma,
                            const double dt,
                            MASS  & mass,
                            STIFF & stiff, 
                            const double d1 = 0., 
                            const double d2 = 0.) 
        : factorM_( 1. / ( beta * dt * dt ) + 
                    d1 * gamma / (beta*dt) ), //!< Add Rayleigh damping
          factorK_( 1.  + 
                    d2 * gamma / (beta*dt) ), //!< Add Rayleigh damping
          mass_(  mass  ),
          stiff_( stiff )
    { 
        FTL_VERIFY( (gamma > 0.) and (gamma <=1.0) );
        FTL_VERIFY( (beta  > 0.) and (beta  <=0.5) );
    }

protected:
    //--------------------------------------------------------------------------
    //! Constructor which allows direct setting of the factors
    //!
    //! \param[in]   factorM    Factor for mass matrix
    //! \param[in]   factorK    Factor for stiffness matrix
    //! \param[in]   mass       Functor for the mass matrix computation
    //! \param[in]   stiff      Functor for the stiffness matrix computation
    NewmarkMatrixAssembler( const double factorM,
                            const double factorK,
                            MASS &  mass,
                            STIFF & stiff) 
        : factorM_( factorM ),
          factorK_( factorK ),
          mass_(    mass ),
          stiff_(   stiff )
    { }

public:
    //--------------------------------------------------------------------------
    //! Function call applied to element and system solver pointers
    //!
    //! \param[in] ep     Pointer to element
    //! \param[in] target Pointer to system solver
    void operator( ) ( const ELEMENT * ep, SYSSOLVER * target ) const
    { 
        // get element ep's degrees of freedom
        std::vector<unsigned> dofs;
        ep -> getDofIndices( dofs );

        // allocate storage for element ep's dynamic stiffness matrix
        Eigen::MatrixXd stiffMat( dofs.size(), dofs.size() );
        stiffMat.setZero();

        // compute element ep's mass matrix
        mass_( ep, stiffMat );
        stiffMat *= (factorM_ / factorK_);// denominator will be canceled below

        // compute and add element ep's stiffness matrix
        stiff_( ep, stiffMat );
        stiffMat *= factorK_;

        // assemble element ep's stiffness matrix
        target -> insertToMatrix( stiffMat, dofs, dofs );
                                    
        return;
    }

private:
    const double factorM_; //!< Multiplier for the element's mass matrix
    const double factorK_; //!< Multiplier for the element's stiffness matrix

    MASS  &      mass_;    //!< Object to compute/obtain the mass matrix
    STIFF &      stiff_;   //!< Obtain to compute/obtain the stiffness matrix
};

//==============================================================================


//------------------------------------------------------------------------------
/** \brief Inertia forces for the Newmark time stepping algorithm
 *
 *  \details The complete RHS in case of a Newmark time integration scheme reads
 *  \f[
 *     RHS = f^{ext}_{n+1} - f^{int}(x_{n+1}^{(k)})
 *      - M \left[ \frac{1}{\beta \Delta t^2} (x^{(k)}_{n+1} - x_n) +
 *                 \frac{-1}{\beta \Delta t} v_n + (1-1/(2\beta)) a_n \right]
 *      - C \left[ \frac{\gamma}{\beta \Delta t} (x^{(k)}_{n+1} - x_n) +
 *                 (1-\gamma/beta) v_n + \Delta t (1-\gamma/(2\beta)) a_n \right]
 *  \f]
 *  It is assumed the the external and internal forces, \f$f^{ext}\f$ and
 *  \f$ f^{int} \f$ are computed separately. 
 *  \tparam SYSSOLVER  Final storage of the computed result
 *  \tparam ELEMENT Type of element to give the data
 *  \tparam MASS      A functor which gives the element mass matrix
 *  \tparam STIFF     A functor which gives the element's linearised stiffness
 */
template<typename SYSSOLVER, typename ELEMENT, typename MASS, typename STIFF>
class corlib::NewmarkInertiaForceAssembler
    : public boost::function<void (const ELEMENT*, SYSSOLVER*) >
{
public:
    //--------------------------------------------------------------------------
    //! Constructor
    //!
    //! \param[in]   beta    Newmark beta coefficient
    //! \param[in]   gamma   Newmark gamma coefficient
    //! \param[in]   dt      time step size
    //! \param[in]   mass    Functor for the mass matrix computation
    //! \param[in]   stiff   (Optional) stiffness functor (needed for damping)
    //! \param[in]   d1,d2   Rayleigh damping: C = d1 * M + d2 * K
    NewmarkInertiaForceAssembler( const double beta,
                                  const double gamma, 
                                  const double dt, 
                                  MASS  & mass, 
                                  STIFF   stiff = STIFF(), 
                                  const double d1 = 0.,
                                  const double d2 = 0. ) 
        : facVD_(  gamma/beta/dt ),
          facVV_(  1. - gamma/beta ),
          facVA_(  1. - gamma/2./beta ),
          facAD_(  1./(beta*dt*dt) ),
          facAV_( -1./(beta*dt) ),
          facAA_(  1. - 0.5/beta ),
          mass_(  mass ),
          stiff_( stiff ), 
          damp1_( d1 ), 
          damp2_( d2 )
    { 
        FTL_VERIFY( (gamma > 0.) and (gamma <=1.0) );
        FTL_VERIFY( (beta  > 0.) and (beta  <=0.5) );
    }

protected:
    //--------------------------------------------------------------------------
    //! Constructor which allows direct setting of the factors
    //!
    //! \param[in]   facVD    Factor for displacement increment
    //! \param[in]   facVV    Factor for velocity
    //! \param[in]   facVA    Factor for acceleration
    //! \param[in]   facAD    Factor for displacement increment
    //! \param[in]   facAV    Factor for velocity
    //! \param[in]   facAA    Factor for acceleration
    //! \param[in]   mass     Functor for the mass matrix computation
    //! \param[in]   stiff    (Optional) stiffness functor (needed for damping)
    //! \param[in]   d1,d2    Rayleigh damping: C = d1 * M + d2 * K
    NewmarkInertiaForceAssembler( const double facVD,
                                  const double facVV,
                                  const double facVA,
                                  const double facAD,
                                  const double facAV,
                                  const double facAA,
                                  MASS  & mass,
                                  STIFF   stiff = STIFF(),
                                  const double d1 = 0.,
                                  const double d2 = 0. )
        : facVD_( facVD ), facVV_( facVV ), facVA_( facVA ),
          facAD_( facAD ), facAV_( facAV ), facAA_( facAA ),
          mass_(  mass  ), stiff_( stiff ),
          damp1_( d1 ), damp2_( d2 )
    { }
    
public:
    void operator() ( const ELEMENT * ep, SYSSOLVER * target )
    {
        // get element nodal data
        // current modal quantities at target time t_n
        typedef typename ELEMENT::MatDofNF    MatDofNF;
        MatDofNF disInc;  ep -> nodalIncrements( disInc ); // d_{n+1}^k - d_n
        MatDofNF vel;     ep -> nodalVelocities(    vel ); // v_n
        MatDofNF acc;     ep -> nodalAccelerations( acc ); // a_n

        // sizes
        const unsigned dof = ELEMENT::dof;
        const unsigned matSize = dof * disInc.cols( );

        // compute element's Newmark velocities
        Eigen::VectorXd veloc( matSize );
        for ( unsigned v = 0; v < disInc.cols( ); v ++ ) {
            veloc.segment( v*dof, dof ) =
                facVD_ * disInc.col( v ) +
                facVV_ * vel.col( v ) +
                facVA_ * acc.col( v );
        }

        // compute element's Newmark accelerations
        Eigen::VectorXd accel( matSize );
        for ( unsigned v = 0; v < disInc.cols( ); v ++ ) {
            accel.segment( v*dof, dof ) =
                facAD_ * disInc.col( v ) +
                facAV_ * vel.col( v ) +
                facAA_ * acc.col( v );
        }

        // get element ep's degrees of freedom
        std::vector<unsigned> rows;
        ep -> getDofIndices( rows );

        // compute element ep's mass matrix
        Eigen::MatrixXd massMat( rows.size(), rows.size() );
        massMat.setZero();
        mass_( ep, massMat );

        // damping matrix a la Rayleigh: C = damp1 * M + damp2 * K
        Eigen::MatrixXd dampMat = damp1_ * massMat;
        {
            Eigen::MatrixXd stiffMat( rows.size(), rows.size() );
            stiffMat.setZero();
            stiff_( ep, stiffMat );
            dampMat += damp2_ * stiffMat;
        }

        // element ep's (negative) inertial force vector
        const Eigen::VectorXd dynForce =
            - ( massMat * accel )
            - ( dampMat * veloc );

        // assemble inertia forces to system  RHS 
        target -> insertToRhs( dynForce, rows );

        return;
    }

private:
    //! @name Contribution to velocity
    //@{
    const double facVD_;
    const double facVV_;
    const double facVA_;
    //@}

    //! @name Contribution to acceleration
    //@{
    const double facAD_;
    const double facAV_;
    const double facAA_;
    //@}

    //! @name Functors for element matrices
    //@{
    MASS  &      mass_; //! Mass      matrix providing functor
    STIFF        stiff_;//! Stiffness matrix providing functor
    //@}

    //! @name Parameters of the Rayleigh damping
    //@{
    const double damp1_; //!< Mass      matrix contribution
    const double damp2_; //!< Stiffness matrix contribution
    //@}
};

//==============================================================================

//------------------------------------------------------------------------------
/** \brief Update the nodal velocities and accelerations a la Newmark
 *  \details
 *  The update rules of the Newmark time integration scheme are
 *  \f[
 *      v_{n+1} = \frac{\gamma}{\beta \Delta t} (x_{n+1} - x_n) +
 *                (1-\gamma/beta) v_n + \Delta t (1-\gamma/(2\beta)) a_n 
 *  \f]
 *  and
 *  \f[
 *      a_{n+1} = \frac{1}{\beta \Delta t^2} (x_{n+1} - x_n) +
 *                \frac{-1}{\beta \Delta t} v_n + (1- 1/(2\beta)) a_n
 *  \f]
 *  \tparam NODE    Type of node
 */
template<typename NODE> 
class corlib::NewmarkSolutionDistributor
    : public boost::function<void( NODE* )>
{
public:
    //! Constructor of the distributor
    //!
    //! \param[in]   beta     in (0,1/2]
    //! \param[in]   gamma    in (0,1]
    //! \param[in]   dt       time step size 
    NewmarkSolutionDistributor( const double beta, 
                                const double gamma,
                                const double dt )
        : facVD_( gamma / (beta * dt) ),
          facVV_( 1. - gamma/beta ),
          facVA_( dt * (1. - 0.5 * gamma/beta) ),
          facAD_( 1./(beta * dt * dt) ),
          facAV_( - 1. / (beta * dt) ),
          facAA_( 1. - 0.5 / beta ) 
    { 
        FTL_VERIFY( (gamma > 0.) and (gamma <=1.0) );
        FTL_VERIFY( (beta  > 0.) and (beta  <=0.5) );
    }

protected:
    //--------------------------------------------------------------------------
    //! Constructor which allows direct setting of the factors (Euler backward)
    //!
    //! \param[in] facVD Factor for displacement increment to determine velocity
    //! \param[in] facVV Factor for velocity to determine velocity
    //! \param[in] facVA Factor for acceleration to determine velocity
    //! \param[in] facAD Factor for displacement increment to determine acceleration
    //! \param[in] facAV Factor for velocity to determine acceleration
    //! \param[in] facAA Factor for acceleration to determine acceleration
    NewmarkSolutionDistributor( const double facVD,
                                const double facVV,
                                const double facVA,
                                const double facAD,
                                const double facAV,
                                const double facAA )
        : facVD_( facVD ), facVV_( facVV ), facVA_( facVA ),
          facAD_( facAD ), facAV_( facAV ), facAA_( facAA )
    { }

public:
    //--------------------------------------------------------------------------
    /** Given the nodal displacement increment, compute the new values of the 
     *  velocities and acclerations according to the Newmark algorithm:
     *  \f[
     *       v_{n+1} = \frac{\gamma}{\beta \Delta t} \Delta u
     *               + (1-\frac{\gamma}{\beta}) v_{n}
     *               + \Delta t (1-\frac{\gamma}{2\beta}) a_{n}
     *  \qquad
     *       a_{n+1} = \frac{1}{\beta \Delta t^2} \Delta u 
     *               - \frac{1}{\beta\Delta t} v_{n1}
     *               + (1-\frac{1}{2\beta}) a_{n}
     *  \f]
     *  \param[in]  np      Pointer to a (dynamic) node
     */
    void operator()( NODE * np ) const
    {
        typedef typename NODE::VecDof   VecDof_;

        // get node np's current velocities and accelerations
        VecDof_ V( np -> giveVelocities() );
        VecDof_ A( np -> giveAccelerations() );

        // get node np's displacement increment
        const VecDof_ dIncrement( np -> giveIncrement() );

        // compute node np's new velocities and accelerations
        const VecDof_ temp( V );

        V = facVD_ * dIncrement + facVV_ * temp + facVA_ * A;
        A = facAD_ * dIncrement + facAV_ * temp + facAA_ * A;

        // set node np's velocities and accelerations
        np -> setVelocities( V );
        np -> setAccelerations( A );

        return;
    }

private:
    const double facVD_;  //!< \f$ \gamma /(\beta \Delta t   ) \f$
    const double facVV_;  //!< \f$ 1- \gamma/\beta             \f$
    const double facVA_;  //!< \f$ \Delta t(1-\gamma/(2\beta)) \f$
    const double facAD_;  //!< \f$ 1/(\beta \Delta t^2)        \f$
    const double facAV_;  //!< \f$ - 1/(\beta \Delta t)        \f$
    const double facAA_;  //!< \f$ 1 - 1/ (2\beta)             \f$
};
//------------------------------------------------------------------------------
#endif
