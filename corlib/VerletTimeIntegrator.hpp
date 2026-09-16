// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//
//                 Computational Structural Mechanics Lab
//                         University of Cambridge
//
// This software is copyrighted and all rights are retained by the University of
// Cambridge. All use, disclosure, and/or reproduction of any part not expressly
// authorised by the University of Cambridge is prohibited. (C) 2010
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

//! @file VerletTimeIntegrator.hpp

#ifndef corlib_verlettimeintegrator_h
#define corlib_verlettimeintegrator_h

#include <functional>

//------------------------------------------------------------------------------
namespace corlib{

    template< typename NODE >
    class VerletPredictor;
    
    template< typename NODE >
    class VerletCorrector;
    
    template< typename NODE >
    class VerletInitialisor;
    
    namespace ublas = boost::numeric::ublas;
}

//==============================================================================
/// Predictor of Verlet (explicit Newmark) method
template< typename NODE >
class corlib::VerletPredictor
    : public boost::function<void(NODE*)>
{
private:
    ///
    typedef typename NODE::VecDof                   VecDof_;

public:
    /// Constructor
    VerletPredictor( const double & timeStepSize,
                     const double gamma = 0.5 )
        : timeStepSize_( timeStepSize ),
          facDA_( 0.5 * timeStepSize_ * timeStepSize_ ),
          facVA_( (1.0 - gamma) * timeStepSize_ )
    { }

    /// New displacement increment
    ///\f[
    ///     \Delta\vec{d} = \vec{d}_{n+1} - \vec{d}_{n+1}
    ///                   = \Delta{t} \vec{v}_n + \frac{\Delta{t}^2}{2} \vec{a}_n
    ///\f]
    /// and mid-velocities
    ///\f[
    ///     \vec{v}_{n+1/2} = \vec{v}_n
    ///                     + (1-\gamma)\Delta{t}\vec{a}_{n}
    ///     \quad\mbox{with}\quad \gamma=1/2
    ///\f]
    /// with displacements \f$\vec{d}\f$, velocities \f$\vec{v}\f$
    /// and accelerations \f$\vec{a}\f$ and time step size \f$\Delta{t}=t_{n+1}-t_{n}\f$.
    void operator()( NODE * node )
    {
        const VecDof_ dispIncr = 
            timeStepSize_ * node->giveVelocities( ) +
            facDA_ * node->giveAccelerations( );
        node->setIncrement( dispIncr );

        const VecDof_ midVel =
            node->giveVelocities( ) +
            facVA_ * node->giveAccelerations( );
        node->setVelocities( midVel );
    }

private:
    /// Time step size
    double timeStepSize_;

    /// Factor by which accelerations \f$\vec{a}_n\f$ contribute to displacement
    double facDA_;
    /// Factor by which accelerations \f$\vec{a}_n\f$ contribute to velocities
    double facVA_;
};


//==============================================================================
/// Corrector of Verlet (explicit Newmark) method
template< typename NODE >
class corlib::VerletCorrector
    : public boost::function<void(NODE*)>
{
private:
    ///
    typedef typename NODE::VecDof                   VecDof_;

public:
    /// Constructor
    VerletCorrector( const double & timeStepSize,
                     const double gamma = 0.5,
                     const double damping = 0.0 )
        : timeStepSize_( timeStepSize ),
          facVA_( gamma * timeStepSize_ ),
          damping_( damping )
    { }

    /// Update of state variables
    /// 
    /// The new accelerations are determined by inversion with nodal point mass.
    ///\f[
    ///     \vec{a}_{n+1} = -\frac{1}{m}\Big(
    ///                       + f_{int}(d_{n+1}) 
    ///                       - f_{ext}(t_{n+1})
    ///                       + c * v_{n+1/2}
    ///                   \Big)
    ///\f]
    /// with nodal point mass \f$m\f$ and damping \f$c\f$.
    ///
    /// The displacements
    ///\f[
    ///     \vec{d}_{n+1} = \vec{d}_n + \Delta\vec{d}
    ///\f]
    /// and velocities
    ///\f[
    ///     \vec{v}_{n+1} = \vec{v}_{n+1/2} + \gamma \Delta{t} \vec{a}_{n+1}
    ///\f]
    void operator()( NODE * node )
    {
        const double nodeMass = node->giveMass( );
        const VecDof_ accNew = ( node->giveForce( ) -
                                 damping_ * node->giveVelocities( ) ) / nodeMass;
        const VecDof_ velNew = node->giveVelocities( ) + facVA_ * accNew;

        // update state
        node->updateDisplacements( );
        node->setVelocities( velNew );
        node->setAccelerations( accNew );
    }

private:
    /// Time step size
    double timeStepSize_;

    /// Factor by which accelerations \f$\vec{a}_{n+1}\f$ contribute to velocities
    double facVA_;

    /// Node-wise diagonal damping factor
    double damping_;
};

//==============================================================================
/// Initialisor of Verlet (explicit Newmark) method
template< typename NODE >
class corlib::VerletInitialisor
    : public boost::function<void(NODE*)>
{
private:
    ///
    typedef typename NODE::VecDof                   VecDof_;

public:
    /// Constructor
    VerletInitialisor( const double damping = 0.0 )
        : damping_( damping )
    { }

    
    /// Update of accelerations
    /// 
    /// The new accelerations are determined by inversion with nodal point mass.
    ///\f[
    ///     \vec{a}_{n+1} = - \frac{1}{m}\Big(
    ///                       + f_{int}(d_{n+1}) 
    ///                       - f_{ext}(t_{n+1})
    ///                       + c * v_{n+1/2}
    ///                   \Big)
    ///\f]
    /// with nodal point mass \f$m\f$ and damping \f$c\f$.
    void operator()( NODE * node )
    {
        const double nodeMass = node->giveMass( );
        const VecDof_ accNew = ( node->giveForce( ) -
                                 damping_ * node->giveVelocities( ) ) / nodeMass;
        node->setAccelerations( accNew );
    }

private:
    /// Node-wise diagonal damping factor
    double damping_;
};

#endif
