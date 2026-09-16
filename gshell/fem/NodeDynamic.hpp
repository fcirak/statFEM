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
//! @date   2011

#ifndef gshell_fem_nodedynamic_h
#define gshell_fem_nodedynamic_h

#include <gshell/fem/NodeStatic.hpp>

//------------------------------------------------------------------------------
namespace gshell{
    namespace fem{

        template< typename BNODE, unsigned DOF >
        class NodeDynamic;

    }
}

//------------------------------------------------------------------------------
//! Extension of NodeStatic to dynamic shell problems
//!
//! This object stores, in addition to NodeStatic, nodal velocities and
//! accelerations. Both quantities can be queried. Additional functionality is
//! the computation of new velocities and acclerations given the current 
//! displacement increment.
//!
//! \tparam DOF  Number of degrees of freedom
template< typename BNODE, unsigned DOF >
class gshell::fem::NodeDynamic : public gshell::fem::NodeStatic< BNODE, DOF >
{
protected:
    typedef gshell::fem::NodeStatic< BNODE, DOF >       NodeStatic_;

public:
    static const unsigned dim = NodeStatic_::dim;

    typedef typename NodeStatic_::VecDim                VecDim;
    typedef typename NodeStatic_::VecDof                VecDof;

public:
    //! Constructor
    NodeDynamic() :
        NodeStatic_(),
        veloc_( Eigen::VectorXd::Zero( NodeStatic_::dof ) ),
        accel_( Eigen::VectorXd::Zero( NodeStatic_::dof ) ),
        mass_( 0. )
    { }

    //! @name Repeated functions to enable binding
    //!
    //! <b>NOTE:</b> The problem arises at instantiation of #corlib::distributorFun
    //!              with #NodeStatic_::copyDofArray and #NodeDynamic_::setAccelerations, i.e.
    //!              'mixed' member functions
    //@{
    void copyDofArray( std::vector< unsigned > & dofIndices ) const { this->NodeStatic_::copyDofArray( dofIndices ); }

    void addToIncrement( const VecDof & inc ) { this->NodeStatic_::addToIncrement( inc ); }

    void setIncrement( const VecDof & inc ) { this->NodeStatic_::setIncrement( inc ); }

    VecDof giveForce( ) const { return this->NodeStatic_::giveForce( ); }
    //@}

    //! @name Accessors
    //@{
    //! Query nodal velocities
    VecDof giveVelocities( ) const { return veloc_; }
    //! Query sub-range of nodal velocities
    VecDim giveVelocitiesRange( const Eigen::VectorXi & r ) const
    {
        FTL_VERIFY( r.size( ) <= dim );
        const VecDof vel = this->giveVelocities( );
        //VecDim v = vel( r );
        VecDim v;
        for ( unsigned i = 0; i < r.size( ); ++i )
            v( i ) = vel( r( i ) );

        return v;
    }
    
    //! Query nodal accelerations
    VecDof giveAccelerations( ) const { return accel_; }
    //! Query sub-range of nodal velocities
    VecDim giveAccelerationsRange( const Eigen::VectorXi & r ) const
    {
        const VecDof acc = this->giveAccelerations( );
        FTL_VERIFY( r.size( ) <= dim );
        //VecDim a = acc( r );
        VecDim a;
        for ( unsigned i = 0; i < r.size( ); ++i )
            a( i ) = acc( r( i ) );
        return a;
    }
    //@}

    //! @name Mutators
    //@{
    //! Set nodal velocities
    void setVelocities( const VecDof & v ) { veloc_ = v; return; }

    //! Set nodal accelerations
    void setAccelerations( const VecDof & a ) { accel_ = a; return; }
    //@}

    //! @name Methods for concentrated mass
    //@{
    //! Set point mass
    void setMass( const double & mass ) { mass_ = mass; return; }

    //! Add to point mass
    void addMass( const double & mass ) { mass_ += mass; return; }

    //! Return point mass
    double giveMass( ) const { return mass_; }
    //@}

protected:
    using NodeStatic_::disp_; 
    using NodeStatic_::increment_;

    VecDof  veloc_;      //!< nodal velocities
    VecDof  accel_;      //!< nodal acclerations

    //! Point mass of node
    //!
    //! Note: This concentrated point mass is used with explicit time integration only.
    double mass_;

};

//------------------------------------------------------------------------------
#endif
