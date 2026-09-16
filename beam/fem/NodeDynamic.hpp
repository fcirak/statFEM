// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//
//                 Computational Structural Mechanics Lab
//                         University of Cambridge
//
// This software is copyrighted and all rights are retained by the University of
// Cambridge. All use, disclosure, and/or reproduction of any part not expressly
// authorised by the University of Cambridge is prohibited. (C) 2010
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

//! @file NodeDynamic.hpp

#ifndef beam_fem_nodedynamic_h
#define beam_fem_nodedynamic_h

#include <Eigen/Core>
#include <beam/fem/NodeStatic.hpp>
//------------------------------------------------------------------------------
namespace beam{
    namespace fem{

        template< unsigned DIM, unsigned DOF >
        class NodeDynamic;

    }
}

//------------------------------------------------------------------------------
/** \brief Extension of NodeStatic for elastodynamic problems
 *  \details This object stores, in addition to NodeStatic, nodal velocities and
 *  accelerations. Both quantities can be queried. Additional functionality is
 *  the computation of new velocities and acclerations given the current 
 *  displacement increment.
 *  \tparam DIM  Spatial dimension of the problem
 *  \tparam DOF  Number of degrees of freedom
 */
template<unsigned DIM, unsigned DOF>
class beam::fem::NodeDynamic : public beam::fem::NodeStatic<DIM,DOF>
{
private:
    typedef beam::fem::NodeStatic< DIM, DOF >           NodeStatic_;

public:
    typedef typename NodeStatic_::VecDof                VecDof;

    NodeDynamic( ) : veloc_( Eigen::VectorXd::Zero( NodeStatic_::dof ) ),
                     accel_( Eigen::VectorXd::Zero( NodeStatic_::dof ) ) { }

    /** @name Accessors */
    //@{
    //! Query nodal velocities
    VecDof giveVelocities(     ) const { return veloc_; }
    
    //! Query nodal accelerations
    VecDof giveAccelerations(  ) const { return accel_; }
    //@}

    /** @name Mutators */
    //@{
    //! Set nodal velocities
    void setVelocities(    const VecDof & v ) { veloc_ = v; return; }

    //! Set nodal accelerations
    void setAccelerations( const VecDof & a ) { accel_ = a; return; }

    //@}

protected:
    VecDof  veloc_;      //!< nodal velocities
    VecDof  accel_;      //!< nodal acclerations
};

//------------------------------------------------------------------------------
#endif
