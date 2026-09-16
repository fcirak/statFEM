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

#ifndef solid_fem_nodedynamic_h
#define solid_fem_nodedynamic_h

#include <solid/fem/NodeStatic.hpp>
//------------------------------------------------------------------------------
namespace solid{
    namespace fem{
        
        template<typename BNODE, unsigned DOF> class NodeDynamic;
    }
}

//------------------------------------------------------------------------------
/** \brief Extension of NodeStatic for elastodynamic problems
 *  \details This object stores, in addition to NodeStatic, nodal velocities and
 *  accelerations. Both quantities can be set and queried. 
 *  \tparam BNODE Type of basis node
 *  \tparam DOF   Number of degrees of freedom
 */
template<typename BNODE, unsigned DOF = BNODE::dim>
class solid::fem::NodeDynamic : public solid::fem::NodeStatic<BNODE,DOF>
{
public:
    //--------------------------------------------------------------------------
    //! @name Convenience typedefs
    //@{
    typedef solid::fem::NodeStatic<BNODE,DOF> NodeStatic;
    typedef typename NodeStatic::VecDof           VecDof;
    //@}

    //! Empty constructor, initializes fields
    NodeDynamic( )
        : veloc_( ublas::zero_vector<double>( NodeStatic::dof ) ),
          accel_( ublas::zero_vector<double>( NodeStatic::dof ) ) 
    { }

    //--------------------------------------------------------------------------
    //! @name Accessors 
    //@{

    //! Query nodal velocities
    VecDof giveVelocities(     ) const { return veloc_; }
    
    //! Query nodal accelerations
    VecDof giveAccelerations(  ) const { return accel_; }
    //@}

    //--------------------------------------------------------------------------
    //! @name Mutators 
    //@{

    //! Set nodal velocities
    void setVelocities(    const VecDof & v ) { veloc_ = v; return; }
    
    //! Set nodal accelerations
    void setAccelerations( const VecDof & a ) { accel_ = a; return; }
    //@}

private:
    VecDof  veloc_;      //!< nodal velocities
    VecDof  accel_;      //!< nodal acclerations
};
//------------------------------------------------------------------------------
#endif
