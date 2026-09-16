// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//
//                 Computational Structural Mechanics Lab
//                         University of Cambridge
//
// This software is copyrighted and all rights are retained by the University of
// Cambridge. All use, disclosure, and/or reproduction of any part not expressly
// authorised by the University of Cambridge is prohibited. (C) 2010
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

//! @file NodeStatic.hpp

#ifndef beam_fem_nodestatic_h
#define beam_fem_nodestatic_h
//------------------------------------------------------------------------------
#include <vector>
#include <utility>
#include <Eigen/Core>
#include <corlib/eigenX.hpp>
#include <corlib/NodeBasic.hpp>

//------------------------------------------------------------------------------
namespace beam {
    namespace fem {

        template< unsigned DIM, unsigned DOF >
        class NodeStatic;

        namespace eigenX = corlib::eigenX;

    }
}

//------------------------------------------------------------------------------
/** \brief Node for static calculations 
 *
 *  Stores the nodal displacements (or similar physical quantity)
 *  \tparam DIM the dimension
 *  \tparam DOF number of degrees of freedom
 */
template< unsigned DIM, unsigned DOF>
class beam::fem::NodeStatic : public corlib::NodeBasic< DIM >
{
public:
    static const unsigned dof     = DOF;        //!< number of dofs
    static const unsigned dim     = DIM;        //!< embedding space
    typedef eigenX::VectorSd<dof>             VecDof;

protected:
    typedef std::array< unsigned, dof >                    IndexArray_;
    typedef std::vector< std::pair<unsigned,double> >      ConstraintVec_;
    
public:
    //! Empty Constructor 
    NodeStatic( ) 
        : disp_(      Eigen::VectorXd::Zero( dof ) ),
          increment_( Eigen::VectorXd::Zero( dof ) ),
          force_(     Eigen::VectorXd::Zero( dof ) )
    {}

    //! assign numbers to the degrees of freedom
    void numberDOFs( unsigned & counter );

    //! copy the dof array given an iterator
    void copyDofArray( std::vector<unsigned> & dofIndices ) const;

    //! @name Actions on displacements
    //@{
    //! clear displacement vector
    void clearDisplacements() { disp_.setZero(); return; }
    //! return the current displacements
    VecDof giveDisplacements( ) const { return (disp_ + increment_); }
    //! Add increment onto displacements and clear increment
    void updateDisplacements( );
    //! Set last converged displacements
    void setDisplacements( const VecDof & disp ) { disp_ = disp; return; }

    //! Clear increments
    void clearIncrement( ) { increment_.setZero(); return; }
    //! Return displacement increments
    VecDof giveIncrement( ) const { return increment_; }
    /** Set the nodal displacement increment
     *  \param[in] inc  New value of the nodal increment  */
    void setIncrement( const VecDof & inc ) { increment_ = inc; return; }
    /** Add to the the nodal displacement increment
     *  \param[in] inc  Value added to the nodal increment  */
    void addToIncrement( const VecDof & inc ) { increment_ += inc; return; }
    //@}
    
    //! @name Actions on forces
    //@{
    //! clear the array of forces
    void clearForce() { force_.setZero(); return; }
    //! return the forces
    VecDof giveForce() const { return force_; }
    //! add to the nodal force given some iterators
    void addToForce( const VecDof & inc ) { force_ += inc; return; }
    //@}

    //! clear constraints
    void clearConstraints() { constraints_.clear(); return; }

    //! store the constraints in a <dof,value>-pair
    void storeConstraint( const unsigned & component, const double & value )
    {
        if (component < dof)
            constraints_.push_back( std::make_pair( component, value ) );
        return;
    }

    //! return constraints
    void giveConstraints( ConstraintVec_ & globalConstraints ) const
    {
        ConstraintVec_::const_iterator cIter = constraints_.begin();
        ConstraintVec_::const_iterator cEnd  = constraints_.end( );
        for( ; cIter != cEnd; cIter ++ ) {
            const unsigned index = dofIndices_[ cIter -> first ];
            globalConstraints.push_back( std::make_pair( index,
                                                         cIter -> second ) );
        }
        return;
    }
    //@}

protected:
    VecDof          disp_;           //!< converged displacement components
    VecDof          increment_;      //!< current increment of displacements
    VecDof          force_;          //!< nodal force
    IndexArray_     dofIndices_;     //!< indices of the degrees of freedom
    ConstraintVec_  constraints_;    //!< <dof,value>-pairs for prescribed dofs
};




#include "NodeStatic.ipp"

#endif
