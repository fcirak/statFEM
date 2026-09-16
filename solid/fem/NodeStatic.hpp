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

#ifndef solid_fem_nodestatic_h
#define solid_fem_nodestatic_h
//------------------------------------------------------------------------------
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>
#include <vector>
#include <utility>

//------------------------------------------------------------------------------
namespace solid {
    namespace fem {
        template<typename BNODE, unsigned DOF> class NodeStatic;
        namespace ublas = boost::numeric::ublas;
    }
}

//------------------------------------------------------------------------------
/** \brief Node for static calculations of elastic solids
 *  \details Stores the nodal displacements 
 *  \tparam BNODE Type of Basis node
 *  \tparam DOF   number of degrees of freedom
 */
template<typename BNODE, unsigned DOF = BNODE::dim>
class solid::fem::NodeStatic : public BNODE
{
public:
    //! Typedef of basis class
    typedef          BNODE BasisNode;

    //! @name Basic attributes
    //@{
    static const unsigned dim = BasisNode::dim;
    static const unsigned dof     = DOF; //!< number of dofs
    //@}

    //! Convenience typedef
    typedef ublas::bounded_vector<double,dof> VecDof;
 
    //! Empty constructor, initializes fields
    NodeStatic( ) 
        : disp_(      ublas::zero_vector<double>( dof ) ),
          increment_( ublas::zero_vector<double>( dof ) ),
          force_(     ublas::zero_vector<double>( dof ) ) 
    { }

    //--------------------------------------------------------------------------
    //! @name Dof numbering stuff

    //@{
    //! assign numbers to the degrees of freedom
    void numberDOFs( unsigned & counter );

    //! query number of degrees of freedom
    unsigned giveNumDofs() const { return dof; };

    //! copy the dof array given an iterator
    void copyDofArray( std::vector<unsigned> & dofIndices ) const;
    //@}

    //--------------------------------------------------------------------------
    //! @name Actions on displacements and increment
    //@{
    //! return the current displacements
    VecDof giveDisplacements( ) const { return (disp_ + increment_); }

    //! Add increment onto displacements and clear increment
    void updateDisplacements( );

    //! Return displacement increments
    VecDof giveIncrement( ) const { return increment_; }

    //! Return displacement increments
    void setDisplacements( const VecDof & disp ) { disp_ = disp; return; }

    //! Clear increments
    void clearIncrement( ) { increment_.clear(); return; }

    /** Set the nodal displacement increment
     *  \param[in] inc  New value of the nodal increment  */
    void setIncrement( const VecDof & inc ) { increment_ = inc; return; }

    /** Add to the the nodal displacement increment
     *  \param[in] inc  Value added to the nodal increment  */
    void addToIncrement( const VecDof & inc ) { increment_ += inc; return; }
    //@}

    //--------------------------------------------------------------------------    
    //! @name Actions on forces
    //@{

    //! clear the array of forces
    void clearForce() { force_.clear(); return; }

    //! return the forces
    VecDof giveForce() const { return force_; }

    //! add to the nodal force given some iterators
    void addToForce( const VecDof & inc ) { force_ += inc; return; }
    //@}

    //--------------------------------------------------------------------------
    //! @name Actions on constraints
private:
    typedef std::vector< std::pair<unsigned,double> > ConstraintVec_;
    
public:
    //@{

    //! clear constraints
    void clearConstraints() { constraints_.clear(); return; }

    //! store the constraints in a <dof,value>-pair
    void storeConstraint( const unsigned & component, const double & value );

    //! return constraints
    void giveConstraints( ConstraintVec_ & globalConstraints ) const;
    //@}

protected:
    //--------------------------------------------------------------------------
    VecDof          disp_;           //!< converged displacement components
    VecDof          increment_;      //!< current increment of displacements
    VecDof          force_;          //!< nodal force (internal & external)

    boost::array<unsigned,dof> dofIndices_;  //!< indices of the d. o. f.
    ConstraintVec_             constraints_; //!< <dof,value>-pairs for prescribed dofs
};


#include "NodeStatic.ipp"

#endif
