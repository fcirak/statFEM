// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//
//                 Computational Structural Mechanics Lab
//                         University of Cambridge
//
// This software is copyrighted and all rights are retained by the University of
// Cambridge. All use, disclosure, and/or reproduction of any part not expressly
// authorised by the University of Cambridge is prohibited. (C) 2010
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

//! @file NodePotential.hpp

#ifndef del2_fem_nodepotential_h
#define del2_fem_nodepotential_h
//------------------------------------------------------------------------------
#include <vector>
#include <utility>
#include <corlib/eigenX.hpp>

//------------------------------------------------------------------------------
namespace del2 {
    namespace fem {

        template<typename BNODE> class NodePotential;

        namespace eigenX = corlib::eigenX;
    }
}

//------------------------------------------------------------------------------
/**  \brief Finite element node for Poisson equation and colleagues
 *   \details This object has one (scalar) degree of freedom, it manages its
 *   storage, stores an increment for non-linear computations, holds the nodal
 *   force corresponding to that dof and can store a constraint.
 */
template<typename BNODE>
class del2::fem::NodePotential : public BNODE
{
public:
    //! Typedef of basis class
    typedef          BNODE BasisNode;

    //! @name Basic attributes
    //@{
    static const unsigned dim = BasisNode::dim;
    static const unsigned dof = 1; //!< number of dofs
    //@}

    //! Convenience typedef
    typedef eigenX::VectorSd<dof> VecDof;

    //! Empty constructor, initializes fields
    NodePotential( )
        : pot_( VecDof::Zero() ),
          increment_( VecDof::Zero() ),
          force_( VecDof::Zero() ),
          dofIndex_( std::numeric_limits<unsigned>::max() )
    { }

    //--------------------------------------------------------------------------
    //! @name Dof numbering functions
    //@{

    //! Asign numbers to the degrees of freedom
    void numberDOFs( unsigned & counter ) { dofIndex_ = counter++; }

    //! Query number of degrees of freedom
    unsigned giveNumDofs() const { return dof; };

    //! copy the index of the degree of freedom
    void copyDofArray( std::vector<unsigned> & dofIndices ) const
    {
        dofIndices.resize( 1 );
        dofIndices[0] = dofIndex_;
        return;
    }
    //@}

    //--------------------------------------------------------------------------
    //! @name Potential and increment

    //@{
    //! Set the nodal potential
    void setPotential( const VecDof & pot ) { pot_ = pot; return; }

    //! Add increment onto potential
    void addToPotential() { pot_ += increment_; }

    //! Return the nodal potential
    VecDof getPotential() const { return pot_; }

    //! Set potential equal to increment
    void updatePotential() { pot_ = increment_; increment_.setZero(); }

    /** Set the nodal displacement increment
     *  \param[in] inc  New value of the nodal increment  */
    void setIncrement( const VecDof & inc ) { increment_ = inc; return; }

    /** Add to the the nodal displacement increment
     *  \param[in] inc  Value added to the nodal increment  */
    void addToIncrement( const VecDof & inc ) { increment_ += inc; return; }
    //@}

    //--------------------------------------------------------------------------
    //! @name Nodal force

    //@{
    //! clear the array of forces
    void clearForce() { force_.setZero() ; return; }

    //! add to the nodal force given some iterators
    void addToForce( const VecDof & inc ) { force_ += inc; return; }

    //! return the forces
    VecDof getForce() const { return force_; }
    //@}

    //--------------------------------------------------------------------------
    //! @name Constraint handling routines
    //@{
private:
    typedef std::vector< std::pair<unsigned,double> > ConstraintVec_;

public:
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
        if ( !constraints_.empty() ) {
            std::pair<unsigned, double> cstr( dofIndex_, constraints_[0].second );
            globalConstraints.push_back( cstr );
        }
        return;
    }
    //@}

protected:
    //--------------------------------------------------------------------------
    VecDof          pot_;            //!< converged potential
    VecDof          increment_;      //!< current increment
    VecDof          force_;          //!< nodal force
    unsigned        dofIndex_;       //!< indices of the degrees of freedom
    ConstraintVec_  constraints_;    //!< <dof,value>-pairs for prescribed dofs

};

#endif
