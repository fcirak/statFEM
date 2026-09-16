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
//! @date   2010

#ifndef gshell_fem_nodestatic_h
#define gshell_fem_nodestatic_h

//------------------------------------------------------------------------------
#include <vector>
#include <utility>
#include <array>

#include <Eigen/Core>
#include <corlib/eigenX.hpp>

#include <gshell/fem/Config.hpp>

//------------------------------------------------------------------------------
namespace gshell {
    namespace fem {

        template< typename BNODE, unsigned DOF >
        class NodeStatic;

        namespace eigenX = corlib::eigenX;
    }
}

//------------------------------------------------------------------------------
//! Generalised-shell node for static calculations 
//!
//! <h4>In case of <i>generalised</i> (shear-flexible) shell:</h4>
//! Stores the nodal displacements of the mid-surface \f$(u^1,u^2,u^3)\f$ and the shear
//! direction components \f$(w^1,w^2)\f$. In total these are 5 unknows per node. These
//! are stored in the #disp_ vector such that \f$(u^1,u^2,u^3,w^1,w^2)\f$. The
//!  same pattern must be followed implicitly by the nodal test diplacements
//! \f$(\delta{u}^1,\delta{u}^2,\delta{u}^3,\delta{w}^1,\delta{w}^2)\f$
//! 
//! <h4>In case of <i>Kirchhoff-Love</i> (shear-rigid) shell:</h4>
//! Stores 3 nodal displacements of the mid-surface  \f$(u^1,u^2,u^3)\f$.
//!
//! \tparam   DOF   Number of DOFs per node, i.e. 3 or 5 (see comment above)
template< typename BNODE, unsigned DOF >
class gshell::fem::NodeStatic : public BNODE
{
public:
    typedef BNODE                                          BasisNode;

    //! Dimension of embedding space (always 3)
    static const unsigned dim     = BasisNode::dim;
    //! Total number of DOFs per node
    static const unsigned dof     = DOF;  // must be 3 or 5
    //! Size of force vector on node
    static const unsigned vecSize = dof;

    //! Vector of DOF-length
    typedef eigenX::VectorSd<dof>           VecDof;
    typedef eigenX::VectorSd<dim>           VecDim;

protected: 
    typedef NodeStatic< BNODE, DOF >        NodeStatic_;

public: 
    //! constructor with a given id (passed to NodeBasic)
    NodeStatic() :
        BasisNode(),
        disp_(      Eigen::VectorXd::Zero( dof ) ),
        increment_( Eigen::VectorXd::Zero( dof ) ),
        force_(     Eigen::VectorXd::Zero( dof ) ),
        forceExt_(  Eigen::VectorXd::Zero( dof ) ),
        gamma_(     Eigen::VectorXd::Zero( dim ) ),
        aGrad_(     Eigen::VectorXd::Zero( dim ) )
    {
        static_assert( ( dof == 3 ) or ( dof == 5 ) );
    }

    //! Assign numbers to the degrees of freedom
    void numberDOFs( unsigned & counter );

    //! Query number of degrees of freedom
    unsigned giveNumDofs() const { return dof; };

private:
    typedef std::array< unsigned, dof >                    IndexArray_;
    typedef std::vector< std::pair< unsigned, double > >   ConstraintVec_;
    
public:
    //! Copy the dof array given an iterator
    //!
    //! \param[in]  dofIndices Storage of dof indices
    //!             The vector must come _allocated_ to required size.
    void copyDofArray( std::vector< unsigned > & dofIndices ) const;

    //! @name Actions on displacements
    //@{
    //! Clear displacement vector
    void clearDisplacements() { disp_.setZero(); return; }

    //! return the current displacements
    VecDof giveDisplacements( ) const { return ( disp_ + increment_ ); }
    //! Return range of current "displacements" (ie total DOFs per node)
    VecDim giveDisplacementsRange( const Eigen::VectorXi & r ) const
    {
        const VecDof dis = this->giveDisplacements( );
        FTL_VERIFY( r.size( ) <= dim );
        //VecDim d = dis( r );
        VecDim d;
        for( unsigned i = 0; i < r.size( ); ++i )
            d( i ) = dis( r( i ) );

        return d;
    }

    //! Add increment onto displacements and clear increment
    void updateDisplacements( );

    //! Return displacement increments
    VecDof giveIncrement( ) const { return increment_; }

    //! Return displacement increments
    void setDisplacements( const VecDof & disp ) { disp_ = disp; return; }

    //! Clear increments
    void clearIncrement( ) { increment_.setZero(); return; }

    //! Set the nodal displacement increment
    //! \param[in] inc  New value of the nodal increment
    void setIncrement( const VecDof & inc ) { increment_ = inc; return; }

    //! Add to the the nodal displacement increment
    //!  \param[in] inc  Value added to the nodal increment
    void addToIncrement( const VecDof & inc ) { increment_ += inc; return; }
    //@}
    
    //! @name Actions on forces
    //@{
    //! clear the array of forces
    void clearForce() { force_.setZero(); return; }

    //! return the forces
    VecDof giveForce( ) const { return force_; }
    //! Return range of nodal "forces"
    VecDim giveForceRange( const Eigen::VectorXi & r ) const
    {
        const VecDof frc = this->giveForce( );
        FTL_VERIFY( r.size( ) <= dim );
        //VecDim f = frc( r );
        VecDim f;
        for( unsigned i = 0; i < r.size( ); ++i )
            f( i ) = frc( r( i ) );

        return f;
    }

    //! add to the nodal force given some iterators
    void addToForce( const VecDof & inc ) { force_ += inc; return; }

    //! clear array of external forces
    void copyForceToForceExt( ) { forceExt_ = force_; return; }

    //! Return external forces
    VecDof giveForceExt( ) const { return forceExt_; }

    //@}

    //! @name Actions on constraints
    //@{
    //! clear constraints
    void clearConstraints( ) { constraints_.clear(); }

    //! store the constraints in a <dof,value>-pair
    void storeConstraint( const unsigned & component, const double & value );

    //! return constraints
    //!
    //! \param[in,out]   globalConstraints   Vector of <globalDofId,value> pairs
    void giveConstraints( ConstraintVec_ & globalConstraints ) const;

    //! Overwrite values in displacement vector by prescribed values
    void copyConstraintsToIncrements( );
    //@}

    //! @name Actions on tangents (reference and current configuration)
    //@{
protected:
    //! Return non-unit tangent vector
    //!
    //! \param[in]  conf     Spatial configuration
    //! \param[in]  dir      First or second tangent
    //! \param[out] tangent  Tangent vector
    void sumNonUnitTangent_( const enum gshell::fem::config conf,
                             const unsigned dir,
                             VecDim & tangent ) const;

public:
    //! Return unit tangent vector
    //!
    //! \param[in]  dir  First or second tangent
    //! \return          Tangent vector
    VecDim getTangent( const enum gshell::fem::config conf,
                       const unsigned dir ) const;

    //! Return unit normal vector
    //!
    //! \return          Tangent vector
    VecDim getNormal( const enum gshell::fem::config conf ) const;

    //! Return tangent coefficient scaled w.r.t. to unit vector in configuration
    double getTangentCoefficient( const enum gshell::fem::config conf,
                                  NodeStatic_ * nodePtr,
                                  const unsigned dir ) const;

    // ! Return shape gradient
    VecDim getGamma() const { return gamma_; }
    // ! Set shape gradient
    void setGamma(const VecDim & gamma) { gamma_ = gamma; }

    // ! Return area gradient
    VecDim getAreaGrad() const { return aGrad_; }
    // ! Set shape gradient
    void setAreaGrad(const VecDim & grad) { aGrad_ = grad; }
    //@}

protected:

    VecDof          disp_;           //!< converged displacement components
    VecDof          increment_;      //!< current increment of displacements

    VecDof          force_;          //!< nodal force
    VecDof          forceExt_;       //!< nodal external force

    IndexArray_     dofIndices_;     //!< indices of the degrees of freedom

    ConstraintVec_  constraints_;    //!< <dof,value>-pairs for prescribed DOFs

    VecDim          gamma_;          //!< shape gradient
    VecDim          aGrad_;          //!< area gradient

    using BasisNode::tangCoeffs_;
};

#include "NodeStatic.ipp"

#endif
