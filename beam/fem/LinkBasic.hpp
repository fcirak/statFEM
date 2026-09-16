// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//
//                 Computational Structural Mechanics Lab
//                         University of Cambridge
//
// This software is copyrighted and all rights are retained by the University of
// Cambridge. All use, disclosure, and/or reproduction of any part not expressly
// authorised by the University of Cambridge is prohibited. (C) 2010
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

//! @file LinkBasic.hpp

#ifndef beam_fem_linkbasic_h
#define beam_fem_linkbasic_h
//------------------------------------------------------------------------------
#include <iostream>
#include <cassert>
#include <algorithm>
#include <Eigen/Core>
#include <corlib/eigenX.hpp>

namespace beam {
    namespace fem {

        template< typename NODE >
        class LinkBasic;

        namespace eigenX = corlib::eigenX;
    }
}

//------------------------------------------------------------------------------
/** \brief Basic link template
 *
 *  \details A link object establishes a linear combination
 *  of several DOFs \f$u_{j(i)}^{k(i)}\f$ (attached to nodes \f$k\f$) subject to factors
 *  \f$\alpha_i\f$ and a right hand side \f$v\f$ (the #value_):
 *  \f[
 *       v = \sum_i \big( \alpha_i u^{k(i)}_{j(i)} \big)
 *  \f]
 *
 *  The link is enforced implicitly
 *  with a Lagrange multiplier whose DOF index is stored on #dofIndex_.
 *  The respective coupling matrices and associated DOF indices
 *  can be retrieved. The link class essentially behaves like
 *  an element, but it does not require a quadrature.
 *
 *  \tparam NODE   Type of node as vertex
 */
template< typename NODE >
class beam::fem::LinkBasic
{
public:
    typedef NODE                                              Node;

    //! Number of primary DOFs: the one Lagrange multiplier
    static const unsigned dof = 1;

    //! Vector of primary DOFs: the one Lagrange multiplier
    typedef eigenX::VectorSd< dof >              VecDof;

protected:
    //! DOF index array (it's of length 1)
    typedef std::array< unsigned, dof >                       IndexArray_;
    //! The tuple storing a summand in the link connection.
    //! It is made up of: pointer to node, the (local) node DOF, and the factor/coefficient
    typedef std::tuple< const Node *, unsigned, double >      NodeDofCoeff_;
    //! List/Vector of all summands of link
    typedef std::vector< NodeDofCoeff_ >                      NodeDofCoeffVec_;
    typedef typename NodeDofCoeffVec_::iterator               NodeDofCoeffVecIt_;
    typedef typename NodeDofCoeffVec_::const_iterator         NodeDofCoeffVecConstIt_;

public:
    //! Default constructor which does nothing
    LinkBasic()
        : dofIndex_(),
          value_( Eigen::VectorXd::Zero( dof ) ),
          multiplier_( Eigen::VectorXd::Zero( dof ) ),
          increment_( Eigen::VectorXd::Zero( dof ) )
    {
        dofIndex_.fill( 0 );
    }

    //! Set the vertex pointer index to node
    //!
    //! \param[in]    node     pointer to node \f$k(i)\f$
    //! \param[in]    dof      local node DOF \f$j(i)\f$
    //! \param[in]    coeff    factor \f$\alpha_i\f$
    void addNodeDofCoeff( const Node * node, const unsigned dof, const double coeff )
    {
        nodeDofCoeffs_.push_back( std::make_tuple( node, dof, coeff ) );
        return;
    }

    //! Set the i-th (node, dof, coeff) tuple of this link
    //!
    //! \param[in]    i        tuple index 
    //! \param[in]    node     pointer to node \f$k(i)\f$
    //! \param[in]    dof      local node DOF \f$j(i)\f$
    //! \param[in]    coeff    factor \f$\alpha_i\f$
    void setNodeDofCoeff( const unsigned i,
                          const Node * node, const unsigned dof, const double coeff )
    {
        nodeDofCoeffs_[ i ] = std::make_tuple( node, dof, coeff );
        return;
    }

    //! Get the i-th (node, dof, coeff) tuple of this link
    //!
    //! \param[in]     i        tuple index 
    //! \param[in,out] node     pointer to node \f$k(i)\f$
    //! \param[in,out] dof      local node DOF \f$j(i)\f$
    //! \param[in,out] coeff    factor \f$\alpha_i\f$
    void getNodeDofCoeff( const unsigned i,
                          const Node * & node, unsigned & dof, double & coeff ) const
    {
        node = std::get< 0 >( nodeDofCoeffs_[ i ] );
        dof = std::get< 1 >( nodeDofCoeffs_[ i ] );
        coeff = std::get< 2 >( nodeDofCoeffs_[ i ] );
        return;
    }


    //! Set RHS of link
    //!
    //! \param[in]    value    right-hand-side value \f$v\f$
    void setValue( const VecDof value )
    {
        value_ = value;
        return;
    }

    //! Get RHS of link
    VecDof getValue( ) const
    {
        return value_;
    }


    //! Virtual destructor to please the compiler
    virtual ~LinkBasic( )
    { 
        nodeDofCoeffs_.clear();
        return;
    }

    //! Assign numbers to the degrees of freedom
    void numberDOFs( unsigned & counter )
    {
        dofIndex_[ 0 ] = counter++;
        return;
    }

    //! Give global DOF indices of linked DOF of stored nodes
    void getDofIndices( std::vector<unsigned> & dofIndices ) const;

    //! Give global DOF index of Lagrange multiplier
    void getDofIndicesP( std::vector<unsigned> & dofIndex ) const;

    //! Set Lagrange multiplier
    void setMultiplier( const VecDof & mult ) { multiplier_ = mult; }

    //! Return the associated Lagrange multiplier of link
    VecDof giveMultiplier( ) const { return ( multiplier_ + increment_ ); }

    //! Clear Lagrange multiplier increment
    void clearIncrement( ) { increment_.setZero( ); }

    //! Return Lagrange multiplier increment
    VecDof getIncrement() const { return increment_; }

    //! Set Lagrange multiplier increment
    void setIncrement( const VecDof & inc ) { increment_ = inc; }

    //! Update Lagrange multiplier increment
    void addToIncrement( const VecDof & inc ) { increment_ += inc; }

    //! Update Lagrange multiplier and clear its increment
    void updateMultiplier( )
    {
        multiplier_ += increment_;
        increment_.setZero( );
        return;
    }

    //! Return link matrix containing DOF coupling
    void giveCouplingMatrix( Eigen::MatrixXd & result ) const;

    //! Return transposed link matrix containing DOF coupling
    void giveCouplingMatrixT( Eigen::MatrixXd & result ) const
    { 
        Eigen::MatrixXd tmp( result.cols(), result.rows() );
        this -> giveCouplingMatrix( tmp );
        result = tmp.transpose( );
    }

    //! Return negative RHS contribution onto main field
    Eigen::VectorXd giveMultiplierRhs( ) const;

    //! Return negative residuum of link
    Eigen::VectorXd giveLinkResiduum( ) const;

    //! Return negative residuum of link with scaled prescribed value
    Eigen::VectorXd giveLinkResiduumScaled( const double factor ) const;

    //! Return RHS of link
    VecDof giveValue() const
    {
        return value_;
    }

    //! Print contents
    std::ostream & print( std::ostream & os ) const;

    //! Print connectivity 
    std::ostream & printDofs( std::ostream & os ) const;

protected:
    NodeDofCoeffVec_    nodeDofCoeffs_;    //!< pointers to linked nodes
    IndexArray_         dofIndex_;         //!< index of Lagrange multiplier
    VecDof              value_;            //!< value of link (RHS)
    VecDof              multiplier_;       //!< Lagrange multiplier
    VecDof              increment_;        //!< Lagrange multiplier increment
};


//------------------------------------------------------------------------------
#include "LinkBasic.ipp"


#endif
