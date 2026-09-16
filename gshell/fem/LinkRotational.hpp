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

#ifndef gshell_fem_linkrotational_h
#define gshell_fem_linkrotational_h

//------------------------------------------------------------------------------
#include <iostream>
#include <cassert>
#include <array>

#include <gshell/fem/LinkPenaltyLMultiplier.hpp>
#include <gshell/fem/Config.hpp>

namespace gshell {
    namespace fem {

        template< typename ELEMENT >
        class LinkRotationalExt;

        template< typename ELEMENT >
        class LinkRotationalInt;

        namespace eigenX = corlib::eigenX;
    }
}

//==============================================================================
/// Rotational link constraining a normal to a external (prescribed) vector.
/// The normal is evaluated at a node on the surface.
///
/// \tparam  ELEMENT  Element type, which describes the surface and has a tangent
template< typename ELEMENT >
class gshell::fem::LinkRotationalExt :
    public gshell::fem::LinkPenaltyLMultiplier< typename ELEMENT::Node, 2 >
{
public:
    typedef ELEMENT                                           Element;
    typedef typename Element::Node                            Node;

    /// Constraint mode
    enum mode {
        PENALTY,   //!< penalty method
        LAGMULT,   //!< Lagrange multiplier method
        AUGLAG     //!< Lagrange multipler method + penalty method
    };

private:
    /// Base class type
    typedef LinkPenaltyLMultiplier< Node, 2 >                 LinkPenaltyLMultiplier_;
    /// Its own type
    typedef LinkRotationalExt< Element >                      LinkRotationalExt_;

public:
    /// Number of primary DOFs: the one Lagrange multiplier
    static const unsigned dof = LinkPenaltyLMultiplier_::dof;
    
    /// Vector of primary DOFs: the one Lagrange multiplier
    typedef eigenX::VectorSd<dof>                     VecDof;

protected:
    /// Spatial vector
    typedef typename Node::VecDim                     VecDim_;
    /// General matrix
    typedef Eigen::MatrixXd                           Mat_;
    /// General vector
    typedef Eigen::VectorXd                           Vec_;

public:
    /// Default constructor
    ///
    /// \param[in]  node       Pointer to node whose tangent is constraint
    /// \param[in]  elememt    Element of surface at node
    LinkRotationalExt( Node * node, Element * element );

    /// Virtual destructor to please the compiler
    virtual ~LinkRotationalExt( ) { }

    /// Return pointer to normal's node
    Node * giveNode( ) const { return node_; }

    /// Compute constraint residuals
    void computeResidual( );
    /// Compute Jacobian matrix of constraint residula w.r.t. DOFs
    void computeMatrix( );
    /// Compute constraint residuals and Jacobian matrix
    void computeResidualAndMatrix( );

    /// Print contents
    ///
    /// \param[in,out]   os    Output stream
    /// \return                Output stream
    std::ostream & write( std::ostream & os ) const;

private:
    /// Compute local index of node (vertex) w.r.t. element
    static unsigned nodeLocalIndex_( Element * element,
                                     Node * node )
    {
        unsigned i;
        for ( i = 0; i < Element::numNodes; ++i )
            if ( element->giveNodePtr( i ) == node )
                break;
        FTL_VERIFY( i < Element::numNodes );
        return i;
    }

protected:
    using LinkPenaltyLMultiplier_::nodes_; 
    using LinkPenaltyLMultiplier_::penalty_;          ///< penalty factor
    using LinkPenaltyLMultiplier_::multiplier_;       ///< Lagrange multiplier

    using LinkPenaltyLMultiplier_::resDisp_;          ///< Residual at nodal DOFs (displacements)
    using LinkPenaltyLMultiplier_::resMult_;          ///< Residual at Lagrange multiplier
    using LinkPenaltyLMultiplier_::matDispDisp_;      ///< derivative of #resDisp_ w.r.t. nodal DOFs (displacements)
    using LinkPenaltyLMultiplier_::matMultDisp_;      ///< derivative of #resMult_ w.r.t. nodal DOFs (displacements)

    Node *      node_;             ///< pointer to node
    Element *   element_;          ///< pointer to element

    VecDim_     normal_;           ///< reference normal
};


//==============================================================================
/// Rotational link constraining two normals of two adjacent surfaces.
/// The normals are evaluated at a node on the surfaces.
///
/// \tparam  ELEMENT  Element type, which describes the surfaces and tangents
template< typename ELEMENT >
class gshell::fem::LinkRotationalInt :
    public gshell::fem::LinkPenaltyLMultiplier< typename ELEMENT::Node, 1 >
{
public:
    typedef ELEMENT                                           Element;
    typedef typename Element::Node                            Node;

private:
    /// Base class type
    typedef LinkPenaltyLMultiplier< Node, 1 >                 LinkPenaltyLMultiplier_;
    /// Its own type
    typedef LinkRotationalInt< Element >                      LinkRotationalInt_;

public:
    /// Number of primary DOFs: the one Lagrange multiplier
    static const unsigned dof = LinkPenaltyLMultiplier_::dof;
    
    /// Vector of primary DOFs: the one Lagrange multiplier
    typedef eigenX::VectorSd<dof>              VecDof;

protected:
    /// Spatial vector
    typedef typename Node::VecDim              VecDim_;
    /// General matrix
    typedef Eigen::MatrixXd                    Mat_;
    /// General vector
    typedef Eigen::VectorXd                    Vec_;

public:
    /// Default constructor
    ///
    /// \param[in]  node       Pointer to node whose tangents are constraint
    /// \param[in]  elememt0   Element of 1st surface at node
    /// \param[in]  elememt1   Element of 2nd surface at node
    LinkRotationalInt( Node * node, Element * element0, Element * element1 );

    /// Virtual destructor to please the compiler
    virtual ~LinkRotationalInt( )
    { 
        return;
    }

    /// Return pointer to normal's node
    Node * giveNode( ) const { return node_; }

    /// Set element
    void setElement( const unsigned elemLocId, Element * elem )
    {
        if ( elemLocId == 0 ) element0_ = elem;
        else                  element1_ = elem;
        this->computeNormal_( elemLocId );
        return;
    }

    /// Retrieve an element
    Element * getElement( const unsigned elemLocId ) const
    {
        if ( elemLocId == 0 ) return element0_;
        else                  return element1_;
    }

    /// Compute constraint residuals
    void computeResidual( );
    /// Compute Jacobian matrix of constraint residula w.r.t. DOFs
    void computeMatrix( );
    /// Compute constraint residuals and Jacobian matrix
    void computeResidualAndMatrix( );

    /// Print contents
    ///
    /// \param[in,out]   os    Output stream
    /// \return                Output stream
    std::ostream & write( std::ostream & os ) const;

protected:
    /// (Re-)Compute 1st or 2nd normal
    void computeNormal_( const unsigned elemLocId );

private:
    /// Compute local index of node (vertex) w.r.t. element
    static unsigned nodeLocalIndex_( Element * element,
                                     Node * node )
    {
        unsigned i;
        for ( i = 0; i < Element::numNodes; ++i )
            if ( element->giveNodePtr( i ) == node )
                break;
        FTL_VERIFY( i < Element::numNodes );
        return i;
    }

protected:
    using LinkPenaltyLMultiplier_::nodes_; 
    using LinkPenaltyLMultiplier_::penalty_;          ///< penalty factor
    using LinkPenaltyLMultiplier_::multiplier_;       ///< Lagrange multiplier

    using LinkPenaltyLMultiplier_::resDisp_;          ///< Residual at nodal DOFs (displacements)
    using LinkPenaltyLMultiplier_::resMult_;          ///< Residual at Lagrange multiplier
    using LinkPenaltyLMultiplier_::matDispDisp_;      ///< derivative of #resDisp_ w.r.t. nodal DOFs (displacements)
    using LinkPenaltyLMultiplier_::matMultDisp_;      ///< derivative of #resMult_ w.r.t. nodal DOFs (displacements)

    Node *      node_;             ///< pointer to node common to both elements
    Element *   element0_;         ///< pointer to 1st element
    Element *   element1_;         ///< pointer to 2nd element

    VecDim_     normal0_;          ///< reference normal of 1st surface
    VecDim_     normal1_;          ///< reference normal of 2nd surface
};


//------------------------------------------------------------------------------
#include "LinkRotational.ipp"


#endif
