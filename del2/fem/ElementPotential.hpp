// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//
//                 Computational Structural Mechanics Lab
//                         University of Cambridge
//
// This software is copyrighted and all rights are retained by the University of
// Cambridge. All use, disclosure, and/or reproduction of any part not expressly
// authorised by the University of Cambridge is prohibited. (C) 2010
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

//! @file ElementPotential.hpp

#ifndef del2_fem_elementpotential_h
#define del2_fem_elementpotential_h
//------------------------------------------------------------------------------
#include <Eigen/Core>
#include <corlib/eigenX.hpp>

//------------------------------------------------------------------------------
namespace del2
{
    namespace fem
    {
        template<typename BELEMENT> class ElementPotential;
        namespace eigenX = corlib::eigenX;
    }
}

//------------------------------------------------------------------------------
/** \brief Element for the solution of Poisson/Laplace equation 
 *  \details This element addresses the scalar (potential) equation 
 *  \f[ 
 *          - \mu \nabla^2 u = f
 *  \f] 
 *  with the conductivity \f$ \mu \f$ and provides the stiffness matrix
 *  for a standard finite element discretization thereof.  Additionally, the
 *  dof handling and some field evaluations are provided.
 *  \tparam FELEMENT  Type of basis element 
 */
template<typename BELEMENT>
class del2::fem::ElementPotential: public BELEMENT
{
public:

    //----------------------------------------------------------------------
    //! @name Basic typedefs
    //@{ 
    typedef          BELEMENT                BasisElement;
    typedef typename BasisElement::Node      Node;
    //@}

    //----------------------------------------------------------------------
    //! @name Basic attributes
    //@{
    static const unsigned      dim            = Node::dim;
    static const unsigned      dof            = Node::dof;
    static const unsigned      localDim       = BasisElement::localDim;
    static const unsigned      numNodes       = BasisElement::numNodes;
    FTL_STATIC_ASSERT_MSG((dof==1), "Implementation requires dof=1");
    //@}

    //----------------------------------------------------------------------
    //! @name Convenience typedefs 
    //@{
    typedef typename BasisElement::VecDim              VecDim;
    typedef typename BasisElement::VecLDim             VecLDim;
    typedef typename BasisElement::VecNN               VecNN;
    typedef typename BasisElement::MatDimNN            MatDimNN;
    typedef eigenX::MatrixSd<dof, numNodes>            MatDofNN;
    typedef eigenX::VectorSd<dof>                      VecDof;
    typedef eigenX::VectorSd<3>                        Vec3;
    //@

protected:
    //! used for external evaluation functors
    typedef del2::fem::ElementPotential<BELEMENT>  SelfType_;

public:
    //! Give dof indices of nodes
    void getDofIndices( std::vector<unsigned> & dofIndices ) const;

    //! Give constant of Elasticity
    double getConductivity() const { return conductivity_; }

    //! Evaluate stiffness kernel at xi, weigh result and store it
    void stiffnessIntegrand( const VecLDim & xi, const double & weight,
                             Eigen::MatrixXd & result ) const;

    //! Evaluate the co-normal derivative
    void nitscheTraction( const VecLDim & xi, const VecDim & normal,
                          eigenX::MatrixSd<dof,numNodes*dof> & result ) const;

    //! Set value of the conductivity
    void setConductivity( const double & conductivity )
    {
        conductivity_ = conductivity;
    }

    //! Compute the potential at any interior point
    eigenX::VectorSd<Node::dof>
    potential( const VecLDim & xi ) const;

    //! Compute and return internal flux at xi
    Vec3 internalFlux( const VecLDim & xi ) const;

    //! Adapter for norm-error
    eigenX::MatrixSd<Node::dim,Node::dof>
    potentialGradient( const VecLDim & xi ) const
    {
        Vec3 flux = this -> internalFlux( xi );
        const VecDim flux2 = flux.head(dim);
        eigenX::VectorSd<Node::dim> result;
        result.col(0) = flux2;
        return result;
    }

    //! Body force function
    template<typename FUNC>
    void bodyForce( const VecLDim & xi, const double & weight,
                    FUNC func, const double & factor );

    //! Equation residual
    template<typename FUNC>
    typename Node::VecDof equationResidual( const VecLDim & xi, FUNC func ) const;

    //! Give element stiffness derivative with respect to conductivity
    void stiffnessDerivIntegrand( const VecLDim & xi, const double & weight,
                                  Eigen::MatrixXd & result ) const;

protected:
    double    conductivity_; //!< Conductivity
};
//------------------------------------------------------------------------------
#include "ElementPotential.ipp"
//------------------------------------------------------------------------------
#endif
