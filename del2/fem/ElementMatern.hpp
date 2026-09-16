// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//
//                 Computational Structural Mechanics Lab
//                         University of Cambridge
//
// This software is copyrighted and all rights are retained by the University of
// Cambridge. All use, disclosure, and/or reproduction of any part not expressly
// authorised by the University of Cambridge is prohibited. (C) 2010 --
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

//! @file ElementMatern.hpp

#ifndef del2_fem_elementmatern_h
#define del2_fem_elementmatern_h
//------------------------------------------------------------------------------
#include <Eigen/Core>
#include <corlib/eigenX.hpp>

namespace eigenX = corlib::eigenX;

//------------------------------------------------------------------------------
namespace del2
{
    namespace fem
    {
        template<typename BELEMENT> class ElementMatern;
        namespace eigenX = corlib::eigenX;
    }
}

//------------------------------------------------------------------------------
/** \brief Element for the solution of stochastic partial differential equation
 *  \details
 *
 *  - This element addresses the stochastic partial differential equation
 *    \f[
 *          ( \kappa^2 u - \nabla \cdot ( H \nabla u ) )^\beta = f
 *    \f]
 *
 *  - Provides finite element matrices for the case \f$ \beta = 1 \f$.
 *
 *  - The diffusivity matrix \f$ H \in \mathbb{R}^{d \times d} \f$ must be
 *    symmetric positive definite.
 *
 *  - The constant coefficient \f$ \kappa \f$ is related to the smoothness
 *    \f$ \nu \f$ and length scale \f$ \ell \f$ parameters of the Matern
 *    kernel such that \f$  \kappa = \frac{ \sqrt{ 2 \nu } }{ \ell } \f$.
 *
 *  \tparam BELEMENT  Type of basis element
 */
template<typename BELEMENT>
class del2::fem::ElementMatern: public BELEMENT
{
public:

    //----------------------------------------------------------------------
    //! @name Basic typedefs
    ///@{
    typedef          BELEMENT                BasisElement;
    typedef typename BasisElement::Node      Node;
    ///@}

    //----------------------------------------------------------------------
    //! @name Basic attributes
    ///@{
    static const unsigned dim      = Node::dim;
    static const unsigned dof      = Node::dof;
    static const unsigned localDim = BasisElement::localDim;
    static const unsigned numNodes = BasisElement::numNodes;
    FTL_STATIC_ASSERT_MSG((dof==1), "Implementation requires dof=1");
    ///@}

    //----------------------------------------------------------------------
    //! @name Convenience typedefs 
    ///@{
    typedef typename BasisElement::VecLDim             VecLDim;
    typedef typename BasisElement::VecNN               VecNN;
    typedef typename BasisElement::MatDimNN            MatDimNN;
    typedef eigenX::MatrixSd<dof, numNodes>            MatDofNN;
    typedef eigenX::MatrixSd<dim>                      MatDimDim;
    ///@}

public:

    //! Give dof indices of nodes
    void getDofIndices ( std::vector<unsigned> & dofIndices ) const;

    //! Give kappa
    double getKappa() const { return kappa_; }

    //! Give diffusivity
    void getDiffusivity( MatDimDim & result ) const { result = diffusivity_; }

    //! Evaluate Matern integrand at xi, weigh result and store it
    void maternIntegrand ( const VecLDim & xi, const double & weight,
                           Eigen::MatrixXd & result ) const;

    //! Evaluate unit lumped mass integrand at xi, weigh result and pass to nodes
    void lumpedMassIntegrand ( const VecLDim & xi, const double & weight );

    //! Set value of kappa
    void setKappa ( const double & kappa ) { kappa_ = kappa; }

    //! Set value of diffusivity
    void setDiffusivity ( const MatDimDim & diffusivity ) {

        diffusivity_ = diffusivity;
        return;
    }

private:
    double      kappa_;       //!< Coefficient of u
    MatDimDim   diffusivity_; //!< Diffusivity

};
//------------------------------------------------------------------------------
#include "ElementMatern.ipp"
//------------------------------------------------------------------------------
#endif
