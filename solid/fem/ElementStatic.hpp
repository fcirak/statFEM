// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//
//                 Computational Structural Mechanics Lab
//                         University of Cambridge
//
// This software is copyrighted and all rights are retained by the University of
// Cambridge. All use, disclosure, and/or reproduction of any part not expressly
// authorised by the University of Cambridge is prohibited. (C) 2010
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

//! @file ElementStatic.hpp

#ifndef solid_fem_elementstatic_h
#define solid_fem_elementstatic_h
//------------------------------------------------------------------------------
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <solid/fem/NodeStatic.hpp>

//------------------------------------------------------------------------------
namespace solid{
    namespace fem{
        
        template<typename BELEMENT,typename MAT> class ElementStatic;
        namespace ublas = boost::numeric::ublas;
    }
}

//------------------------------------------------------------------------------
/** \brief Element for standard solids, extends corlib::ElementBasic
 * \details Provides accessors to nodal quantities, evaluates shape functions
 * and their derivatives, computes kinematics quantities, and evaluates the
 * integral kernels for stiffness and mass matrices and for internal forces.
 * \tparam BELEMENT Type of basis element
 * \tparam MAT      Type of material the element has access to
 */
template<typename BELEMENT, typename MAT>
class solid::fem::ElementStatic : public BELEMENT
{
public:
    //! @name Basic typedefs
    //@{ 
    typedef          BELEMENT                BasisElement;
    typedef          MAT                     Material;
    typedef typename BasisElement::Node      Node;
    //@}

    //! @name Basic attributes
    //@{
    static const unsigned      dim            = Node::dim;
    static const unsigned      dof            = Node::dof;
    static const unsigned      localDim       = BasisElement::localDim;
    static const unsigned      numNodes       = BasisElement::numNodes;

    //! @name Convenience typedefs 
    //@{
    typedef typename BasisElement::VecDim              VecDim;
    typedef typename BasisElement::VecLDim             VecLDim;
    typedef typename BasisElement::VecNN               VecNN;
    typedef typename BasisElement::MatDimNN            MatDimNN;
    typedef ublas::bounded_matrix<double,dof,numNodes> MatDofNN;
    typedef ublas::bounded_vector<double,dof>          VecDof;
    typedef ublas::bounded_matrix<double, 3, 3>        Mat3x3;
    //! Conformity typedef for corlib::NewmarkAssembler
    typedef          MatDofNN                          MatDofNF;
    //@

protected:
    typedef ublas::bounded_matrix<double, dim, dim>                 MatDimDim_;
    typedef ublas::bounded_vector<MatDimDim_, dim>                  MatDimDimDim_;


protected:
    //! Give general nodal quantity accessed by NodeFunctor op
    template<typename OP> void nodalQuantity_( MatDofNN & Q, OP op ) const;

public:
    /** @name Accessors for nodal quantity */
    //@{
    //! Return nodal displacements
    void nodalDisplacements( MatDofNN & U ) const;

    //! Return nodal increments
    void nodalIncrements( MatDofNN & deltaU ) const;

    void nodalCoordinates( MatDimNN & X ) const;

    //! Give dof indices of nodes
    void getDofIndices( std::vector<unsigned> & dofIndices ) const;
    //@}

    //! Give constant of Elasticity
    double giveYoungsModulus() const { return material_->estimatedElasticityConst(); }

protected:
    //! Compute the full set of kinematic quantities: detG, F, DphiDX
    double kinematics_( const VecLDim & xi, MatDimNN & dphi_dX, Mat3x3   & F ) const;

    //! Compute global derivative of deformation gradient \f$F_{iJ,K}\f$
    double deformationGradientGlobalDerivatives_( const VecLDim & xi,
                                                  MatDimDimDim_ & fGrad ) const;


public:
    /** @name Element stiffness matrix  */
    //@{
    //! Evaluate stiffness kernel at xi, weigh result and store it
    void stiffnessIntegrand( const VecLDim & xi, const double & weight,
                             ublas::matrix<double> & result ) const;
    //@}

    /** @name Element strain energy */
    //@{
    //! Clear strain energy
    void clearEnergy() { energy_ = 0.; return; }

    //! Return the element's internal strain energy
    double giveStrainEnergy( ) const { return energy_; }

    //! Evaluate strain energy density at xi
    void strainEnergyIntegrand( const VecLDim & xi, const double & weight );
    //@}

    //! Evaluate internal force kernel at xi (store results in nodes)
    void internalForceIntegrand( const VecLDim & xi, const double & weight );

    //! Evaluate the traction term used for the Nitsche technique at xi
    void nitscheTraction( const VecLDim & xi, const VecDim & normal,
    		              ublas::bounded_matrix<double,dof,
    		              numNodes*dof> & result ) const;

    //! Compute nodal forces due to body force
    template< typename FUNC > 
    void bodyForce( const VecLDim & xi, const double & weight, FUNC func,
                    const double & factor );

    //! Compute and return first Piola-Kirchhoff stress tensor at xi
    Mat3x3 firstPiolaKirchhoff( const VecLDim & xi ) const;

    //! Compute and return second Piola-Kirchhoff stress tensor at xi
    Mat3x3 secondPiolaKirchhoff( const VecLDim & xi ) const;

    //! Compute and return Cauchy stress tensor at xi
    Mat3x3 cauchyStress( const VecLDim & xi ) const;

    //! Cache material
    void cacheMaterial( MAT * m ) { material_ = m; return; }

    //! Evaluate stress residual at given parametric point
    template<typename FUNC>
    void stressResidual( const VecLDim & xi, FUNC func, VecDim & stressRes );

protected:
    MAT*      material_;  //!< Pointer to material object
    double    energy_;    //!< Element's strain energy
};


//------------------------------------------------------------------------------
#include "ElementStatic.ipp"
//------------------------------------------------------------------------------
#endif
