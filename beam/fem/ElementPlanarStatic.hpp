// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//
//                 Computational Structural Mechanics Lab
//                         University of Cambridge
//
// This software is copyrighted and all rights are retained by the University of
// Cambridge. All use, disclosure, and/or reproduction of any part not expressly
// authorised by the University of Cambridge is prohibited. (C) 2010
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

//! @file   ElementPlanarStatic.hpp

#ifndef beam_fem_elementplanarstatic_h
#define beam_fem_elementplanarstatic_h
//------------------------------------------------------------------------------
#include <cmath>

#include <Eigen/Core>

#include <corlib/eigenX.hpp>
#include <corlib/ElementBasic.hpp>
#include <corlib/verify.hpp>

#include <solid/fem/NodeStatic.hpp>
#include <solid/material/MaterialBase.hpp>


//------------------------------------------------------------------------------
//! Objects related to beams
namespace beam {
    //! Finite element methods for beams
    namespace fem {
        
        template< typename NODE, typename SFUN, typename MAT >
        class ElementPlanarStatic;

        namespace eigenX = corlib::eigenX;

    }
}

//------------------------------------------------------------------------------
/** \defgroup elements Group of Element classes for beams                    */
//------------------------------------------------------------------------------


//------------------------------------------------------------------------------
/** \ingroup elements
 *
 * \brief  Geometrically non-linear planar Bernoulli-beam
 *
 * \details  Provides accessors to nodal quantities, evaluates shape functions
 * and their derivatives, computes kinematics quantities, and evaluates the
 * integral kernels for stiffness and mass matrices and for internal forces.
 *
 * \tparam NODE  Type of nodes this element has as vertices
 * \tparam SFUN  Type of shape function for the geometry and displacements
 * \tparam MAT   Type of material the element has access to
 */
template< typename NODE, typename SFUN, typename MAT >
class beam::fem::ElementPlanarStatic : 
    public corlib::ElementBasic<NODE, corlib::Shapefun<corlib::LINE,2> >
{
public:
    typedef SFUN                                                      ShapeFun;
    typedef NODE                                                      Node;
    typedef MAT                                                       Material;

public:
    //! Topology of element, ie integration cell
    static const corlib::shape myShape       = corlib::LINE;
    //! Number of nodes to read in per element
    static const unsigned     numNodes       = 2;
    //! dimension of embedding space
    static const unsigned     dim            = NODE::dim;
    //! number of degree-of-freedom per node
    static const unsigned     dof            = NODE::dof;
    //! dimension of beam manifold, always 1
    static const unsigned     localDim       = corlib::ShapeTraits< myShape >::dim;
    //! number of vertices of element/cell, always 2
    static const unsigned     numVertices    = corlib::ShapeTraits< myShape >::numVertices;
    //! number of nodes necessary to interpolation, i.e. the "support neighbourhood"
    static const unsigned     numNodesSN     = ShapeFun::numFunctions;
    //! total number of DOFs per element
    static const unsigned     matSize        = numNodesSN * dof;

    //! @name Frequently used types inside and outside this element
    //@{
    typedef eigenX::VectorSd< numNodesSN >               VecNF;
    typedef eigenX::VectorSd< dim >                      VecDim;
    typedef eigenX::VectorSd< dof >                      VecDof;
    typedef eigenX::VectorSd< localDim >                 VecLDim;
    typedef eigenX::MatrixSd< dim, numNodesSN >          MatDimNF;
    typedef eigenX::MatrixSd< dof, numNodesSN >          MatDofNF;
    typedef eigenX::MatrixSd< localDim, numNodesSN >     MatLDimNF;
    typedef Eigen::Matrix< VecNF, localDim, localDim >   MatVecNFLDimLDim;
    typedef eigenX::VectorSd< 3 >                        Vec3;
    typedef eigenX::MatrixSd< 3, 3 >                     Mat3x3;
    //@}

protected:
    typedef Eigen::Matrix< Mat3x3, dim, numNodesSN >                       Mat3x3xDimxNF_;
    typedef corlib::ElementBasic<NODE,corlib::Shapefun<myShape,numNodes> > ElementBasic_;
    typedef eigenX::MatrixSd< matSize, matSize >                           ElemMat_;

protected:
    /** @name Evaluation of shape function and gradient */
    //@{
    //! Return values of shape function at xi
    void sFun_( const VecLDim& xi, VecNF& phi ) const;

    //! Return gradient of shape function at xi
    void sFunGrad_( const VecLDim& xi, MatLDimNF& dphi_dxi ) const;

    //! Return 2nd derivative of shape function at xi
    void sFunGradGrad_( const VecLDim& xi, MatVecNFLDimLDim& dphi_dxi ) const;
    //@}

public:
    /** @name Set-up of element */
    //@{
    //! Constructor
    ElementPlanarStatic() :
        material_( NULL ),
        area_( 0.0 ),
        secMomArea_( 0.0 ),
        energy_( 0.0 ),
        elemStiff_( Eigen::MatrixXd::Zero( matSize, matSize ) ),
        volume_( 0.0 )
    {
        nodes_.assign( NULL );
        FTL_STATIC_ASSERT_MSG( ((numNodesSN-1)%2==1), "Only odd degrees are implemented" );
        FTL_STATIC_ASSERT_MSG( (dof<=dim), "DOF per node must be lower or equal to dim" );
    }

    //! Set specific node
    void setNode( const unsigned index, Node* node )
    {
        FTL_VERIFY_DESCRIPTIVE( (index<numNodesSN), "index is not admissible" );
        nodes_[ index ] = node;
    }

    //! Cache material
    void cacheMaterial( Material * m ) { material_ = m; return; }

    //! Set area of cross section
    void setArea( const double area ) { area_ = area; return; }

    //! Set 2nd moment of area of cross section around z-axis
    void setSecMomArea( const double secMomArea ) { secMomArea_ = secMomArea; return; }

    //! Set thickness of the beam (rectangular cross-section assumed)
    void setThickness ( const double thickness ) 
    {
        setArea( thickness );
        setSecMomArea( std::pow(thickness, 3)/12. );
        return;
    }; 

    //! Membrane and bending elasticity moduli, ie 'EA' and 'EI'
    void elasticityRes( const Mat3x3& baseRef,
                        double& elastModMemb,
                        double& elastModBend ) const;
    
    //@}

    //! Give dof indices of nodes
    void getDofIndices( std::vector<unsigned> & dofIndices ) const;

    //! Evaluate a quantity at given global coordinate. The global co-ordinate is
    //! translated to local co-ordinate by looking at the line formed by
    //! the left and right vertex.
    template< typename OP > 
    typename OP::ReturnType interpolateNodalQuantity( const OP& op,
                                                      const VecDim& coorCur ) const;

    //! Return interpolated reference coordinate at evaluation point
    VecDim giveReferenceCoordinate( const VecLDim& xi ) const;

    //! Return interpolated displacement at evaluation point
    VecDim giveDisplacement( const VecLDim& xi ) const;

    //! Return unit normal to current configuration
    VecDim giveUnitNormalCur( const VecLDim& xi ) const;

    //! Return scaled normal to current configuration whose length carried an
    //! area-scale to the reference configuration
    VecDim giveScaledNormal( const VecLDim& xi ) const;

protected:
    //! Give general nodal quantity accessed by NodeFunctor op
    template< typename OP >
    void nodalQuantity_( MatDofNF& Q, OP op ) const;

public:
    /** @name Accessors for nodal quantity */
    //@{
    //! Build matrix of nodal co-ordinates (ref. configuration) by looping nodes
    void nodalCoordinatesRef( MatDimNF& coorRef ) const;

    //! Build matrix of nodal co-ordinates (current configuration) by looping nodes
    void nodalCoordinatesCur( MatDimNF& coorCur ) const;

    //! Build matrix of nodal co-ordinates by looping nodes
    void nodalDisplacements( MatDofNF& disp ) const;

    //! Return nodal displacement increments
    void nodalIncrements( MatDofNF& dispIncr ) const;
    //@}

protected:
    /** @name Kinematics  */
    //@{
    //!
    //! Determine co-variant reference base vectors at centroid
    //! in global Cartesian components
    void baseRef_( const VecLDim& xi,
                   Mat3x3& base ) const;

    //! Determine co-variant current base vectors at centroid
    //! in global Cartesian components
    void baseCur_( const VecLDim& xi, 
                   Mat3x3& base ) const;

    //! Derivative of reference base with respect to axial/local parameter
    void baseRefGradLoc_( const VecLDim& xi,
                          const Mat3x3& base,
                          Mat3x3& baseGradLoc ) const;

    //! Derivative of current base with respect to axial/local parameter
    void baseCurGradLoc_( const VecLDim & xi,
                          const Mat3x3& base,
                          Mat3x3& BaseGradLoc ) const;

    //! Determine current base vectors gradient
    //! with respect to current configuration
    //! in global Cartesian components
    void baseCurGradCur_( const VecLDim& xi,
                          const Mat3x3& base,
                          Mat3x3xDimxNF_& baseGradient ) const;

    //! Determine current co-variant base vector gradient gradient
    //! with respect to current configuration
    void baseCurGradLocGradCur_( const VecLDim & xi,
                                 const Mat3x3& base,
                                 const Mat3x3& baseGradLoc,
                                 const Mat3x3xDimxNF_& baseGradCur,
                                 Mat3x3xDimxNF_& baseGradLocGradCur ) const;

    //! Compute strain resultants
    void strainRes_( const VecLDim & xi, 
                     const Mat3x3& baseRef,
                     const Mat3x3& baseCur,
                     const Mat3x3& baseRefGradLoc,
                     const Mat3x3& baseCurGradLoc,
                     double& strainMemb,
                     double& strainBend ) const;

    //! Gradient of strain resultants with respect to current
    //! nodal co-ordinates
    void strainResGradCur_( const VecLDim& xi, 
                            const Mat3x3& baseRef,
                            const Mat3x3& baseCur,
                            const Mat3x3& baseCurGradLoc,
                            const Mat3x3xDimxNF_& baseCurGradCur,
                            const Mat3x3xDimxNF_& baseCurGradLocGradCur,
                            MatDimNF& strainMembGradCur,
                            MatDimNF& strainBendGradCur ) const;

    //! Second gradient of strain resultants with respect to current
    //! nodal co-ordinates
    void strainResGradGradCur_( const VecLDim& xi, 
                                const Mat3x3& baseRef,
                                const Mat3x3& baseCur,
                                const Mat3x3& baseCurGradLoc,
                                const Mat3x3xDimxNF_& baseCurGradCur,
                                const Mat3x3xDimxNF_& baseCurGradLocGradCur,
                                ElemMat_& strainMembGradGradCur,
                                ElemMat_& strainBendGradGradCur ) const;

    //! Compute stress resultants
    void stressRes_( const VecLDim & xi, 
                     const Mat3x3& baseRef,
                     const double& strainMemb,
                     const double& strainBend,
                     double& stressMemb,
                     double& stressBend ) const;

    //! Give Jacobian determinant material versus local config
    double jac_( const VecLDim& xi, const bool currentConfig ) const;
    //@}

public:
    /** @name Element strain energy */
    //@{
    //! Clear strain energy
    void clearEnergy() { energy_ = 0.; return; }

    //! Return the element's internal strain energy
    double giveStrainEnergy() const { return energy_; }

    //! Evaluate strain energy density at xi
    void strainEnergyIntegrand( const VecLDim & xi, const double & weight );
    //@}

    //! Evaluate internal force kernel at xi (store results in nodes)
    void internalForceIntegrand( const VecLDim & xi, const double & weight );

    //! Compute nodal forces due to body force (detailed access)
    template< typename FUNC > 
    void bodyForceDetailed( const VecLDim& xi,
                            const double& weight,
                            FUNC func,
                            const double& factor,
                            const bool currentConfig,
                            const bool isAreaIntegrated );

    //! Compute nodal forces due to body force
    template< typename FUNC > 
    void bodyForce( const VecLDim& xi,
                    const double& weight,
                    FUNC func,
                    const double& factor );

    //! compute equivalent nodal forces due to concentacted forces and torques
    void concentratedLoad( const VecLDim& xi,
                           const Vec3& forces,
                           const Vec3& moments,
                           const double& factor );

    //! @name Element stiffness matrix */
    //@{
    //! Clear the element stiffness matrix
    void clearStiffness() { elemStiff_.setZero(); return; }

    //! Add element's stiffness matrix
    void addStiffnessMatrix( Eigen::MatrixXd & matrix ) const { matrix += elemStiff_; }

    //! Evaluate stiffness kernel at xi, weigh result and store it
    void stiffnessIntegrand( const VecLDim & xi, const double & weight );
    //@}

    //! @name Element volume
    //@{
    //! Evaluate infinitesimal volume of element
    //!
    //! Volume is obtained by (Gauss) summation of the infinitesimal volumes
    //! which consist of infinitesimal arc-length times cross section area
    void volumeIntegrand( const VecLDim & xi, const double & weight );

    //! add volume contribution of element to provided variable
    void addVolume( double* volume ) const
    {
        *volume += volume_;
        return;
    }
    //@}

    //! Return stress resultants at evaluation point
    void giveStressResultants( const VecLDim & locCoor,
                               double & stressMemb,
                               double & stressBend ) const;

    //! @name Additional iterator access to the node array stored in this element
    //@{
    typename boost::array<Node*,numNodesSN>::iterator nodesBegin2() { return nodes_.begin(); }
    typename boost::array<Node*,numNodesSN>::iterator nodesEnd2()   { return nodes_.end(); }
    //@}

protected:
    //! Shape function
    ShapeFun      shapeFun_;

    //! @name Material/Geometry
    //@{
    Material*     material_;  //!< Pointer to material object

    double        area_;      //!< beam cross section area
    double        secMomArea_;//!< beam 2nd moment of area of cross section
    //@}

    //! pointers to its nodes
    boost::array<Node*, numNodesSN> nodes_;

    double        energy_;    //!< Element's strain energy
    ElemMat_      elemStiff_; //!< Element stiffness matrix
    double        volume_;    //!< element volume in reference configuration
};


//------------------------------------------------------------------------------
#include "ElementPlanarStatic.ipp"
//------------------------------------------------------------------------------
#endif
