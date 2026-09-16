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
//! @date   2011

#ifndef gshell_fem_elementklstatic_h
#define gshell_fem_elementklstatic_h

//------------------------------------------------------------------------------
#include <functional>
#include <array>

#include <Eigen/Core>

#include <corlib/eigenX.hpp>
#include <corlib/Accessor.hpp>
#include <corlib/ElementBasic.hpp>

#include <gshell/fem/linalg.hpp>
#include <gshell/fem/NodeStatic.hpp>
#include <gshell/fem/SurfaceGeometer.hpp>
#include <gshell/fem/Config.hpp>


//------------------------------------------------------------------------------
namespace gshell {
    namespace fem {
        
        template< typename BELEMENT, typename MAT >
        class ElementKLStatic;

        namespace eigenX = corlib::eigenX;
    }
}

//------------------------------------------------------------------------------
/** \defgroup elements Group of Element classes for g-shells                  */
//------------------------------------------------------------------------------


//------------------------------------------------------------------------------
/** \ingroup elements
 * \brief Element for Kirchhoff-Love shell
 *        extends corlib::ElementBasic
 *
 * \details Provides accessors to nodal quantities, evaluates shape functions
 * and their derivatives, computes kinematics quantities, and evaluates the
 * integral kernels for stiffness and mass matrices and for internal forces.
 *
 * About indices in formulas:
 * - \f$A,B,C,D,...=1,2,3\f$ are used for the global Cartesian basis
 *                         \f$\{\vec{e}_A\}\f$
 *                         and its co-ordinates \f$(x^A)\f$
 * - \f$i,j,k,...=1,2,3\f$ are used for the local curvi-linear bases 
 *                         \f$\{\vec{g}_i\}\f$ (co-variant current),
 *                         \f$\{\vec{g}^i\}\f$ (contra-variant current),
 *                         \f$\{\bar{\vec{g}}_i\}\f$ (co-variant reference),
 *                         \f$\{\bar{\vec{g}}^i\}\f$ (contra-variant refence)
 *                         and its co-ordinates \f$(\xi^i)\f$
 * - \f$\alpha,\beta,...=1,2\f$ are used for the local curvi-linear mid-plane bases 
 *                         \f$\{\vec{a}_\alpha\}\f$ (co-variant current),
 *                         \f$\{\vec{a}^\alpha\}\f$ (contra-variant current),
 *                         \f$\{\bar{\vec{a}}_\alpha\}\f$ (co-variant reference),
 *                         \f$\{\bar{\vec{a}}^\alpha\}\f$ (contra-variant refence),
 *                         and its co-ordinates \f$(\xi^\alpha)\f$
 * - \f$\mathcal{A},\mathcal{B},...=1,2,3\f$ are used for symmetric Voigt-indices of 
 *                         \f$(\alpha,\beta)\f$-indexed quantities, see also #voigtForward
 * - \f$K,L,...=1,...,NF\f$ are the local node indices attached to the
 *                         shape functions \f$N^K(\xi^\alpha)\f$ whose support
 *                         intersects with the element domain
 *
 * \tparam  BELEMENT   Basic element for geometry-related and shape-function methods
 * \tparam  MAT        Type of material the element has access to
 */
template< typename BELEMENT, typename MAT >
class gshell::fem::ElementKLStatic : public BELEMENT
{
public:
    typedef BELEMENT                            BasisElement;
    typedef MAT                                 Material;
    //! Type of node carrying the degrees of freedom
    typedef typename BasisElement::Node         Node;
    //! Type of shape function for the geometry and displacements
    typedef typename BasisElement::ShapeFun     ShapeFun;

public:
    //! Enumerator for description of coordinate system
    enum coordSys {
        LOCAL,           //!< local coordinates
        GLOBAL           //!< global coordinates
    };

    //! In-/Out-of-plane mode
    enum plane {
        INPLANE,         //!< in-plane, i.e. membrane
        OUTOFPLANE       //!< out-of-plane, i.e. bending
    };

    //! Strain measures
    enum strains {
        LOCALGL,         //!< Green-Lagrange in local co-ordinates
        GREENLAGRANGE,   //!< Green-Lagrange in global co-ordinates
        EULERALMANSI     //!< Euler-Almansi in global co-ordinates
    };

    //! Stress measures
    enum stresses {
        LOCALPK2,        //!< 2nd Piola-Kirchhoff in local co-ordinates
        PIOLAKIRCHHOFF2, //!< 2nd Piola-Kirchhoff in global co-ordinates
        CAUCHY           //!< Cauchy in global co-ordinates
    };

    //! Simplex shape of element -- can be corlib::TRIANGLE or corlib::QUADRILATERAL
    static const corlib::shape myShape   = BasisElement::myShape;
    //! Number of vertices of simplex shape -- can be 3 or 4
    static const unsigned      numVertices = BasisElement::numVertices;  // MIGHT BE CALLED numNodes ???
    //! Dimension of embedding space -- always 3
    static const unsigned      dim       = Node::dim;  // always 3
    //! Total number of DOFs at nodes -- always 3
    static const unsigned      dof       = Node::dof;  // always 3
    //! Dimension of shell manifold -- always 2
    static const unsigned      localDim  = BasisElement::localDim;
    //! Number of effective second derivatives subtracting duplicities due to symmetry -- always 3
    static const unsigned      sDim      = BasisElement::sDim;  // always 3
    //! Map to get Voigt indices from pair of local indices
    //! \f$(\alpha,\beta) \mapsto \mathcal{A}\f$
    static const unsigned      voigtForward[localDim][localDim];// = { { 0, 1 }, { 1, 2 } };

    typedef Eigen::VectorXd                  VecNF;
    typedef eigenX::VectorSd< dim >          VecDim;
    typedef eigenX::VectorSd< dof >          VecDof;
    typedef eigenX::VectorSd< localDim >     VecLDim;
    typedef eigenX::MatrixSd< 3, 3 >         Mat3x3;
    typedef Eigen::MatrixXd                  MatDimNF;
    typedef Eigen::MatrixXd                  MatDofNF;
    typedef Eigen::MatrixXd                  MatLDimNF;
    typedef Eigen::MatrixXd                  MatSDimNF;
    typedef Eigen::MatrixXd                  Mat;

private:
    typedef ElementKLStatic< BELEMENT, MAT >                       ElementKLStatic_;

public:
    //! Constructor
    //!
    //! Note:
    //!    - Pointer to element's shape function has to be set
    //!    - The BasisElement::suppportNodes (nodes carrying DOFs) have to assigned
    ElementKLStatic() :
        BasisElement(),
        material_( NULL ),
        thickness_( 0.0 )
    {
        static_assert( dim == 3 );
        static_assert( dof == 3 );
        static_assert( localDim == 2 );
        return;
    }

    //! Destructor
    virtual ~ElementKLStatic() { return; }

    //! Set thickness
    //!
    //! \param[in]  thickness   Provided thickness value, constant on element
    void setThickness( const double thickness ) { thickness_ = thickness; return; }

    //! Return the indices of the nodal degrees of freedom
    //!
    //! \param[in]  dofIndices Storage vector which some DOFs can be written.
    //!                        The vector must come _allocated_ to required size.
    void getDofIndices( std::vector< unsigned > & dofIndices ) const;

protected:
    //! Collect general nodal quantitiy accessed by the passed operator op
    //!
    //! \tparam     OP     Node functor to access the desired nodal quantity
    //!
    //! \param[out] quanU  Matrix (dim x number of nodes) containing the result
    //! \param[in]  op     Node functor for access of nodal quantity
    template< typename OP >
    void nodalQuantity_( MatDimNF & quanU, OP op ) const;

    //! Update general nodal quantitiy accessed by the passed operator op
    //!
    //! \tparam     OP     Node functor to access the desired nodal quantity
    //!
    //! \param[out] quanU  Matrix (dim x number of nodes) additively updated by result
    //! \param[in]  fact   Factor to scale with
    //! \param[in]  op     Node functor for access of nodal quantity
    template< typename OP >
    void updateNodalQuantity_( MatDimNF & quanU, const double fact, OP op ) const;

public:
    //! Return nodal displacements
    //!
    //! \param[out]   u   Matrix with the nodal mid-surface displacements as columns
    //!                   \f$[u^A{}_K]\f$
    void nodalDisplacements( MatDimNF & u ) const;

    //! Return nodal displacement increments
    //!
    //! \param[out]  deltaU    Matrix with nodal mid-surface displacement increments
    //!                        \f$[\delta{u}^A{}_K]\f$
    void nodalIncrements( MatDimNF & deltaU ) const;

    //! Interpolate current displacement at local point
    //!
    //! \param[in]      xi     Local point
    //! \return                Interpolated displacement
    eigenX::VectorSd< BELEMENT::dim >
    giveDisplacement( const VecLDim & xi ) const;

    //! Return scaled normal to current configuration whose length carries an
    //! area (and thickness) scale to the reference configuration
    eigenX::VectorSd< BELEMENT::dim >
    giveScaledNormal( const VecLDim & xi ) const;

    //! Interpolate current displacement at element vertex
    //!
    //! \param[in]      index  Local node ID of vertex
    //! \return                Interpolated displacement
    eigenX::VectorSd< BELEMENT::dim >
    giveDisplacementAtNode( const unsigned index ) const;

    //! Compute tangent to limit surface at node
    //!
    //! Note: These tangents are element-specific, and may differ
    //!       from tangents of neighbouring element.
    //!
    //! \param[in]   conf          Configuration : reference or current 
    //! \param[in]   ind           Local node index 0<=ind<numVertices
    //! \param[out]  tang0         1st tangent vector
    //! \param[out]  tang1         2nd tangent vector
    //! \param[out]  tang0GradDof  Derivative of 1st tangent w.r.t. DOFs
    //! \param[out]  tang1GradDof  Derivative of 2nd tangent w.r.t. DOFs
    void computeTangentsAtNode( const enum gshell::fem::config conf,
                                const unsigned ind,
                                VecDim & tang0,
                                VecDim & tang1,
                                Mat * tang0GradDof,
                                Mat * tang1GradDof );

    //! Compute normal to limit surface at node
    //!
    //! Note: The tangent is element-specific, and may differ
    //!       from normal of neighbouring element (non-manifold shell)
    //!
    //! \param[in]   conf           Configuration : reference or current 
    //! \param[in]   ind            Local node index
    //! \param[out]  normal         Normal
    //! \param[out]  normalGradDof  Derivative of 2nd tangent w.r.t. DOFs
    void computeNormalAtNode( const enum gshell::fem::config conf,
                              const unsigned ind,
                              VecDim & normal,
                              Mat * normalGradDof );

    //@}

protected:

    //! @name Kinematics
    //@{

    //! Resize and (rows x NF)-matrix of ENTRYs and optionally set to zero
    //!
    //! \tparam         ENTRY        sub-entry type
    //!                              (can be fixed-sized matrix, fixed-sized vector, double)
    //!
    //! \param[in,out]  mat          (rows x NF)-matrix
    //! \param[in]      rows         Number of rows
    //! \param[in]      setToZero    Initialise to zero if true
    template< typename ENTRY >
    void resizeMatrixNF_( Eigen::Matrix<ENTRY, Eigen::Dynamic, Eigen::Dynamic> & mat,
                          const unsigned rows,
                          bool setToZero ) const;

    //! Write (rows x NF)-matrix of ENTRYs to stream
    //!
    //! \tparam         ENTRY        sub-entry type
    //!                              (can be fixed-sized matrix, fixed-sized vector, double)
    //!
    //! \param[in,out]  os           Output stream
    //! \param[in]      mat          (rows x NF)-matrix
    //! \param[in]      desc         Descriptive text, will be written in front of matrix
    //! \return                      Output stream
    template< typename ENTRY >
    std::ostream & writeMatrixNF_( std::ostream & os,
                                   const Eigen::Matrix<ENTRY, Eigen::Dynamic, Eigen::Dynamic> & mat,
                                   const std::string desc = "" ) const;

    //! Compute mid-surface vectors and their derivatives w.r.t. nodal DOFs
    //!
    //! \param[in]  conf       Configuration
    //! \param[in]  xi         Local evaluation point
    //!                          \f$\xi^\alpha\f$
    //! \param[in]  phi        Shape functions
    //!                          \f$[\phi^K(\xi^\alpha)]\f$ 
    //! \param[in]  dPhiDXi    Local 1st derivative of shape functions
    //!                          \f$[\phi^K{}_{,\alpha}(\xi^\gamma)]\f$ 
    //! \param[in]  ddPhiDDXi  Local wnd derivative of shape functions
    //!                          \f$[\phi^K{}_{,\alpha\beta}(\xi^\gamma)]\f$ 
    //! \param[out] aCur       Co-variant base in mid-surface
    //!                           \f$[a^A{}_i]\f$
    //! \param[out] bCur       Co-variant local derivative of base
    //!                          \f$[b^A{}_{j\alpha}]=[a^A{}_{j,\alpha}]\f$
    //! \param[out] pCur       Vector \f$\vec{p}=\vec{a}_1\times\vec{a}_2\f$
    //!                           \f$[p^A{}]\f$
    //! \param[out] pLGrad     Local gradient of \f$\vec{p}\f$
    //!                           \f$[p^A{}_{,\beta}]\f$
    void vectors_( const enum gshell::fem::config conf,
                   const VecLDim & xi,
                   const VecNF & phi,
                   const MatLDimNF & dPhiDXi,
                   const MatSDimNF & ddPhiDDXi,
                   double aCur[dim][dim],
                   double bCur[dim][localDim][dim],
                   double pCur[dim],
                   double pLGrad[localDim][dim] ) const;

    //! Compute 1st derivatives of mid-surface vectors w.r.t. nodal DOFs
    //!
    //! \param[in]  d            Index \f$D\f$ and/or \f$\delta\f$
    //! \param[in]  l            Node index \f$L\f$
    //! \param[in]  xi           Local evaluation point
    //!                            \f$\xi^\alpha\f$
    //! \param[in]  phi          Shape functions
    //!                            \f$[\phi^K(\xi^\alpha)]\f$ 
    //! \param[in]  dPhiDXi      Local 1st derivative of shape functions
    //!                            \f$[\phi^K{}_{,\alpha}(\xi^\gamma)]\f$ 
    //! \param[in]  ddPhiDDXi    Local wnd derivative of shape functions
    //!                            \f$[\phi^K{}_{,\alpha\beta}(\xi^\gamma)]\f$ 
    //! \param[in]  aRef         Reference co-variant base in mid-surface
    //!                            \f$[\bar{a}^A{}_i]\f$
    //! \param[in]  aCur         Current co-variant base in mid-surface
    //!                            \f$[a^A{}_i]\f$
    //! \param[in]  bCur         Co-variant local derivative of base
    //!                            \f$[b^A{}_{j\alpha}]=[a^A{}_{j,\alpha}]\f$
    //! \param[in]  pCur         Vector \f$\vec{p}=\vec{a}_1\times\vec{a}_2\f$
    //!                            \f$[p^A{}]\f$
    //! \param[in]  pLGrad       Local gradient of \f$\vec{p}\f$
    //!                            \f$[p^A{}_{,\beta}]\f$
    //! \param[out] aDU_dl       1st derivative
    //!                            \f$[a^A{}_{j,D}{}^L]\f$
    //! \param[out] bDU_dl       1st derivative
    //!                            \f$[a^A{}_{j,\beta,D}{}^L]\f$
    //! \param[out] pDU_dl       1st derivative
    //!                            \f$[p^A{}_{,D}{}^L]\f$
    //! \param[out] pLGradDU_dl  1st derivative
    //!                            \f$[p^A{}_{,\beta,D}{}^L]\f$
    void vectorsGradCur_( const unsigned & d,
                          const unsigned & l,
                          const VecLDim & xi,
                          const VecNF & phi,
                          const MatLDimNF & dPhiDXi,
                          const MatSDimNF & ddPhiDDXi,
                          const double aRef[dim][dim],
                          const double aCur[dim][dim],
                          const double bCur[dim][localDim][dim],
                          const double pCur[dim],
                          const double pLGrad[localDim][dim],
                          double aDU_dl[dim][dim],
                          double bDU_dl[dim][localDim][dim],
                          double pDU_dl[dim],
                          double pLGradDU_dl[localDim][dim] ) const;

    //! Computes the Green-Lagrange strain resultants
    //! \f$E_{ij} = \alpha_{ij} + \xi^3\beta_{ij}\f$
    //! (components in reference curvi-linear system \f$\{\bar{\vec{g}}^i\}\f$)
    //!
    //! \param[in]  aRef     Co-variant base in reference configuration
    //!                        \f$[\bar{a}^A{}_{j}]\f$
    //! \param[in]  cRef     Local gradient of thickness director in reference configuration
    //!                        \f$[\bar{c}^A{}_{\alpha}]=[\bar{a}^A{}_{3,\alpha}]\f$
    //! \param[in]  aCur     Co-variant base in current configuration
    //!                        \f$[a^A{}_{j}]\f$
    //! \param[in]  cCur     Local gradient of thickness director in current configuration
    //!                        \f$[c^A{}_{\alpha}]=[a^A{}_{3,\alpha}]\f$
    //! \param[out] alpha    Membrane strain resultant (sym)
    //!                        \f$[\alpha_{ij}]\f$
    //! \param[out] beta     Bending strain resultant (sym)
    //!                        \f$[\beta_{ij}]\f$
    void strainRes_( const double aRef[dim][dim],
                     const double cRef[localDim][dim],
                     const double aCur[dim][dim],
                     const double cCur[localDim][dim],
                     double alpha[dim][dim],
                     double beta[dim][dim] ) const;

    //! Computes the Green-Lagrange strain resultants
    //! and their derivatives w.r.t. nodal DOFs
    //!
    //! The derivatives of the membrane strains w.r.t. displacements are
    //!\f[
    //!     \eta_{ij}^\mathrm{u}{}_D{}^L
    //!     = \frac{\partial \alpha_{ij}}{\partial u^D{}_L}
    //!\f]
    //!
    //! The derivatives of the bending strains w.r.t. displacements are
    //!\f[
    //!     \zeta_{ij}^\mathrm{u}{}_D{}^L
    //!     = \frac{\partial \beta_{ij}}{\partial u^D{}_L}
    //!\f]
    //! (components in reference curvi-linear system \f$\{\bar{\vec{g}}^i\}\f$)
    //!
    //! \param[in]  d        Dimension index \f$D\f$ and/or \f$\delta\f$
    //! \param[in]  l        Support node index \f$L\f$
    //! \param[in]  aCur     Co-variant base in current configuration
    //!                        \f$[a^A{}_{j}]\f$
    //! \param[in]  aDU_dl   1st derivative of co-varient base w.r.t. nodal displacements
    //!                        \f$[a^A{}_{j,D}{}^L]\f$
    //! \param[in]  cCur     Local gradient of director
    //!                        \f$[c^A{}_{\beta}]=[a^A_{2,\beta}]\f$
    //! \param[in]  cDU_dl   1st derivative of local gradient of director
    //!                      w.r.t. nodal displacements
    //!                        \f$[c^A{}_{\beta,D}{}^L]=[a^A_{2,\beta,D}{}^L]\f$
    //! \param[out] etaU_dl  Membrane strain derivative w.r.t. mid-plane displacements
    //!                        \f$[\eta_{ij}^\mathrm{u}{}_D{}^L]\f$
    //! \param[out] zetaU_dl Bending strain derivative w.r.t. nodal mid-plane displacements
    //!                        \f$[\zeta_{ij}^\mathrm{u}{}_D{}^L]\f$
    void strainResGradCur_( const unsigned & d,
                            const unsigned & l,
                            const double aCur[dim][dim],
                            const double aDU_dl[dim][dim],
                            const double cCur[localDim][dim],
                            const double cDU_dl[localDim][dim],
                            double etaU_dl[dim][dim],
                            double zetaU_dl[dim][dim] ) const;

    //@}

    //! @name Material
    //@{

    //! Computes the elasticity 4-tensor \f$H^{ijkl}\f$
    //! (components in co-variant reference base \f$\{\bar{\vec{g}}_i\}\f$)
    //!
    //! <h4>Hints:</h4>
    //! - The elasticity tensor is the constant 4-tensor of
    //!   St Venant-Kirchhoff material.
    //! - The result assumes \f$\bar{g}^{ij}=\bar{a}^{ij}\f$. Strictly, this
    //!   is only exactly satisfied, if the initial shell is flat.
    //! - The parameters are based on the plane stress \f$S^{33}=0\f$
    //!   version of St Venant-Kirchhoff material. The parameters
    //!   are computed by eliminating \f$S^{33}\f$. This affects the factor
    //!   <i>fact</i>.
    //!
    //! \param[in]  aRef      Co-variant base in reference configuration
    //!                         \f$[\bar{a}^A{}_{j}]\f$
    //! \param[out] hRef      Elasticity tensor
    //!                         \f$H^{ijkl}\f$
    void elasticityTensorVK_( const double aRef[dim][dim],
                              double hRef[dim][dim][dim][dim] ) const;

    //! Computes the 2nd Piola-Kirchhoff stress resultants
    //! \f$n^{ij}\f$ and \f$m^{ij}\f$
    //! (components in co-variant reference base \f$\{\bar{\vec{g}}_i\}\f$)
    //!
    //! \param[in]  aRef     Co-variant base in reference configuration
    //!                        \f$[\bar{a}^A{}_{j}]\f$
    //! \param[in]  alphaCur Membrane strain resultant (sym)
    //!                        \f$[\alpha_{ij}]\f$
    //! \param[in]  betaCur  Bending strain resultant (sym)
    //!                        \f$[\beta_{ij}]\f$
    //! \param[out] nCur     Membrane stress resultants
    //!                        \f$[n^{ij}]\f$
    //! \param[out] mCur     Bending stress resultants
    //!                        \f$[m^{ij}]\f$
    void stressResVK_( const double aRef[dim][dim],
                       const double alphaCur[dim][dim],
                       const double betaCur[dim][dim],
                       double nCur[dim][dim],
                       double mCur[dim][dim] ) const;

    //@}

public:
    //! @name Mechanical forces and stiffness
    //@{

    //! Compute nodal forces due to an applied body force.
    //!
    //! This function is given
    //! the coordinates and weight of a quadrature point. Moreover, a function object
    //! is passed with the operator of type VecDim = func(VecDim), where the 
    //! argument refers to the reference coordinate of the integration point. 
    //!
    //! \tparam     FUNC     Type of function object
    //!
    //! \param[in]  xi       Local coordinate of quadrature point
    //! \param[in]  weight   Weight of quadrature point
    //! \param[in]  func     Body load function (physical coord => force density)
    //! \param[in]  factor   Scalar multiplier
    //! \param[in]  conf     Define to which configuration the load density is expressed
    //! \param[in]  isThicknessIntegrated  If true, then the load function is an area density
    //!                                    otherwise it is a volumetric density.
    template< typename FUNC >
    void bodyForceDetailed( const VecLDim & xi, const double & weight,
                            FUNC func, const double & factor,
                            const enum gshell::fem::config conf,
                            const bool isThicknessIntegrated );

    //! Compute nodal forces due to an applied body load density.
    //!
    //! Note: - The body load density is with respect to the reference configuration.
    //!       - It's a volumetric density, so it will be integrated in thickness too.
    //!
    //! \tparam     FUNC     Type of function object
    //!
    //! \param[in]  xi       Local coordinate of quadrature point
    //! \param[in]  weight   Weight of quadrature point
    //! \param[in]  func     Material density load function (physical coord => force density N/m^3)
    //! \param[in]  factor   Scalar multiplier
    template< typename FUNC >
    void bodyForce( const VecLDim & xi, const double & weight,
                    FUNC func, const double & factor )
    {
        this -> bodyForceDetailed( xi, weight, func, factor,
                                   gshell::fem::REFERENCE, false );
    }

    template< typename FUNC >
    void bodyForceLin( const VecLDim & xi, const double & weight,
                       FUNC func, const double & factor );

    template< typename FUNC >
    void bodyForceLinStiffness( const VecLDim & xi, const double & weight,
                                FUNC func, const double & factor,
                                Eigen::MatrixXd & elemStiff );

    //! Compute nodal forces due to an applied (normal) pressure load.
    //!
    //! This function is given
    //! the coordinates and weight of a quadrature point. Moreover, a function object
    //! is passed with the operator of type double = func(VecDim), where the 
    //! argument refers to the reference coordinate of the integration point. 
    //!
    //! \tparam     FUNC     Type of function object
    //!
    //! \param[in]  xi       Local coordinate of quadrature point
    //! \param[in]  weight   Weight of quadrature point
    //! \param[in]  func     The pressure function object
    //! \param[in]  factor   Scalar multiplier
    template< typename FUNC >
    void externalPressure( const VecLDim & xi, const double & weight,
                           FUNC func, const double & factor );

    template< typename FUNC >
    void externalPressureStiffness( const VecLDim & xi, const double & weight,
                                    FUNC func, const double & factor,
                                    Eigen::MatrixXd & elemStiff );

protected:
    //! Computation of internal force due to internal stress. Evaluates the kernel
    //! function at given quadrature point \f$\xi^\alpha\f$, 
    //! multiplies the result with the weight and passes it to the nodes.
    //!
    //! Nodal internal force vector due to discretised test functions 
    //! in terms of nodal mid-plane displacements \f$(\delta{u}^1,\delta{u}^2,\delta{u}^3)\f$. 
    //! These are stored on the <i>first \f$\mathit{dim}=3\f$ entries</i> of the nodal
    //! force vector, i.e. \f$[f^\mathrm{int}_1,f^\mathrm{int}_2,f^\mathrm{int}_3]^L\f$
    //!\f[
    //!      [f^\mathrm{int}_D]^L
    //!      = \int_{(\xi^\alpha)}\Big(
    //!          n^{ij} \frac{ \partial \alpha_{ij} }{ \partial u^D{}_L } +
    //!          m^{ij} \frac{ \partial \beta_{ij} }{ \partial u^D{}_L }
    //!      \Big) \sqrt{ \det[\bar{a}_{\alpha\beta}] } d\xi^\alpha
    //!\f]
    //! where \f$\mathrm{dof}\f$ refers to the total number of DOFs per node.
    //!
    //! \param[in] xi      Local coordinate at which this function is evaluated
    //! \param[in] weight  Corresponding quadrature weight
    //! \param[in] conf    Configuration in which to integrate
    //! \param[in] factor  Additional factor to scale resulting force vectors
    void internalForceIntegrand_( const VecLDim & xi, const double & weight,
                                  const enum gshell::fem::config & conf,
                                  const double & factor );

public:
    //! Computation of internal force. Result is stored on nodes.
    void internalForceIntegrand( const VecLDim & xi, const double & weight )
    {
        this->internalForceIntegrand_( xi, weight, CURRENT, 1. );
        return;
    }

    //! The integral kernel for the computation of the element stiffness matrix.
    //!
    //! Evaluates the kernel function at the local coordinate xi, multiplies the 
    //! result with the weight and stores the result in a matrix.
    //!\f[
    //!      [k^\mathrm{int}_{(L*\mathit{dof}+D)(K*\mathit{dof}+C)}]
    //!      = \int_{(\xi^\alpha)}\Big(
    //!          \frac{ \partial \alpha_{ij} }{ \partial u^D{}_L }
    //!          H^{ijkl}
    //!          \frac{ \partial \alpha_{kl} }{ \partial u^C{}_K }
    //!          +
    //!          n^{ij} \frac{ \partial^2 \alpha_{ij} }{ \partial u^D{}_L \partial u^C{}_K }
    //!          +
    //!          \frac{ \partial \beta_{ij} }{ \partial u^D{}_L }
    //!          H^{ijkl} \frac{ \partial \beta_{kl} }{ \partial u^C{}_K }
    //!          +
    //!          m^{ij} \frac{ \partial^2 \beta_{ij} }{ \partial u^D{}_L \partial u^C{}_K }
    //!      \Big) \sqrt{ \det[\bar{a}_{\alpha\beta}] } d\xi^\alpha
    //!\f]
    //! where \f$\mathit{dof}\f$ refers to the total number of DOFs per node.
    //!
    //! \param[in]      xi         Local coordinate at which to evaluate
    //!                              \f$\xi^\alpha_\mathrm{GP}\f$
    //! \param[in]      weight     Corresponding quadrature weight
    //!                              \f$w_\mathrm{GP}\f$
    //! \param[in,out]  elemStiff  Element stiffness updated with contribution
    //!                              of quadrature point
    void stiffnessIntegrand( const VecLDim & xi, const double & weight,
                             Eigen::MatrixXd & elemStiff ) const;
    //@}

    //! @name Direct material operations
    //@{
    //! Cache material
    void cacheMaterial( Material * m ) { material_ = m; return; }
    //! Give constant of Elasticity
    double giveYoungsModulus( ) const { return material_->estimatedElasticityConst( ); }
    //! Give Poisson's ratio
    double givePoissonRatio( ) const { return material_->estimatedPoissonRatio( ); }
    //@}

public:
    //! @name Output and information
    //@{
    
    //! Compute strain resultants -- for output
    //!
    //! \param[in]   planeMode    In-plane (membrane), or out-of-plane (bending)
    //! \param[in]   strainMode   Type of strain measure
    //! \param[in]   xi           Local evaluation point
    //! \return                   Strain tensor in global Cartesian components
    Mat3x3 giveStrainResultant( const enum plane planeMode, 
                                const enum strains strainMode,
                                const VecLDim & xi ) const;
    
    //! Compute stress resultants -- for output
    //!
    //! \param[in]   planeMode    In-plane (membrane), or out-of-plane (bending)
    //! \param[in]   stressMode   Type of stress measure
    //! \param[in]   xi           Local evaluation point
    //! \return                   Stress tensor in global Cartesian components
    Mat3x3 giveStressResultant( const enum plane planeMode, 
                                const enum stresses stressMode,
                                const VecLDim & xi ) const;

    //@}

protected:
    //! Pointer to material object
    Material *    material_;
    //! Shell thickness
    double        thickness_;

};


//------------------------------------------------------------------------------
#include "ElementKLStatic.ipp"
//------------------------------------------------------------------------------
#endif
