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

#ifndef gshell_fem_elementstatic_h
#define gshell_fem_elementstatic_h

//------------------------------------------------------------------------------
#include <functional>
#include <array>

#include <Eigen/Core>

#include <corlib/eigenX.hpp>
#include <corlib/linalg.hpp>
#include <corlib/Accessor.hpp>
#include <corlib/ElementBasic.hpp>

#include <gshell/fem/linalg.hpp>
#include <gshell/fem/Config.hpp>
#include <gshell/fem/NodeStatic.hpp>
#include <gshell/fem/SurfaceGeometer.hpp>

//------------------------------------------------------------------------------
namespace gshell {
    namespace fem {
        
        template< typename BELEMENT, typename MAT >
        class ElementStatic;

        namespace eigenX = corlib::eigenX;
    }
}

//------------------------------------------------------------------------------
/** \defgroup elements Group of Element classes for g-shells                   */
//------------------------------------------------------------------------------


//------------------------------------------------------------------------------
/** \ingroup elements
 * \brief Element for generalised-shells (i.e. shear-flexible),
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
class gshell::fem::ElementStatic : public BELEMENT
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
        INPLANE,        //!< in-plane, i.e. membrane
        OUTOFPLANE      //!< out-of-plane, i.e. bending
    };

    //! Strain measures
    enum strains {
        LOCALGL,        //!< Green-Lagrange in local co-ordinates
        GREENLAGRANGE,  //!< Green-Lagrange in global co-ordinates
        EULERALMANSI    //!< Euler-Almansi in global co-ordinates
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
    //! Total number of DOFs at nodes -- always 5
    static const unsigned      dof       = Node::dof;  // always 5
    //! Dimension of shell manifold -- always 2
    static const unsigned      localDim  = BasisElement::localDim;
    //! Number of effective second derivatives subtracting duplicities due to symmetry -- always 3
    static const unsigned      sDim      = BasisElement::sDim;  // always 3
    //! Map to get Voigt indices from pair of local indices
    //! \f$(\alpha,\beta) \mapsto \mathcal{A}\f$
    //by \f$[localDim*\alpha+\beta]\f$ for any \f$\alpha,\beta=1,2\f$
    static const unsigned      voigtForward[localDim][localDim];// = { { 0, 1 }, { 1, 2 } };

    typedef typename BasisElement::VecDim           VecDim;
    typedef eigenX::VectorSd<dof>                   VecDof;
    typedef typename BasisElement::VecLDim          VecLDim;
    typedef eigenX::MatrixSd<3, 3>                  Mat3x3;

    // Note: The following type definitions are all for the same dynamic matrix
    //       type. The idea is to easier understand the dimensions (and thus
    //       content) of the matrices. However, the actual size has to be
    //       defined explicitly.
    typedef typename BasisElement::VecNF                   VecNF;
    typedef typename BasisElement::MatDimNF                MatDimNF;
    typedef Eigen::MatrixXd                                MatDofNF;
    typedef typename BasisElement::MatLDimNF               MatLDimNF;
    typedef typename BasisElement::MatSDimNF               MatSDimNF;
    typedef Eigen::MatrixXd                                Mat;

private:
    typedef ElementStatic< BELEMENT, MAT >                 ElementStatic_;

public:
    //! Constructor
    //!
    //! Note:
    //!    - Pointer to element's shape function has to be set
    //!    - The BasisElement::suppportNodes (nodes carrying DOFs) have to assigned
    ElementStatic() :
        BasisElement(),
        material_( NULL ),
        thickness_( 0.0 )
    {
        static_assert( dim == 3 );
        static_assert( dof == 5 );
        static_assert( localDim == 2 );
        return;
    }

    //! Destructor
    ~ElementStatic( ) { return; }

    //! Set thickness
    //! \param[in]  thickness   Provided thickness value, constant on element
    void setThickness( const double thickness ) { thickness_ = thickness; return; }

protected:
    //! Accessor to a general nodal quantitiy accessed by the pass operator
    //!
    //! \tparam     OP     Node functor which access the desired nodal quantity
    //!
    //! \param[out] quan0  Matrix (dim x number of nodes) containing the result
    //! \param[out] quan1  Matrix (dim x number of nodes) containing the result
    //! \param[in]  op     Node functor for access of nodal quantity
    template< typename OP >
    void nodalQuantity_( MatDimNF & quanU, MatDimNF & quanW,
                         OP op ) const;

public:
    //! @name Accessors of nodal quantities
    //@{

    //! Accessor to the indices of the nodal degrees of freedom
    //!
    //! \param[in]  dofIndices Storage vector for indices
    void getDofIndices( std::vector<unsigned> & dofIndices ) const;

    //! Return nodal displacements 
    //!
    //! \param[out]   u   Matrix with the nodal mid-surface displacements as columns
    //!                   \f$[u^A{}_K]\f$
    void nodalDisplacements( MatDimNF & u ) const;

    //! Return nodal displacements and shears
    //!
    //! \param[out]   u   Matrix with the nodal mid-surface displacements as columns
    //!                   \f$[u^A{}_K]\f$
    //! \param[out]   w   Matrix with the nodal shear vectors as columns
    //!                   \f$[w^A{}_K]\f$
    void nodalDisplacementsAndShears( MatDimNF & u, MatDimNF & w ) const;

    //! Return nodal displacement increments
    //!
    //! \param[out]  deltaU    Matrix with nodal mid-surface displacement increments
    //!                        \f$[\delta{u}^A{}_K]\f$
    //! \param[out]  deltaW    Matrix with nodal shear vector increments
    //!                        \f$[\delta{w}^A{}_K]\f$
    void nodalIncrements( MatDimNF & deltaU, MatDimNF & deltaW ) const;

    //! Interpolate current displacement at local point
    //!
    //! \param[in]      xi     Local point
    //! \return                Interpolated displacement
    eigenX::VectorSd<BELEMENT::dim>
    giveDisplacement( const VecLDim & xi ) const;

    //! Interpolate current displacement at element vertex
    //!
    //! \param[in]      index  Local node ID of vertex
    //! \return                Interpolated displacement
    eigenX::VectorSd<BELEMENT::dim>
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

    //! Transform convected nodal shear vector to nodal Cartesian shear vector
    //!
    //! \param[in,out]  wNF          Nodal shear vector
    void nodalCartesianShearsCur_( MatDimNF & wNF ) const;

    //! Determine coefficients (derivatives) of tangents at element nodes w.r.t. node #l
    //!
    //! \param[in]      l            Element node index
    //! \param[out]     tangCoeff_l  Derivative of tangents at element nodes w.r.t. node #l
    void tangentsGradCur_( const unsigned & l, MatLDimNF & tangCoeff_l ) const;

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

    //! Compute vectors and their derivatives w.r.t. nodal DOFs
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
    //! \param[out] wCur       Shear vector
    //!                          \f$[w^A]\f$
    //! \param[out] zCur       Local gradient of shear vector
    //!                          \f$[z^A{}_{\alpha}]=[w^A{}_{,\alpha}]\f$
    //! \param[out] pCur       Vector \f$\vec{p}=\vec{a}_1\times\vec{a}_2\f$
    //!                          \f$[p^A{}]\f$
    //! \param[out] pLGrad     Local gradient of \f$\vec{p}\f$
    //!                          \f$[p^A{}_{,\beta}]\f$
    //! \param[out] qCur       Vector \f$\vec{q}=\vec{p}+\vec{w}\f$
    //!                          \f$[q^A]\f$
    //! \param[out] qLGrad     Local gradient of \f$\vec{q}\f$
    //!                          \f$[q^A{}_{,\beta}]\f$
    //! \param[out] dCur       Thickness director
    //!                          \f$[d^A]\f$
    //! \param[out] cCur       Local gradient of thickness director
    //!                          \f$[c^A{}_{\alpha}]=[d^A{}_{,\alpha}]\f$
    void vectors_( const enum gshell::fem::config conf,
                   const VecLDim & xi,
                   const VecNF & phi,
                   const MatLDimNF & dPhiDXi,
                   const MatSDimNF & ddPhiDDXi,
                   double aCur[dim][dim],
                   double bCur[dim][localDim][dim],
                   double wCur[dim],
                   double zCur[localDim][dim],
                   double pCur[dim],
                   double pLGrad[localDim][dim],
                   double qCur[dim],
                   double qLGrad[localDim][dim],
                   double dCur[dim],
                   double cCur[localDim][dim] ) const;

    //! Compute vectors and their derivatives w.r.t. nodal DOFs
    //!
    //! Nodal DOFs are \f$x^D{}_L\f$ and \f$w^\delta{}_L\f$
    //!
    //! \param[in]  d            Index \f$D\f$ and/or \f$\delta\f$
    //! \param[in]  l            Node index \f$L\f$
    //! \param[in]  xi           Local evaluation point
    //!                            \f$\xi^\alpha\f$
    //! \param[in]  phi          Shape functions
    //!                            \f$[\phi^K(\xi^\alpha)]\f$ 
    //! \param[in]  dPhiDXi      Local 1st derivative of shape functions
    //!                            \f$[\phi^K{}_{,\alpha}(\xi^\gamma)]\f$ 
    //! \param[in]  ddPhiDDXi    Local 2nd derivative of shape functions
    //!                            \f$[\phi^K{}_{,\alpha\beta}(\xi^\gamma)]\f$ 
    //! \param[in]  tangCoeff_l  Derivative of nodal tangent (effective components)
    //!                            \f$[\ell_{\alpha K}{}^L\f$
    //! \param[in]  aRef         Reference co-variant base in mid-surface
    //!                            \f$[\bar{a}^A{}_i]\f$
    //! \param[in]  aCur         Current co-variant base in mid-surface
    //!                            \f$[a^A{}_i]\f$
    //! \param[in]  bCur         Co-variant local derivative of base
    //!                            \f$[b^A{}_{j\alpha}]=[a^A{}_{j,\alpha}]\f$
    //! \param[in]  wCur         Shear vector
    //!                            \f$[w^A]\f$
    //! \param[in]  zCur         Local gradient of shear vector
    //!                            \f$[z^A{}_{\alpha}]=[w^A{}_{,\alpha}]\f$
    //! \param[in]  pCur         Vector \f$\vec{p}=\vec{a}_1\times\vec{a}_2\f$
    //!                            \f$[p^A{}]\f$
    //! \param[in]  pLGrad       Local gradient of \f$\vec{p}\f$
    //!                            \f$[p^A{}_{,\beta}]\f$
    //! \param[in]  qCur         Vector \f$\vec{q}=\vec{p}+\vec{w}\f$
    //!                            \f$[q^A]\f$
    //! \param[in]  qLGrad       Local gradient of \f$\vec{q}\f$
    //!                            \f$[q^A{}_{,\beta}]\f$
    //! \param[in]  dCur         Thickness director
    //!                            \f$[d^A]\f$
    //! \param[in]  cCur         Local gradient of thickness director
    //!                            \f$[c^A{}_{\alpha}]=[d^A{}_{,\alpha}]\f$
    //! \param[out] pDU_dl       1st derivative
    //!                            \f$[p^A{}_{,D}{}^L]\f$
    //! \param[out] pLGradDU_dl  1st derivative
    //!                            \f$[p^A{}_{,\beta,D}{}^L]\f$
    //! \param[out] qDU_dl       1st derivative
    //!                            \f$[q^A{}_{,D}{}^L]\f$
    //! \param[out] qDW_dl       1st derivative 
    //!                            \f$[q^A{}_{,\delta}{}^L]\f$
    //! \param[out] qLGradDU_dl  1st derivative
    //!                            \f$[q^A{}_{,\beta,D}{}^L]\f$
    //! \param[out] qLGradDW_dl  1st derivative
    //!                            \f$[q^A{}_{,\beta,\delta}{}^L]\f$
    //! \param[out] dDU_dl       1st derivative
    //!                            \f$[d^A{}_{,D}{}^L]\f$
    //! \param[out] dDW_dl       1st derivative
    //!                            \f$[d^A{}_{,\delta}{}^L]\f$
    //! \param[out] cDU_dl       1st derivative
    //!                            \f$[d^A{}_{,\beta,D}{}^L]\f$
    //! \param[out] cDW_dl       1st derivative
    //!                           \f$[d^A{}_{,\beta,\delta}{}^L]\f$
    void vectorsGradCur_( const unsigned & d,
                          const unsigned & l,
                          const VecLDim & xi,
                          const VecNF & phi,
                          const MatLDimNF & dPhiDXi,
                          const MatSDimNF & ddPhiDDXi,
                          const MatLDimNF & tangCoeff_l,
                          const double aRef[dim][dim],
                          const double aCur[dim][dim],
                          const double bCur[dim][localDim][dim],
                          const double wCur[dim],
                          const double zCur[localDim][dim],
                          const double pCur[dim],
                          const double pLGrad[localDim][dim],
                          const double qCur[dim],
                          const double qLGrad[localDim][dim],
                          const double dCur[dim],
                          const double cCur[localDim][dim],
                          double pDU_dl[dim],
                          double pLGradDU_dl[localDim][dim],
                          double qDU_dl[dim],
                          double qDW_dl[dim],
                          double qLGradDU_dl[localDim][dim],
                          double qLGradDW_dl[localDim][dim],
                          double dDU_dl[dim],
                          double dDW_dl[dim],
                          double cDU_dl[localDim][dim],
                          double cDW_dl[localDim][dim] ) const;

    //! Computes the Green-Lagrange strain resultants
    //! \f$E_{ij} = \alpha_{ij} + \xi^3\beta_{ij}\f$
    //! (components in reference curvi-linear system \f$\{\bar{\vec{g}}^i\}\f$)
    //!
    //! \param[in]  aRef     Co-variant base in reference configuration
    //!                        \f$[\bar{a}^A{}_{j}]\f$
    //! \param[in]  dRef     Thickness director in reference configuration
    //!                        \f$[\bar{d}^A]\f$
    //! \param[in]  cRef     Local gradient of thickness director in reference configuration
    //!                        \f$[\bar{c}^A{}_{\alpha}]=[d^A{}_{,\alpha}]\f$
    //! \param[in]  aCur     Co-variant base in current configuration
    //!                        \f$[a^A{}_{j}]\f$
    //! \param[in]  dCur     Thickness director
    //!                        \f$[d^A]\f$
    //! \param[in]  cCur     Local gradient of thickness director in current configuration
    //!                        \f$[c^A{}_{\alpha}]=[d^A{}_{,\alpha}]\f$
    //! \param[out] alpha    Membrane strain resultant (sym)
    //!                        \f$[\alpha_{ij}]\f$
    //! \param[out] beta     Bending strain resultant (sym)
    //!                        \f$[\beta_{ij}]\f$
    void strainRes_( const double aRef[dim][dim],
                     const double dRef[dim],
                     const double cRef[localDim][dim],
                     const double aCur[dim][dim],
                     const double dCur[dim],
                     const double cCur[localDim][dim],
                     double alpha[dim][dim],
                     double beta[dim][dim] ) const;


    //! Computes Green-Lagrange strain resultants and their derivatives w.r.t. nodal DOFs
    //!
    //! The derivatives of the membrane strains w.r.t. displacements
    //! and shear vectors are
    //!\f[
    //!     \eta_{ij}^\mathrm{u}{}_D{}^L
    //!     = \frac{\partial \alpha_{ij}}{\partial u^D{}_L}
    //!     \qquad
    //!     \eta_{ij}^\mathrm{w}{}_D{}^L
    //!     = \frac{\partial \alpha_{ij}}{\partial w^D{}_L}
    //!\f]
    //!
    //! The derivatives of the bending strains w.r.t. displacements
    //! and shear vectors are
    //!\f[
    //!     \zeta_{ij}^\mathrm{u}{}_D{}^L
    //!     = \frac{\partial \beta_{ij}}{\partial u^D{}_L}
    //!     \qquad
    //!     \zeta_{ij}^\mathrm{w}{}_D{}^L
    //!     = \frac{\partial \beta_{ij}}{\partial w^D{}_L}
    //!\f]
    //! (components in reference curvi-linear system \f$\{\bar{\vec{g}}^i\}\f$)
    //!
    //! \param[in]  d        Dimension index \f$D\f$ and/or \f$\delta\f$
    //! \param[in]  l        Support node index \f$L\f$
    //! \param[in]  aCur     Co-variant base in current configuration
    //!                        \f$[a^A{}_{j}]\f$
    //! \param[in]  dCur     Thickness director
    //!                        \f$[d^A]\f$
    //! \param[in]  dDU_dl   1st derivative of director
    //!                      w.r.t. nodal displacements
    //!                        \f$[d^A_{,D}{}^L]\f$
    //! \param[in]  dDW_dl   1st derivative of director
    //!                      w.r.t. nodal shears
    //!                        \f$[d^A_{,\delta}{}^L]\f$
    //! \param[in]  cCur     Local gradient of director
    //!                        \f$[c^A_{\beta}]=[d^A{}_{,\beta}]\f$
    //! \param[in]  cDU_dl   1st derivative of local gradient of director
    //!                      w.r.t. nodal displacements
    //!                        \f$[c^A_{\beta,D}{}^L]\f$
    //! \param[in]  cDW_dl   1st derivative of local gradient of director
    //!                      w.r.t. nodal shears
    //!                        \f$[c^A_{\beta,\delta}{}^L]\f$
    //! \param[in]  dPhiDXi  \f$[\phi_{,\beta}{}^K]\f$
    //! \param[out] alphaDU_dl  Membrane strain derivative w.r.t. mid-plane displacements
    //!                         \f$[\alpha_{ij}^\mathrm{u}{}_{,D}{}^L]\f$
    //! \param[out] alphaDW_dl  Membrane strain derivative w.r.t. nodal shear vectors
    //!                         \f$[\alpha_{ij}^\mathrm{w}{}_{,\delta}{}^L]\f$
    //! \param[out] betaDU_dl   Bending strain derivative w.r.t. nodal mid-plane displacements
    //!                         \f$[\beta_{ij}^\mathrm{u}{}_{,D}{}^L]\f$
    //! \param[out] betaDW_dl   Bending strain derivative w.r.t. nodal shear vectors
    //!                         \f$[\beta_{ij}^\mathrm{w}{}_{,\delta}{}^L]\f$
    void strainResGradCur_( const unsigned & d,
                            const unsigned & l,
                            const double aCur[dim][dim],
                            const double dCur[dim],
                            const double dDU_dl[dim],
                            const double dDW_dl[dim],
                            const double cCur[localDim][dim],
                            const double cDU_dl[localDim][dim],
                            const double cDW_dl[localDim][dim],
                            const MatLDimNF & dPhiDXi,
                            double alphaDU_dl[dim][dim],
                            double alphaDW_dl[dim][dim],
                            double betaDU_dl[dim][dim],
                            double betaDW_dl[dim][dim] ) const;

    //@}

    //! @name Material
    //@{

    //! Computes the elasticity 4-tensor \f$H^{ijkl}\f$
    //!
    //! <h4>Hints:</h4>
    //! - components in co-variant reference base \f$\{\bar{\vec{g}}_i\}\f$
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

    //! Computes 2nd Piola-Kirchhoff stress resultants \f$n^{ij}\f$ and \f$m^{ij}\f$
    //!
    //! The components are in co-variant reference base \f$\{\bar{\vec{g}}_i\}\f$
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

    //! Compute nodal forces due to body force
    //!
    //! Compute nodal forces due to an applied body force. This function is given
    //! the coordinates and weight of a quadrature point. Moreover, a function object
    //! is passed with the operator of type VecDim = func(VecDim), where the 
    //! argument refers to the reference coordinate of the integration point. 
    //!
    //! \tparam     FUNC     Type of function object
    //!
    //! \param[in]  xi       Local coordinate of quadrature point
    //! \param[in]  weight   Weight of quadrature point
    //! \param[in]  func     The force function object
    //! \param[in]  factor   Scalar multiplier
    template< typename FUNC >
    void bodyForce( const VecLDim & xi, const double & weight, FUNC func,
                    const double & factor );

    template< typename FUNC >
    void bodyForceLin( const VecLDim & xi, const double & weight, FUNC func,
                       const double & factor );

    template< typename FUNC >
    void bodyForceLinStiffness( const VecLDim & xi, const double & weight, FUNC func,
                                const double & factor, Eigen::MatrixXd & elemStiff );

    //! Evaluate internal force kernel at xi (store results in nodes)
    //!
    //! Computation of internal force due to internal stress. Evaluates the kernel
    //! function at given quadrature point \f$\xi^\alpha\f$, multiplies the result with the weight
    //! and passes it to the nodes.
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
    //! The <i>second \f$\mathit{localDim}=2\f$ entries</i> of the nodal
    //! force vector, i.e. \f$[f^\mathrm{int}_4,f^\mathrm{int}_5]^L\f$,
    //! are filled by the nodal internal force vector due to discretised test
    //! functions in terms of shear vectors \f$(\delta{w}^1,\delta{w}^2)\f$.
    //!\f[
    //!      [f^\mathrm{int}_{(\mathit{dim}+\delta)}]^L
    //!      = \int_{(\xi^\alpha)}\Big(
    //!          n^{ij} \frac{ \partial \alpha_{ij} }{ \partial w^\delta{}_L } +
    //!          m^{ij} \frac{ \partial \beta_{ij} }{ \partial w^\delta{}_L }
    //!      \Big) \sqrt{ \det[\bar{a}_{\alpha\beta}] } d\xi^\alpha
    //!\f]
    //! where \f$\mathrm{dof}\f$ refers to the total number of DOFs per node.
    //!
    //! \param[in] xi      Local coordinate at which this function is evaluated
    //! \param[in] weight  Corresponding quadrature weight
    void internalForceIntegrand( const VecLDim & xi, const double & weight );

    //! Evaluate stiffness kernel at xi, weigh result and store it on #elemStiff_
    //!
    //! The integral kernel for the computation of the element stiffness matrix. 
    //! Evaluates the kernel function at the local coordinate xi, multiplies the 
    //! result with the weight and stores the result in a matrix.
    //!
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
    //! and
    //!\f[
    //!      [k^\mathrm{int}_{(L*\mathit{dof}+\mathit{dim}+\delta)(K*\mathit{dof}+\mathit{dim}+\gamma)}]
    //!      = \int_{(\xi^\alpha)}\Big(
    //!          \frac{ \partial \alpha_{ij} }{ \partial w^\delta{}_L }
    //!          H^{ijkl}
    //!          \frac{ \partial \alpha_{kl} }{ \partial w^\gamma{}_K }
    //!          +
    //!          n^{ij} \frac{ \partial^2 \alpha_{ij} }{ \partial w^\delta{}_L \partial w^\gamma{}_K }
    //!          +
    //!          \frac{ \partial \beta_{ij} }{ \partial w^\delta{}_L }
    //!          H^{ijkl} \frac{ \partial \beta_{kl} }{ \partial w^\gamma{}_K }
    //!          +
    //!          m^{ij} \frac{ \partial^2 \beta_{ij} }{ \partial w^\delta{}_L \partial w^\gamma{}_K }
    //!      \Big) \sqrt{ \det[\bar{a}_{\alpha\beta}] } d\xi^\alpha
    //!\f]
    //! and
    //!\f[
    //!      [k^\mathrm{int}_{(L*\mathit{dof}+D)(K*\mathit{dof}+\mathit{dim}+\gamma)}]
    //!      = \int_{(\xi^\alpha)}\Big(
    //!          \frac{ \partial \alpha_{ij} }{ \partial u^D{}_L }
    //!          H^{ijkl}
    //!          \frac{ \partial \alpha_{kl} }{ \partial w^\gamma{}_K }
    //!          +
    //!          n^{ij} \frac{ \partial^2 \alpha_{ij} }{ \partial u^D{}_L \partial w^\gamma{}_K }
    //!          +
    //!          \frac{ \partial \beta_{ij} }{ \partial u^D{}_L }
    //!          H^{ijkl} \frac{ \partial \beta_{kl} }{ \partial w^\gamma{}_K }
    //!          +
    //!          m^{ij} \frac{ \partial^2 \beta_{ij} }{ \partial u^D{}_L \partial w^\gamma{}_K }
    //!      \Big) \sqrt{ \det[\bar{a}_{\alpha\beta}] } d\xi^\alpha
    //!\f]
    //! etc
    //! where \f$\mathit{dof}\f$ refers to the total number of DOFs per node.
    //!
    //! \param[in]  xi      Local coordinate at which this function is evaluated \f$\xi^\alpha_\mathrm{GP}\f$
    //! \param[in]  weight  Corresponding quadrature weight \f$w_\mathrm{GP}\f$
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

public:
    //! @name Debug methods
    //@{
    
    // //! Prepare Gnuplot files containing shape functions and tangents
    // void gnuPlotShapeFunctions( const std::string & baseName,
    //                             const VecLDim & centre,
    //                             const double & size,
    //                             const VecLDim & direction,
    //                             const unsigned & numSteps,
    //                             const bool & scaleToOne );

    // //! Prepare Gnuplot files containing shape derivatives and tangents
    // void gnuPlotShapeDeriv1( const std::string & baseName,
    //                          const VecLDim & centre,
    //                          const double & size,
    //                          const VecLDim & direction,
    //                          const unsigned & numSteps,
    //                          const bool & scaleToOne );

    //! Compute mid-surface vectors (like bases, director, etc) and their
    //! 1st derivatives w.r.t. nodal DOFs. The vectors are computed for a range
    //! of values and a tangent is computed at the middle of the range as well.
    //! The vectors can be visualised with Gnuplot.
    void gnuPlotVector( const std::string vectorName,
                        const unsigned a,
                        const unsigned be,
                        const std::string fileBaseName,
                        const double mid,
                        const double width,
                        const unsigned steps,
                        const unsigned dofLow,
                        const unsigned dofUp );

    //! Compute mid-surface vectors (like bases, director, etc) and their
    //! 1st or 2nd derivatives w.r.t. nodal DOFs. The result of the derivatives
    //! is compared to finite difference counterpart
    void finiteDifferenceVector( std::ostream & out,
                                 const std::string vectorName,
                                 const unsigned a,
                                 const unsigned be,
                                 const double mid,
                                 const double width,
                                 const unsigned secondDeriv );

    //! Write nodal tangents
    void writeNodalTangents( std::ostream & out );

    //@}

protected:
    //! Pointer to material object
    Material *    material_;
    //! Shell thickness
    double        thickness_;

};


//------------------------------------------------------------------------------
#include "ElementStatic.ipp"
#include "ElementStaticDebug.ipp"
//------------------------------------------------------------------------------
#endif
