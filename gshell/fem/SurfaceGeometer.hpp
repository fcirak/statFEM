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

#ifndef gshell_fem_surfacegeometer_h
#define gshell_fem_surfacegeometer_h

#include <functional>

#include <Eigen/Core>
#include <corlib/eigenX.hpp>

#include <gshell/fem/linalg.hpp>

//------------------------------------------------------------------------------
namespace gshell {
    namespace fem {
        
        class SurfaceGeometer;

        namespace eigenX = corlib::eigenX;

    }
}



//==============================================================================
//! A bunch of methods to interpolate geometric quantities of a surface.
//! The object is only a container dedicated to offer the methods rather than
//! storing any data.
//!
//! About indices in formulas:
//! - \f$A,B,C,D,...=1,2,3\f$ are used for the global Cartesian basis
//!                         \f$\{\vec{e}_A\}\f$
//!                         and its co-ordinates \f$(x^A)\f$
//! - \f$i,j,k,...=1,2,3\f$ are used for the local curvi-linear bases 
//!                         \f$\{\vec{g}_i\}\f$ (co-variant current),
//!                         \f$\{\vec{g}^i\}\f$ (contra-variant current),
//!                         \f$\{\bar{\vec{g}}_i\}\f$ (co-variant reference),
//!                         \f$\{\bar{\vec{g}}^i\}\f$ (contra-variant refence)
//!                         and its co-ordinates \f$(\xi^i)\f$
//! - \f$\alpha,\beta,...=1,2\f$ are used for the local curvi-linear bases, e.g.
//!                         \f$\{\vec{a}_\alpha\}\f$ (co-variant current),
//!                         \f$\{\vec{a}^\alpha\}\f$ (contra-variant current),
//!                         \f$\{\bar{\vec{a}}_\alpha\}\f$ (co-variant reference),
//!                         \f$\{\bar{\vec{a}}^\alpha\}\f$ (contra-variant refence),
//!                         and its mid-plane co-ordinates \f$(\xi^\alpha)\f$
//! - \f$\mathcal{A},\mathcal{B},...=1,2,3\f$ are used for symmetric Voigt-indices of 
//!                         \f$(\alpha,\beta)\f$-indexed quantities, see also #voigtForward
//! - \f$K,L,...=1,...,NF\f$ are the local node indices attached to the
//!                         shape functions \f$N^K(\xi^\alpha)\f$ whose support
//!                         intersects with the element domain
class gshell::fem::SurfaceGeometer 
{
public:
    //! Dimension of embedding space
    static const unsigned      dim       = 3;
    //! Dimension of shell manifold
    static const unsigned      localDim  = 2;
    //! Number of effective second derivatives subtracting duplicities due
    //! to symmetry
    static const unsigned      sDim      = 3; //localDim*(localDim+1)/2;
    //! Map to get Voigt indices from pair of local indices
    //! \f$(\alpha,\beta) \mapsto \mathcal{A}\f$
    //by \f$[localDim*\alpha+\beta]\f$ for any \f$\alpha,\beta=1,2\f$
    static const unsigned      voigtForward[2][2]; //= { { 0, 1 }, { 1, 2 } };

public:
    //! Interpolate vector quantity
    //!
    //! \param[in]  numFunctions  Number of shape functions NF
    //! \param[in]  phi  Local shape functions (evaluated at a point)
    //!                    \f$[\phi^K(\xi^\alpha)]\f$,
    //!                    (NF)-vector
    //! \param[in]  qNF  Nodal quantities
    //!                    \f$[q^A{}_K]\f$,
    //!                    (dim)x(NF)-matrix in row-major flat storage
    //! \param[out] q    Interpolated quantity
    //!                    \f$q^A{}(\xi^\alpha)\f$
    static
    void interpolate( const unsigned numFunctions,
                      const double * phi,
                      const double * qNF,
                      double q[dim] );

    //! Interpolate first derivative of vector quantity w.r.t. local co-ordinates
    //!
    //! \param[in]  numFunctions  Number of shape functions NF
    //! \param[in]  dPhiDXi  Local derivatives of shape functions
    //!                        at \f$\xi^\alpha\f$ evaluated
    //!                        \f$[\phi_{,\alpha}^K(\xi^\alpha)]\f$,
    //!                        (localDim)x(NF)-matrix in row-major flat storage
    //! \param[in]  qNF      Nodal quantities
    //!                        \f$[q^A{}_K]\f$,
    //!                        (dim)x(NF)-matrix in row-major flat storage
    //! \param[out] qDeriv1  Interpolated 1st derivative of quantity
    //!                        \f$q^A{}_{,\alpha}(\xi^\alpha)\f$
    static
    void interpolateDeriv1( const unsigned numFunctions,
                            const double * dPhiDXi,
                            const double * qNF,
                            double qDeriv1[localDim][dim] );

    //! Interpolate second derivative of vector quantity w.r.t. local co-ordinates
    //!
    //! \param[in]  numFunctions  Number of shape functions NF
    //! \param[in]  ddPhiDDXi   Local 2nd derivatives of shape functions
    //!                           at \f$\xi^\alpha\f$ evaluated
    //!                           \f$[DD\phi_{\mathcal{A}}^K(\xi^\alpha)]=[\phi_{,\alpha\beta}^K(\xi^\alpha)]\f$,
    //!                           (sDim)x(NF)-matrix in row-major flat storage
    //! \param[in]  qNF         Nodal quantities
    //!                           \f$[q^A{}_K]\f$,
    //!                           (dim)x(NF)-matrix in row-major flat storage
    //! \param[out] qDeriv2     Interpolated 2nd derivative of quantity
    //!                           \f$q^A{}_{,\alpha\beta}(\xi^\alpha)\f$
    static
    void interpolateDeriv2( const unsigned numFunctions,
                            const double * ddPhiDDXi,
                            const double * qNF,
                            double qDeriv2[localDim][localDim][dim] );

    //! Computes the co-variant base vectors in mid-surface
    //! (components in global Cartesian system)
    //!
    //! \param[in]  numFunctions  Number of shape functions NF
    //! \param[in]  dPhiDXi  Local derivatives of shape functions
    //!                        \f$[\phi_{,\alpha}^K]\f$,
    //!                        (localDim)x(NF)-matrix in row-major flat storage
    //! \param[in]  xNF      Nodal positions
    //!                        \f$[x^A{}_K]\f$,
    //!                        (dim)x(NF)-matrix in row-major flat storage
    //! \param[out] a        Co-variant base
    //!                        \f$\{a^A{}_{j}\vec{e}_A\}\f$
    //! \retval              The determinant
    //!                        \f$\det[a^A{}_j]\f$
    static
    double covariantBase( const unsigned numFunctions,
                          const double * dPhiDXi,
                          const double * xNF,
                          double a[dim][dim] );

    //! Computes the derivatives of co-variant base vectors in mid-surface
    //!
    //! (components in global Cartesian system)
    //!
    //! \param[in]  numFunctions  Number of shape functions NF
    //! \param[in]  ddPhiDDXi  Local derivatives of shape functions,
    //!                          \f$[DD\phi_{\mathcal{A}}^K]=[\phi_{,\alpha\beta}^K]\f$,
    //!                          (sDim)x(NF)-matrix in row-major flat storage
    //! \param[in]  xNF        Nodal positions
    //!                          \f$[x^A{}_K]\f$,
    //!                          (dim)x(NF)-matrix in row-major flat storage
    //! \param[in]  a          Co-variant base
    //!                          \f$\{a^A{}_{j}\vec{e}_A\}\f$
    //! \param[out] b          Co-variant local derivative of base
    //!                          \f$[b^A{}_{j\alpha}]=[a^A{}_{j,\alpha}]\f$
    static 
    void covariantBaseDeriv1( const unsigned numFunctions,
                              const double * ddPhiDDXi,
                              const double * xNF,
                              const double a[dim][dim],
                              double b[dim][localDim][dim] );

    //! Computes the thickness director
    //!
    //! (components in global Cartesian system)
    //!
    //! \param[in]  p        Co-variant vector (e.g. \f$\vec{a}_1\times\vec{a}_2\f$)
    //!                        \f$[p^A]\f$
    //! \param[in]  w        Shear vector
    //!                        \f$[w^A]\f$
    //! \param[out] d        Thickness director
    //!                        \f$[d^A]\f$
    static
    void director( const double p[dim],
                   const double w[dim],
                   double d[dim] );

    //! Computes the local gradient of the thickness director
    //!
    //! (components in global Cartesian system)
    //!
    //! \param[in]  a        Co-variant base
    //!                        \f$[a^A{}_{j}]\f$
    //! \param[in]  w        Shear vector
    //!                        \f$[w^A]\f$
    //! \param[in]  d        Thickness director
    //!                        \f$[d^A]\f$
    //! \param[in]  b        Co-variant local derivative of base
    //!                        \f$[b^A{}_{j\alpha}] = [a^A{}_{j,\alpha}]\f$
    //! \param[in]  z        Local gradient of shear vector
    //!                        \f$[z^A{}_{\alpha}]=[w^A{}_{,\alpha}]\f$
    //! \param[out] c        Local gradient of thickness director
    //!                        \f$[c^A{}_{\alpha}]=[d^A{}_{,\alpha}]\f$
    static
    void directorDeriv1( const double a[dim][dim],
                         const double w[dim],
                         const double d[dim],
                         const double b[dim][localDim][dim],
                         const double z[localDim][dim],
                         double c[localDim][dim] );

    //! Compute (local)  derivative of cross-produced vector
    //! \f$\vec{p}=\vec{a}_1\times\vec{a}_2\f$
    //!
    //! \param[in]  a0         Vector \f$\vec{a}_1\f$
    //! \param[in]  a0Deriv1   Derivative
    //!                          \f$[a^A{}_{1,\beta}]\f$
    //! \param[in]  a1         Vector \f$\vec{a}_2\f$
    //! \param[in]  a1Deriv1   Derivative
    //!                          \f$[a^A{}_{1,\beta}]\f$
    //! \param[out] pDeriv1    Derivative
    //!                          \f$[p^A{}_{,\beta}]\f$
    static
    void crossprodVectorDeriv1( const double a0[dim],
                                const double a0Deriv1[localDim][dim],
                                const double a1[dim],
                                const double a1Deriv1[localDim][dim],
                                double pDeriv1[localDim][dim] );

    //! Compute 1st derivative of cross-produced vector
    //! \f$\vec{p}=\vec{a}_1\times\vec{a}_2\f$
    //! with respect to nodal DOFs \f$u^D{}_L\f$
    //!
    //! HINT: The nodal DOFs can be displacements or shears.
    //!
    //! \param[in]  a0        Vector \f$\vec{a}_1\f$
    //! \param[in]  a0DU_dl   1st derivative of \f$\vec{a}_1\f$ w.r.t.
    //!                         nodal DOFs \f$[a^A{}_{1,D}{}^L]\f$
    //! \param[in]  a1        Vector \f$\vec{a}_2\f$
    //! \param[in]  a1DU_dl   1st derivative of \f$\vec{a}_2\f$ w.r.t.
    //!                         nodal DOFs \f$[a^A{}_{2,D}{}^L]\f$
    //! \param[out] pDU_dl    1st derivative of \f$\vec{d}\f$ w.r.t.
    //!                         nodal DOFs \f$[p^A{}_{,D}{}^L]\f$
    static
    void crossprodVectorGradDof( const double a0[dim],
                                 const double a0DU_dl[dim],
                                 const double a1[dim],
                                 const double a1DU_dl[dim],
                                 double pDU_dl[dim] );

    //! Compute 1st derivative of cross-produced vectors
    //! \f$\vec{y}_\beta=(\vec{a}_1\times\vec{a}_2)_{,\beta}\f$
    //! with respect to nodal DOFs \f$u^D{}_L\f$
    //!
    //! HINT: The nodal DOFs can be displacements or shears.
    //!
    //! \param[in]  a0                 Vector \f$\vec{a}_1\f$
    //! \param[in]  a0DU_dl            1st derivative of \f$\vec{a}_1\f$ w.r.t.
    //!                                  nodal DOFs \f$[a^A{}_{1,D}{}^L]\f$
    //! \param[in]  a0Deriv1           Local gradient vector \f$\vec{a}_{1,\beta}\f$
    //! \param[in]  a0Deriv1DU_dl      1st derivative of local gradient \f$\vec{a}_{1,\beta}\f$ w.r.t.
    //!                                  nodal DOFs \f$[a^A{}_{1,\beta,D}{}^L]\f$
    //! \param[in]  a1                 Vector \f$\vec{a}_2\f$
    //! \param[in]  a1DU_dl            1st derivative of \f$\vec{a}_1\f$ w.r.t.
    //!                                  nodal DOFs \f$[a^A{}_{2,D}{}^L]\f$
    //! \param[in]  a1Deriv1           Local gradient vector \f$\vec{a}_{1,\beta}\f$
    //! \param[in]  a1Deriv1DU_dl      1st derivative of local gradient \f$\vec{a}_{2,\beta}\f$ w.r.t.
    //!                                  nodal DOFs \f$[a^A{}_{2,\beta,D}{}^L]\f$
    //! \param[out] pDeriv1DU_dl       1st derivative of \f$\vec{p}\f$ w.r.t.
    //!                                  nodal DOFs \f$[p^A{}_{,\beta,D}{}^L]\f$
    static
    void crossprodVectorDeriv1GradDof( const double a0[dim],
                                       const double a0DU_dl[dim],
                                       const double a0Deriv1[localDim][dim],
                                       const double a0Deriv1DU_dl[localDim][dim],
                                       const double a1[dim],
                                       const double a1DU_dl[dim],
                                       const double a1Deriv1[localDim][dim],
                                       const double a1Deriv1DU_dl[localDim][dim],
                                       double pDeriv1DU_dl[localDim][dim] );

    //! Compute 2nd derivative of cross-produced vector
    //! \f$\vec{p}=\vec{a}_1\times\vec{a}_2\f$
    //! with respect to nodal DOFs \f$u^D{}_L\f$
    //!
    //! \param[in]  a0                 Vector \f$\vec{a}_1\f$
    //! \param[in]  a0DW_dl            1st derivative of \f$\vec{a}_1\f$ w.r.t.
    //!                                  nodal DOFs \f$[a^A{}_{1,D}{}^L]\f$
    //! \param[in]  a0DU_ck            1st derivative of \f$\vec{a}_1\f$ w.r.t.
    //!                                  nodal DOFs \f$[a^A{}_{1,C}{}^K]\f$
    //! \param[in]  a0DWDU_dlck        2nd derivative of \f$\vec{a}_1\f$ w.r.t.
    //!                                  nodal DOFs \f$[a^A{}_{1,D}{}^L{}_{,C}{}^K]\f$
    //! \param[in]  a1                 Vector \f$\vec{a}_2\f$
    //! \param[in]  a1DW_dl            1st derivative of \f$\vec{a}_2\f$ w.r.t.
    //!                                  nodal DOFs \f$[a^A{}_{2,D}{}^L]\f$
    //! \param[in]  a1DU_ck            1st derivative of \f$\vec{a}_2\f$ w.r.t.
    //!                                  nodal DOFs \f$[a^A{}_{2,C}{}^K]\f$
    //! \param[in]  a1DWDU_dlck        2nd derivative of \f$\vec{a}_2\f$ w.r.t.
    //!                                  nodal DOFs \f$[a^A{}_{2,D}{}^L{}_{,C}{}^K]\f$
    //! \param[out] pDWDU_dlck         2nd derivative of \f$\vec{p}\f$ w.r.t.
    //!                                  nodal DOFs \f$[p^A{}_{,D}{}^L{}_{,C}{}^K]\f$
    static
    void crossprodVectorHessDof( const double a0[dim],
                                 const double a0DW_dl[dim],
                                 const double a0DU_ck[dim],
                                 const double a0DWDU_dlck[dim],
                                 const double a1[dim],
                                 const double a1DW_dl[dim],
                                 const double a1DU_ck[dim],
                                 const double a1DWDU_dlck[dim],
                                 double pDWDU_dlck[dim] );

    //! Compute 2nd derivative of cross-produced vectors
    //! \f$\vec{y}_\beta=(\vec{a}_1\times\vec{a}_2)_{,\beta}\f$
    //! with respect to nodal DOFs \f$u^D{}_L\f$
    //!
    //! \param[in]  a0                 Vector \f$\vec{a}_1\f$
    //! \param[in]  a0DW_dl            1st derivative of \f$\vec{a}_1\f$ w.r.t.
    //!                                  nodal DOFs \f$[a^A{}_{1,D}{}^L]\f$
    //! \param[in]  a0DU_ck            1st derivative of \f$\vec{a}_1\f$ w.r.t.
    //!                                  nodal DOFs \f$[a^A{}_{1,C}{}^K]\f$
    //! \param[in]  a0DWDU_dlck        2nd derivative of \f$\vec{a}_1\f$ w.r.t.
    //!                                  nodal DOFs \f$[a^A{}_{1,D}{}^L{}_{,C}{}^K]\f$
    //! \param[in]  a0Deriv1           Local gradient vector \f$\vec{a}_{1,\beta}\f$
    //! \param[in]  a0Deriv1DW_dl      1st derivative of local gradient \f$\vec{a}_{1,\beta}\f$ w.r.t.
    //!                                  nodal DOFs \f$[a^A{}_{1,\beta,D}{}^L]\f$
    //! \param[in]  a0Deriv1DU_ck      1st derivative of local gradient \f$\vec{a}_{1,\beta}\f$ w.r.t.
    //!                                  nodal DOFs \f$[a^A{}_{1,\beta,C}{}^K]\f$
    //! \param[in]  a0Deriv1DWDU_dlck  2nd derivative of local gradient \f$\vec{a}_{1,\beta}\f$ w.r.t.
    //!                                  nodal DOFs \f$[a^A{}_{1,\beta,D}{}^L{}_{,C}{}^K]\f$
    //! \param[in]  a1                 Vector \f$\vec{a}_2\f$
    //! \param[in]  a1DW_dl            1st derivative of \f$\vec{a}_2\f$ w.r.t.
    //!                                  nodal DOFs \f$[a^A{}_{2,D}{}^L]\f$
    //! \param[in]  a1DU_ck            1st derivative of \f$\vec{a}_2\f$ w.r.t.
    //!                                  nodal DOFs \f$[a^A{}_{2,C}{}^K]\f$
    //! \param[in]  a1DWDU_dlck        2nd derivative of \f$\vec{a}_2\f$ w.r.t.
    //!                                  nodal DOFs \f$[a^A{}_{2,D}{}^L{}_{,C}{}^K]\f$
    //! \param[in]  a1Deriv1           Local gradient vector \f$\vec{a}_{2,\beta}\f$
    //! \param[in]  a1Deriv1DW_dl      1st derivative of local gradient \f$\vec{a}_{2,\beta}\f$ w.r.t.
    //!                                  nodal DOFs \f$[a^A{}_{2,\beta,D}{}^L]\f$
    //! \param[in]  a1Deriv1DU_ck      1st derivative of local gradient \f$\vec{a}_{2,\beta}\f$ w.r.t.
    //!                                  nodal DOFs \f$[a^A{}_{2,\beta,C}{}^K]\f$
    //! \param[in]  a1Deriv1DWDU_dlck  2nd derivative of local gradient \f$\vec{a}_{2,\beta}\f$ w.r.t.
    //!                                  nodal DOFs \f$[a^A{}_{2,\beta,D}{}^L{}_{,C}{}^K]\f$
    //! \param[out] pDeriv1DWDU_dlck   2nd derivative of \f$\vec{p}_{,\beta}\f$ w.r.t.
    //!                                  nodal DOFs \f$[p^A{}_{,\beta,D}{}^L{}_{,C}{}^K]\f$
    static
    void crossprodVectorDeriv1HessDof( const double a0[dim],
                                       const double a0DW_dl[dim],
                                       const double a0DU_ck[dim],
                                       const double a0DWDU_dlck[dim],
                                       const double a0Deriv1[localDim][dim],
                                       const double a0Deriv1DW_dl[localDim][dim],
                                       const double a0Deriv1DU_ck[localDim][dim],
                                       const double a0Deriv1DWDU_dlck[localDim][dim],
                                       const double a1[dim],
                                       const double a1DW_dl[dim],
                                       const double a1DU_ck[dim],
                                       const double a1DWDU_dlck[dim],
                                       const double a1Deriv1[localDim][dim],
                                       const double a1Deriv1DW_dl[localDim][dim],
                                       const double a1Deriv1DU_ck[localDim][dim],
                                       const double a1Deriv1DWDU_dlck[localDim][dim],
                                       double pDeriv1DWDU_dlck[localDim][dim] );

    //! Compute (local) derivative of unit vector
    //! \f$\vec{d}=\vec{q}/|\vec{q}|\f$
    //!
    //! \param[in]  qLen          Length \f$|\vec{q}|\f$
    //! \param[in]  qDeriv1       Local gradient vector \f$\vec{q}_{,\beta}\f$
    //! \param[in]  d             Unit vector \f$\vec{d}\f$
    //! \param[out] dDeriv1       Local gradient vector \f$\vec{d}_{,\beta}\f$
    static
    void unitVectorDeriv1 ( const double & qLen,
                            const double qDeriv1[localDim][dim],
                            const double d[dim],
                            double dDeriv1[localDim][dim] );

    //! Compute 1st derivative of unit-vector
    //! \f$\vec{d}=\vec{q}/|\vec{q}|\f$
    //! with respect to nodal DOFs \f$u^D{}_L\f$
    //!
    //! HINT: The nodal DOFs can be displacements or shears.
    //!
    //! \param[in]  qLen   Length \f$|\vec{q}|\f$
    //! \param[in]  qDU_dl 1st derivative of \f$\vec{q}\f$ w.r.t.
    //!                      nodal DOFs \f$[q^A{}_{,D}{}^L]\f$
    //! \param[in]  d      Unit-vector \f$[d^A]\f$
    //! \param[out] dDU_dl 1st derivative of \f$\vec{d}\f$ w.r.t.
    //!                      nodal DOFs \f$[d^A{}_{,D}{}^L]\f$
    static
    void unitVectorGradDof( const double & qLen,
                            const double qDU_dl[dim],
                            const double d[dim],
                            double dDU_dl[dim] );

    //! Compute 1st derivative of local gradient of unit-vectors
    //! \f$\vec{d}_{,\beta}=\big(\vec{q}/|\vec{q}|\big)_{,\beta}\f$
    //! with respect to nodal DOFs \f$u^D{}_L\f$.
    //! Compute 1st derivative of local gradient of unit-vectors
    //! \f$\vec{d}_{,\beta}=\big(\vec{q}/|\vec{q}|\big)_{,\beta}\f$
    //! with respect to nodal DOFs \f$u^D{}_L\f$.
    //!
    //! \param[in]  qLen          Length \f$|\vec{q}|\f$
    //! \param[in]  qDU_dl        1st derivative of \f$\vec{q}\f$ w.r.t.
    //!                             nodal DOFs \f$[q^A{}_{,D}{}^L]\f$
    //! \param[in]  qDeriv1DU_dl  1st derivative of local gradient \f$\vec{q}_{,\beta}\f$ w.r.t.
    //!                             nodal DOFs \f$[q^A{}_{,\beta,D}{}^L]\f$
    //! \param[in]  d             Unit vector \f$[d^A]\f$
    //! \param[in]  dDU_dl        1st derivative of \f$\vec{d}\f$ w.r.t.
    //!                             nodal DOFs \f$[d^A{}_{,D}{}^L]\f$
    //! \param[in]  s             Vector \f$\vec{q}_{,\beta}/|\vec{q}|\f$
    //! \param[out] cDU_dl        1st derivative of \f$\vec{d}_{,\beta}\f$ w.r.t.
    //!                             nodal DOFs \f$[d^A{}_{,\beta,D}{}^L]\f$
    static
    void unitVectorDeriv1GradDof( const double & qLen,
                                  const double qDU_dl[dim],
                                  const double qDeriv1DU_dl[localDim][dim],
                                  const double d[dim],
                                  const double dDU_dl[dim],
                                  const double s[localDim][dim],
                                  double cDU_dl[localDim][dim] );

    //! Compute 2nd derivative of unit-vector
    //! \f$\vec{d}=\vec{q}/|\vec{q}|\f$
    //! with respect to nodal DOFs
    //!
    //! \param[in]  qLen        Length \f$|\vec{q}|\f$
    //! \param[in]  qDW_dl      1st deriv. of \f$\vec{q}\f$ w.r.t. DOFs
    //!                           \f$[q^A{}_{,D}{}^L]\f$
    //! \param[in]  qDU_ck      1st deriv. of \f$\vec{q}\f$ w.r.t. DOFs
    //!                           \f$[q^A{}_{,C}{}^K]\f$
    //! \param[in]  qDWDU_dlck  2nd deriv. of \f$\vec{q}\f$ w.r.t. DOFs
    //!                           \f$[q^A{}_{,D}{}^L{}_{,C}{}^K]\f$
    //! \param[in]  d           Unit vector \f$[d^A]\f$
    //!                           \f$[d^A]\f$
    //! \param[in]  dDW_dl      1st deriv. of director w.r.t. DOFs
    //!                           \f$[d^A{}_{,D}{}^L]\f$
    //! \param[in]  dDU_ck      1st deriv. of director w.r.t. DOFs
    //!                           \f$[d^A{}_{,C}{}^K]\f$
    //! \param[out] dDWDU_dlck  2nd deriv. of director w.r.t. DOFs
    //!                           \f$[d^A{}_{,D}{}^L{}_{,C}{}^K]\f$
    static
    void unitVectorHessDof( const double & qLen,
                            const double qDW_dl[dim],
                            const double qDU_ck[dim],
                            const double qDWDU_dlck[dim],
                            const double d[dim],
                            const double dDW_dl[dim],
                            const double dDU_ck[dim],
                            double dDWDU_dlck[dim] );

    //! Compute 2nd derivative of local gradient of unit-vectors
    //! \f$\vec{d}_{,\beta}=\big(\vec{q}/|\vec{q}|\big)_{,\beta}\f$
    //! with respect to nodal DOFs \f$u^D{}_L\f$.
    //! 
    //! \param[in]  qLen             Length \f$|\vec{q}|\f$
    //! \param[in]  qDW_dl           1st derivative of \f$\vec{q}\f$ w.r.t.
    //!                                nodal DOFs \f$[q^A{}_{,D}{}^L]\f$
    //! \param[in]  qDU_ck           1st derivative of \f$\vec{q}\f$ w.r.t.
    //!                                nodal DOFs \f$[q^A{}_{,C}{}^K]\f$
    //! \param[in]  qDWDU_dlck       2nd deriv. of \f$\vec{q}\f$ w.r.t. DOFs
    //!                                \f$[q^A{}_{,D}{}^L{}_{,C}{}^K]\f$
    //! \param[in]  qDeriv1DW_dl     1st derivative of local gradient \f$\vec{q}_{,\beta}\f$
    //!                                w.r.t. DOFs \f$[q^A{}_{,\beta,D}{}^L]\f$
    //! \param[in]  qDeriv1DU_ck     1st derivative of local gradient \f$\vec{q}_{,\beta}\f$
    //!                                w.r.t. DOFs \f$[q^A{}_{,\beta,C}{}^K]\f$
    //! \param[in]  qDeriv1DWDU_dlck  2nd deriv. of \f$\vec{q}_{,\beta}\f$ w.r.t. DOFs
    //!                                \f$[q^A{}_{,\beta,D}{}^L{}_{,C}{}^K]\f$
    //! \param[in]  d                Unit-vector \f$[d^A]\f$
    //! \param[in]  dDW_dl           1st deriv. of unit-vector w.r.t. DOFs
    //!                                \f$[d^A{}_{,D}{}^L]\f$
    //! \param[in]  dDU_ck           1st deriv. of unit-vector w.r.t. DOFs
    //!                                \f$[d^A{}_{,C}{}^K]\f$
    //! \param[in]  dDWDU_dlck       2nd deriv. of unit-vector w.r.t. DOFs
    //!                                \f$[d^A{}_{,D}{}^L{}_{,C}{}^K]\f$
    //! \param[in]  s                Vector \f$\vec{q}_{,\beta}/|\vec{q}|\f$
    //! \param[out] cDWDU_dlck       2nd deriv. of unit-vector w.r.t. DOFs
    //!                                \f$[d^A{}_{,\beta,D}{}^L{}_{,C}{}^K]\f$
    static
    void unitVectorDeriv1HessDof( const double & qLen,
                                  const double qDW_dl[dim],
                                  const double qDU_ck[dim],
                                  const double qDWDU_dlck[dim],
                                  const double qDeriv1DW_dl[localDim][dim],
                                  const double qDeriv1DU_ck[localDim][dim],
                                  const double qDeriv1DWDU_dlck[localDim][dim],
                                  const double d[dim],
                                  const double dDW_dl[dim],
                                  const double dDU_ck[dim],
                                  const double dDWDU_dlck[dim],
                                  const double s[localDim][dim],
                                  double cDWDU_dlck[localDim][dim] );

private:
    //! Constructor is private
    SurfaceGeometer( ) { }

    //! Copy constructor is private
    SurfaceGeometer( SurfaceGeometer & other ) { }
};

//------------------------------------------------------------------------------
#include "SurfaceGeometer.ipp"
//------------------------------------------------------------------------------
#endif
