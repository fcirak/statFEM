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

//------------------------------------------------------------------------------
// Map to get Voigt indices from pair of local indices
// \f$(\alpha,\beta) \mapsto \mathcal{A}\f$
//
// Tensor indices     \alpha\beta=  11,  12,  21,  22
// C indices                        00,  01,  10,  11
// Access   2*\alpha+\beta           0,   1    2    3
// 3-Voigt C-indices                 0    1    1    2
const unsigned gshell::fem::SurfaceGeometer::voigtForward[2][2] = { { 0, 1 }, { 1, 2 } };

//------------------------------------------------------------------------------
// Interpolate vector quantity \f$q^A\f$
void gshell::fem::SurfaceGeometer::interpolate(
    const unsigned numFunctions,
    const double * phi,  // (NF)-vector
    const double * qNF,  // (dim)x(NF)-array
    double q[dim]
    )
{
    la::productN<dim>( q, numFunctions, qNF, phi );
    return;
}

//------------------------------------------------------------------------------
// Interpolate first derivative of vector quantity w.r.t. local co-ordinates
void gshell::fem::SurfaceGeometer::interpolateDeriv1(
    const unsigned numFunctions,  // NF
    const double * dPhiDXi,  // (localDim)x(NF)-array
    const double * qNF,  // (dim)x(NF)-array
    double qDeriv1[localDim][dim]
    )
{
    la::productNT<localDim,dim>( qDeriv1, numFunctions, dPhiDXi, qNF );
    return;
}

//------------------------------------------------------------------------------
// Interpolate second derivative of vector quantity w.r.t. local co-ordinates
void gshell::fem::SurfaceGeometer::interpolateDeriv2(
    const unsigned numFunctions,  // NF
    const double * ddPhiDDXi,  // (sDim)x(NF)-matrix
    const double * qNF,  // (dim)x(NF)-matrix
    double qDeriv2[localDim][localDim][dim]
    )
{
    // \vec{q}_{,\alpha\beta} in Voigt-matrix
    double ddqDDXi[ sDim ][dim];
    la::productNT<sDim,dim>( ddqDDXi, numFunctions, ddPhiDDXi, qNF );
    // \vec{q}_{,\alpha\beta} in index notation
    for ( unsigned al=0; al<localDim; ++al )
        for ( unsigned be=0; be<localDim; ++be )
            la::assign<dim>( qDeriv2[ al ][ be ], ddqDDXi[ voigtForward[ al ][ be ] ] );
    return;
}

//------------------------------------------------------------------------------
// Computes the co-variant base vectors in mid-surface
double gshell::fem::SurfaceGeometer::covariantBase(
    const unsigned numFunctions,  // NF
    const double * dPhiDXi,  // (localDim)x(NF)-matrix
    const double * xNF,  // (dim)x(NF)-matrix
    double a[dim][dim]
    )
{
    // covariant tangent vectors
    la::productNT<localDim,dim>( a, numFunctions, dPhiDXi, xNF );
    // normal vector
    la::crossProduct( a[ 2 ], a[ 0 ], a[ 1 ] );
    // normalise 
    double len;  la::norm2<dim>( len, a[ 2 ] );
    la::divide<dim>( a[ 2 ], len );
    // return its determinant
    return la::determinant<dim>( a );
}

//------------------------------------------------------------------------------
// Computes the derivatives of co-variant base vectors in mid-surface
void gshell::fem::SurfaceGeometer::covariantBaseDeriv1(
    const unsigned numFunctions,  // NF
    const double * ddPhiDDXi,  // (sDim)x(NF)-matrix
    const double * xNF,  // (dim)x(NF)-matrix
    const double a[dim][dim],
    double b[dim][localDim][dim]
    )
{
    // second derivatives of position w.r.t. local co-ordinates
    //
    // \vec{x}_{,\alpha\beta}
    //
    // = [ \vec{x}_{,11} \vec{x}_{,12} \vec{x}_{,22} ]
    //
    //   [       x_{,11}       x_{,12}       x_{,22} ]
    // = [       y_{,11}       y_{,12}       y_{,22} ]
    //   [       z_{,11}       z_{,12}       z_{,22} ]
    double ddxDDXi[ sDim ][dim];
    la::productNT<sDim,dim>( ddxDDXi, numFunctions, ddPhiDDXi, xNF );
    
    // store \vec{a}_{\alpha,\beta}
    for ( unsigned al=0; al<localDim; ++al ) 
        for ( unsigned be=0; be<localDim; ++be )
            la::assign<dim>( b[ al ][ be ], ddxDDXi[ voigtForward[ al ][ be ] ] );

    // compute \vec{a}_{3,\beta}
#if 1
    double p[dim];  la::crossProduct( p, a[ 0 ], a[ 1 ] );
    double pLen;    la::norm2<dim>( pLen, p );
    double pDeriv1[localDim][dim];
    SurfaceGeometer::crossprodVectorDeriv1( a[ 0 ], b[ 0 ], a[ 1 ], b[ 1 ], pDeriv1 );
    SurfaceGeometer::unitVectorDeriv1( pLen, pDeriv1, a[ 2 ], b[ 2 ] );
#else
    double tmp[dim];
    la::crossProduct( tmp, a[ 0 ], a[ 1 ] );
    double len;
    la::norm2<dim>( len, tmp );
    for ( unsigned al=0; al<localDim; ++al ) {
        la::crossProduct( tmp, b[ 0 ][ al ], a[ 1 ] );
        la::addCrossProduct( tmp, a[ 0 ], b[ 1 ][ al ] );
        la::divide<dim>( tmp, len );
        double fact;
        la::innerProduct<dim>( fact, a[ 2 ], tmp );
        la::update<dim>( tmp, -fact, a[ 2 ] );
        la::assign<dim>( b[ 2 ][ al ], tmp );
    }
#endif

    return;
}

//------------------------------------------------------------------------------
// Computes the thickness director
void gshell::fem::SurfaceGeometer::director(
    const double p[dim],
    const double w[dim],
    double d[dim]
    )
{
    la::assign<dim>( d, p );
    la::add<dim>( d, w );
    double len;
    la::norm2<dim>( len, d );
    la::divide<dim>( d, len );
    return;
}

//------------------------------------------------------------------------------
// Computes the local gradient of the thickness director
void gshell::fem::SurfaceGeometer::directorDeriv1(
    const double a[dim][dim],
    const double w[dim],
    const double d[dim],
    const double b[dim][localDim][dim],
    const double z[localDim][dim],
    double c[localDim][dim]
    )
{
    // local gradient of thickness director
    double tmp[dim];
    la::assign<dim>( tmp, a[ 2 ] );
    la::add<dim>( tmp, w );
    double len;  la::norm2<dim>( len, tmp );
    for ( unsigned al=0; al<localDim; ++al ) {
        la::assign<dim>( tmp, b[ 2 ][ al ] );
        la::add<dim>( tmp, z[ al ] );
        la::divide<dim>( tmp, len );
        double fact;
        la::innerProduct<dim>( fact, w, b[ 2 ][ al ] );
        fact /= len;
        la::addInnerProduct<dim>( fact, d, z[ al ] );
        fact /= len;
        la::update<dim>( tmp, -fact, d );
        la::assign<dim>( c[ al ], tmp );
    }

    // done
    return;
}

//------------------------------------------------------------------------------
// Compute (local)  derivative of cross-produced vector
// \f$\vec{p}=\vec{a}_1\times\vec{a}_2\f$
void gshell::fem::SurfaceGeometer::crossprodVectorDeriv1(
    const double a0[dim],
    const double a0Deriv1[localDim][dim],
    const double a1[dim],
    const double a1Deriv1[localDim][dim],
    double pDeriv1[localDim][dim]
    )
{
    for ( unsigned be=0; be<localDim; ++be ) {
        la::crossProduct(    pDeriv1[ be ], a0Deriv1[ be ], a1 );
        la::addCrossProduct( pDeriv1[ be ], a0, a1Deriv1[ be ] );
    }
    return;
}

//------------------------------------------------------------------------------
// Compute 1st derivative of cross-produced vector
// \f$\vec{p}=\vec{a}_1\times\vec{a}_2\f$
// with respect to nodal DOFs \f$u^D{}_L\f$
void gshell::fem::SurfaceGeometer::crossprodVectorGradDof(
    const double a0[dim],
    const double a0DU_dl[dim],
    const double a1[dim],
    const double a1DU_dl[dim],
    double pDU_dl[dim]
    )
{
    la::crossProduct(    pDU_dl, a0DU_dl, a1 );
    la::addCrossProduct( pDU_dl, a0, a1DU_dl );
    return;
}

//------------------------------------------------------------------------------
// Compute 1st derivative of cross-produced vectors
// \f$\vec{y}_\beta=(\vec{a}_1\times\vec{a}_2)_{,\beta}\f$
// with respect to nodal DOFs \f$u^D{}_L\f$
void gshell::fem::SurfaceGeometer::crossprodVectorDeriv1GradDof(
    const double a0[dim],
    const double a0DU_dl[dim],
    const double a0Deriv1[localDim][dim],
    const double a0Deriv1DU_dl[localDim][dim],
    const double a1[dim],
    const double a1DU_dl[dim],
    const double a1Deriv1[localDim][dim],
    const double a1Deriv1DU_dl[localDim][dim],
    double pDeriv1DU_dl[localDim][dim]
    )
{
    for ( unsigned be=0; be<localDim; ++be ) {
        la::crossProduct(    pDeriv1DU_dl[ be ], a0Deriv1DU_dl[ be ], a1 );
        la::addCrossProduct( pDeriv1DU_dl[ be ], a0Deriv1[ be ], a1DU_dl );
        la::addCrossProduct( pDeriv1DU_dl[ be ], a0DU_dl, a1Deriv1[ be ] );
        la::addCrossProduct( pDeriv1DU_dl[ be ], a0, a1Deriv1DU_dl[ be ] );
    }
    return;
}
//------------------------------------------------------------------------------
// Compute 2nd derivative of cross-produced vector
// \f$\vec{p}=\vec{a}_1\times\vec{a}_2\f$
// with respect to nodal DOFs \f$u^D{}_L\f$
void gshell::fem::SurfaceGeometer::crossprodVectorHessDof(
    const double a0[dim],
    const double a0DW_dl[dim],
    const double a0DU_ck[dim],
    const double a0DWDU_dlck[dim],
    const double a1[dim],
    const double a1DW_dl[dim],
    const double a1DU_ck[dim],
    const double a1DWDU_dlck[dim],
    double pDWDU_dlck[dim]
    )
{
    la::crossProduct(    pDWDU_dlck, a0DWDU_dlck, a1 );
    la::addCrossProduct( pDWDU_dlck, a0DW_dl, a1DU_ck );
    la::addCrossProduct( pDWDU_dlck, a0DU_ck, a1DW_dl );
    la::addCrossProduct( pDWDU_dlck, a0, a1DWDU_dlck );
    return;
}

//------------------------------------------------------------------------------
// Compute 2nd derivative of cross-produced vectors
// \f$\vec{y}_\beta=(\vec{a}_1\times\vec{a}_2)_{,\beta}\f$
// with respect to nodal DOFs \f$u^D{}_L\f$
void gshell::fem::SurfaceGeometer::crossprodVectorDeriv1HessDof(
    const double a0[dim],
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
    double pDeriv1DWDU_dlck[localDim][dim]
    )
{
    for ( unsigned be=0; be<localDim; ++be ) {
        la::crossProduct(    pDeriv1DWDU_dlck[ be ], a0Deriv1DWDU_dlck[ be ], a1 );
        la::addCrossProduct( pDeriv1DWDU_dlck[ be ], a0Deriv1DW_dl[ be ], a1DU_ck );
        la::addCrossProduct( pDeriv1DWDU_dlck[ be ], a0Deriv1DU_ck[ be ], a1DW_dl );
        la::addCrossProduct( pDeriv1DWDU_dlck[ be ], a0Deriv1[ be ], a1DWDU_dlck );
        la::addCrossProduct( pDeriv1DWDU_dlck[ be ], a0DWDU_dlck, a1Deriv1[ be ] );
        la::addCrossProduct( pDeriv1DWDU_dlck[ be ], a0DW_dl, a1Deriv1DU_ck[ be ] );
        la::addCrossProduct( pDeriv1DWDU_dlck[ be ], a0DU_ck, a1Deriv1DW_dl[ be ] );
        la::addCrossProduct( pDeriv1DWDU_dlck[ be ], a0, a1Deriv1DWDU_dlck[ be ] );
    }
    return;
}

//------------------------------------------------------------------------------
// Compute (local) derivative of unit vector
// \f$\vec{d}=\vec{q}/|\vec{q}|\f$
void gshell::fem::SurfaceGeometer::unitVectorDeriv1(
    const double & qLen,
    const double qDeriv1[localDim][dim],
    const double d[dim],
    double dDeriv1[localDim][dim]
    )
{
    la::assign<localDim,dim>( dDeriv1, qDeriv1 );
    for ( unsigned al=0; al<localDim; ++al ) {
        double fact_al;
        la::innerProduct<dim>( fact_al, d, qDeriv1[ al ] );
        la::update<dim>( dDeriv1[ al ], -fact_al, d );
    }
    la::divide<localDim,dim>( dDeriv1, qLen );
    return;
}

//------------------------------------------------------------------------------
// Compute 1st derivative of unit-vector
// \f$\vec{d}=\vec{q}/|\vec{q}|\f$
// with respect to nodal DOFs \f$u^D{}_L\f$
void gshell::fem::SurfaceGeometer::unitVectorGradDof(
    const double & qLen,
    const double qDU_dl[dim],
    const double d[dim],
    double dDU_dl[dim]
    )
{
    double fact_dl;
    la::innerProduct<dim>( fact_dl, d, qDU_dl );
    la::assign<dim>( dDU_dl, qDU_dl );
    la::update<dim>( dDU_dl, -fact_dl, d );
    la::multiply<dim>(  dDU_dl, 1./qLen );
    return;
}

//------------------------------------------------------------------------------
// Compute 1st derivative of local gradient of unit-vectors
// \f$\vec{d}_{,\beta}=\big(\vec{q}/|\vec{q}|\big)_{,\beta}\f$
// with respect to nodal DOFs \f$u^D{}_L\f$.
void gshell::fem::SurfaceGeometer::unitVectorDeriv1GradDof(
    const double & qLen,
    const double qDU_dl[dim],
    const double qDeriv1DU_dl[localDim][dim],
    const double d[dim],
    const double dDU_dl[dim],
    const double s[localDim][dim],
    double cDU_dl[localDim][dim]
    )
{
    double fact_dl;
    la::innerProduct<dim>( fact_dl, d, qDU_dl );

    double sDU_dl[localDim][dim];
    la::assign<localDim,dim>( sDU_dl, qDeriv1DU_dl );
    la::update<localDim,dim>( sDU_dl, -fact_dl, s );
    la::divide<localDim,dim>( sDU_dl, qLen );

    for ( unsigned be=0; be<localDim; ++be ) {
        
        la::innerProduct<dim>(    fact_dl, s[ be ], dDU_dl );
        la::addInnerProduct<dim>( fact_dl, d, sDU_dl[ be ] );

        double fact;
        la::innerProduct<dim>( fact, d, s[ be ] );

        la::assign<dim>( cDU_dl[ be ], sDU_dl[ be ] );
        la::update<dim>( cDU_dl[ be ], -fact_dl, d );
        la::update<dim>( cDU_dl[ be ], -fact, dDU_dl );
    }

    return;
}

//------------------------------------------------------------------------------
// Compute 2nd derivative of unit-vector
// \f$\vec{d}=\vec{q}/|\vec{q}|\f$
// with respect to nodal DOFs
void gshell::fem::SurfaceGeometer::unitVectorHessDof(
    const double & qLen,
    const double qDW_dl[dim],
    const double qDU_ck[dim],
    const double qDWDU_dlck[dim],
    const double d[dim],
    const double dDW_dl[dim],
    const double dDU_ck[dim],
    double dDWDU_dlck[dim]
    )
{
    double fact_ck;
    la::innerProduct<dim>( fact_ck, d, qDU_ck );
    double fact_dl;
    la::innerProduct<dim>( fact_dl, d, qDW_dl );
    double fact_dlck;
    la::innerProduct<dim>( fact_dlck, dDU_ck, qDW_dl );
    la::addInnerProduct<dim>( fact_dlck, d, qDWDU_dlck );

    la::assign<dim>( dDWDU_dlck, qDWDU_dlck );
    la::update<dim>( dDWDU_dlck, -fact_ck, dDW_dl );
    la::update<dim>( dDWDU_dlck, -fact_dl, dDU_ck );
    la::update<dim>( dDWDU_dlck, -fact_dlck, d );
    la::divide<dim>( dDWDU_dlck, qLen );
    return;
}

//------------------------------------------------------------------------------
// Compute 2nd derivative of local gradient of unit-vectors
// \f$\vec{d}_{,\beta}=\big(\vec{q}/|\vec{q}|\big)_{,\beta}\f$
// with respect to nodal DOFs \f$u^D{}_L\f$.
void gshell::fem::SurfaceGeometer::unitVectorDeriv1HessDof(
    const double & qLen,
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
    double cDWDU_dlck[localDim][dim]
    )
{
    double fact_dl;
    la::innerProduct<dim>( fact_dl, d, qDW_dl );
    double fact_ck;
    la::innerProduct<dim>( fact_ck, d, qDU_ck );
    double fact_dlck;
    la::innerProduct<dim>( fact_dlck, dDU_ck, qDW_dl );
    la::addInnerProduct<dim>( fact_dlck,d, qDWDU_dlck );

    double sDW_dl[localDim][dim];
    la::assign<localDim,dim>( sDW_dl, qDeriv1DW_dl );
    la::update<localDim,dim>( sDW_dl, -fact_dl, s );
    la::divide<localDim,dim>( sDW_dl, qLen );

    double sDU_ck[localDim][dim];
    la::assign<localDim,dim>( sDU_ck, qDeriv1DU_ck );
    la::update<localDim,dim>( sDU_ck, -fact_ck, s );
    la::divide<localDim,dim>( sDU_ck, qLen );
    
    double sDWDU_dlck[localDim][dim];
    la::assign<localDim,dim>( sDWDU_dlck, qDeriv1DWDU_dlck );
    la::update<localDim,dim>( sDWDU_dlck, -fact_ck, sDW_dl );
    la::update<localDim,dim>( sDWDU_dlck, -fact_dl, sDU_ck );
    la::update<localDim,dim>( sDWDU_dlck, -fact_dlck, s );
    la::divide<localDim,dim>( sDWDU_dlck, qLen );

    for ( unsigned be=0; be<localDim; ++be ) {

        la::innerProduct<dim>( fact_dl, s[ be ], dDW_dl );
        la::addInnerProduct<dim>( fact_dl, d, sDW_dl[ be ] );

        la::innerProduct<dim>( fact_dlck, sDU_ck[ be ], dDW_dl );
        la::addInnerProduct<dim>( fact_dlck, s[ be ], dDWDU_dlck );
        la::addInnerProduct<dim>( fact_dlck, dDU_ck, sDW_dl[ be ] );
        la::addInnerProduct<dim>( fact_dlck, d, sDWDU_dlck[ be ] );

        la::innerProduct<dim>( fact_ck, dDU_ck, s[ be ] );
        la::addInnerProduct<dim>( fact_ck, d, sDU_ck[ be ] );

        double fact;
        la::innerProduct<dim>( fact, d, s[ be ] );

        la::assign<dim>( cDWDU_dlck[ be ], sDWDU_dlck[ be ] );
        la::update<dim>( cDWDU_dlck[ be ], -fact_dl, dDU_ck );
        la::update<dim>( cDWDU_dlck[ be ], -fact_dlck, d );
        la::update<dim>( cDWDU_dlck[ be ], -fact_ck, dDW_dl );
        la::update<dim>( cDWDU_dlck[ be ], -fact, dDWDU_dlck );

    }

    return;
}

