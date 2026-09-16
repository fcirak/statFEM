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
// Tensor indices     \alpha\beta=  11,  12,  21,  22
// C indices                        00,  01,  10,  11
// Access   2*\alpha+\beta           0,   1    2    3
// 3-Voigt C-indices                 0    1    1    2
template< typename BELEMENT, typename MAT >
const unsigned gshell::fem::ElementKLStatic<BELEMENT,MAT>::
voigtForward[localDim][localDim] = {
    { BELEMENT::voigtForward[0][0], BELEMENT::voigtForward[0][1] },
    { BELEMENT::voigtForward[1][0], BELEMENT::voigtForward[1][1] }
};

//------------------------------------------------------------------------------
// Accessor to the indices of the nodal degrees of freedom
template< typename BELEMENT, typename MAT >
void gshell::fem::ElementKLStatic<BELEMENT,MAT>::getDofIndices(
    std::vector< unsigned > & dofIndices
    ) const
{
    const unsigned numNodes = this->BasisElement::numFunctions( );
    const unsigned matSize = dof * numNodes;
    dofIndices.resize( matSize );
    std::vector< unsigned >::iterator iter = dofIndices.begin();
    for ( unsigned v = 0; v < numNodes; v ++ ) { 
        std::vector< unsigned > aux;
        this->BasisElement::giveSupportNodePtr( v ) -> copyDofArray( aux );
        iter = std::copy( aux.begin(), aux.end(), iter );
    }
    return;
}

//------------------------------------------------------------------------------
// Collect general nodal quantitiy accessed by the passed operator op
template< typename BELEMENT, typename MAT >
template< typename OP >
void gshell::fem::ElementKLStatic<BELEMENT,MAT>::nodalQuantity_(
    MatDimNF & quanU,
    OP op
    ) const
{
    const unsigned numNodes = this->BasisElement::numFunctions( );
    quanU.conservativeResize( dim, numNodes );
    for ( unsigned j = 0; j < numNodes; j ++ ) {
        quanU.col( j ) = op( this->BasisElement::giveSupportNodePtr( j ) );
    }
    return;
}

//------------------------------------------------------------------------------
// Collect general nodal quantitiy accessed by the passed operator op
template< typename BELEMENT, typename MAT >
template< typename OP >
void gshell::fem::ElementKLStatic<BELEMENT,MAT>::updateNodalQuantity_(
    MatDimNF & quanU,
    const double fact,
    OP op
    ) const
{
    const unsigned numNodes = this->BasisElement::numFunctions( );
    quanU.conservativeResize( dim, numNodes );
    for ( unsigned j = 0; j < numNodes; j ++ ) {
        const VecDim q = op( this->BasisElement::giveSupportNodePtr( j ) );
        for ( unsigned d = 0; d < dim; ++d ) {
            quanU( d, j ) += fact * q[ d ];
        }
    }
    return;
}

//------------------------------------------------------------------------------
// Return nodal displacements
template< typename BELEMENT, typename MAT >
void gshell::fem::ElementKLStatic<BELEMENT,MAT>::nodalDisplacements(
    MatDimNF & u
    ) const
{
    this->nodalQuantity_( u, std::mem_fun( &Node::giveDisplacements ) );

    const unsigned nf = this->BasisElement::numFunctions( );
    FTL_VERIFY_DESCRIPTIVE( u.rows( ) == dim,
                            "u.rows( )=%d unequal dim=%d\n", u.rows( ), dim );
    FTL_VERIFY_DESCRIPTIVE( u.cols( ) == nf,
                            "u.cols( )=%d unequal nf=%d\n",  u.cols( ), nf  );
    return;
}

//------------------------------------------------------------------------------
// Return nodal displacement increments
template< typename BELEMENT, typename MAT >
void gshell::fem::ElementKLStatic<BELEMENT,MAT>::nodalIncrements(
    MatDimNF & deltaU
    ) const
{
    this->nodalQuantity_( deltaU, std::mem_fun( &Node::giveIncrement ) );
    return;
}

//------------------------------------------------------------------------------
// Interpolate current displacement at local point
template< typename BELEMENT, typename MAT >
corlib::eigenX::VectorSd< BELEMENT::dim > 
gshell::fem::ElementKLStatic<BELEMENT,MAT>::giveDisplacement(
    const VecLDim & xi
    ) const
{
    // check xi
    FTL_VERIFY_DESCRIPTIVE( corlib::ShapeTraits< myShape >::isInside( xi ),
                            "with xi[0]=%g, xi[1]=%g\n",
                            xi[0], xi[1] );
    // shape functions at xi
    VecNF phi;          this->BasisElement::sfun( xi, phi );
    // nodal displacements (and shears)
    MatDimNF uNF;       this->nodalDisplacements( uNF );
    // interpolate
    const MatDimNF uNFT = uNF.transpose();
    VecDim u;           SurfaceGeometer::interpolate( this->BasisElement::numFunctions(), &(phi(0)),
                                                      &(uNFT(0,0)), &(u(0)) );
    return u;
}

//------------------------------------------------------------------------------
// Interpolate current displacement at local point
template< typename BELEMENT, typename MAT >
corlib::eigenX::VectorSd< BELEMENT::dim > 
gshell::fem::ElementKLStatic<BELEMENT,MAT>::giveScaledNormal(
    const VecLDim & xi
    ) const
{
    // nodal coordinates
    MatDimNF xNF;      this -> supportNodeCoordinates( xNF );
    // nodal displacements
    MatDimNF uNF;      this -> nodalDisplacements( uNF );
    // shape function derivatives
    MatLDimNF dPhiDXi; this -> sfunGrad( xi, dPhiDXi );

    // reference tangents
    const MatDimNF xNFT     = xNF.transpose();
    const MatDimNF dPhiDXiT = dPhiDXi.transpose();
    VecDim tang0Ref;   SurfaceGeometer::interpolate( this->BasisElement::numFunctions(),
                                                     &(dPhiDXiT(0,0)), &(xNFT(0,0)),
                                                     &(tang0Ref(0)) );
    VecDim tang1Ref;   SurfaceGeometer::interpolate( this->BasisElement::numFunctions(),
                                                     &(dPhiDXiT(0,1)), &(xNFT(0,0)),
                                                     &(tang1Ref(0)) );
    // current tangents
    xNF += uNF;
    MatDimNF xNFTcur = xNF.transpose();
    VecDim tang0Cur;   SurfaceGeometer::interpolate( this->BasisElement::numFunctions(),
                                                     &(dPhiDXiT(0,0)), &(xNFTcur(0,0)),
                                                     &(tang0Cur(0)) );
    VecDim tang1Cur;   SurfaceGeometer::interpolate( this->BasisElement::numFunctions(),
                                                     &(dPhiDXiT(0,1)), &(xNFTcur(0,0)),
                                                     &(tang1Cur(0)) );


    // normal
    VecDim normalCur;  la::crossProduct( &(normalCur(0)), &(tang0Cur(0)), &(tang1Cur(0)) );
    VecDim normalRef;  la::crossProduct( &(normalRef(0)), &(tang0Ref(0)), &(tang1Ref(0)) );
    double metricRef;  la::norm2<dim>( metricRef, &(normalRef(0)) );
    la::divide<dim>( &(normalCur(0)), metricRef*thickness_ );
    return normalCur;
}

//------------------------------------------------------------------------------
// Interpolate current displacement at local point
template< typename BELEMENT, typename MAT >
corlib::eigenX::VectorSd< BELEMENT::dim > 
gshell::fem::ElementKLStatic<BELEMENT,MAT>::giveDisplacementAtNode(
    const unsigned index
    ) const
{
    // check xi
    FTL_VERIFY( index < numVertices );
    // shape functions at xi
    VecNF phi;          this->BasisElement::sfun( index, phi );
    // nodal displacements (and shears)
    MatDimNF uNF;       this->nodalDisplacements( uNF );
    // interpolate
    const MatDimNF uNFT = uNF.transpose();
    VecDim u;           SurfaceGeometer::interpolate( this->numFunctions(), &(phi(0)),
                                                      &(uNFT(0,0)), &(u(0)) );
    return u;
}

//------------------------------------------------------------------------------
// Compute tangents to limit surface at node
template< typename BELEMENT, typename MAT >
void gshell::fem::ElementKLStatic<BELEMENT,MAT>::computeTangentsAtNode(
    const enum gshell::fem::config conf,
    const unsigned ind,
    VecDim & tang0,
    VecDim & tang1,
    Mat * tang0GradDof,
    Mat * tang1GradDof
    )
{
    // number of shape functions
    const unsigned numFunctions = this->numFunctions();

    // nodal reference co-ordinates
    MatDimNF xNF;  this->supportNodeCoordinates( xNF );
    // nodal displacements (and shears)
    if ( conf == gshell::fem::CURRENT ) {
        MatDimNF uNF;  this->nodalDisplacements( uNF );
        xNF += uNF;
    }
        
    // get coefficients
    FTL_VERIFY( ind < numVertices );
    MatLDimNF coeff( localDim, numFunctions );
    this->sfunGrad( ind, coeff );
        
    // compute tangents at configuration
    typedef eigenX::MatrixSd< dim, localDim > MatDimLDim;
    const MatDimLDim tang = xNF * coeff.transpose( );
    tang0 = tang.col( 0 );
    tang1 = tang.col( 1 );
        
    // return gradients
    if ( tang0GradDof or tang1GradDof ) {
        if ( tang0GradDof ) { 
            tang0GradDof->resize( dim, numFunctions*dof );
            tang0GradDof->setZero();
        }
        if ( tang1GradDof ) {
            tang1GradDof->resize( dim, numFunctions*dof );
            tang1GradDof->setZero();
        }
        for ( unsigned i=0; i<numFunctions; ++i ) {
            for ( unsigned d=0; d<dim; ++d ) {
                if ( tang0GradDof ) (*tang0GradDof)( d, i*dof+d ) = coeff( 0, i );
                if ( tang1GradDof ) (*tang1GradDof)( d, i*dof+d ) = coeff( 1, i );
            }
        }
    }
}

//------------------------------------------------------------------------------
// Compute normal to limit surface at node
template< typename BELEMENT, typename MAT >
void gshell::fem::ElementKLStatic<BELEMENT,MAT>::computeNormalAtNode(
    const enum gshell::fem::config conf,
    const unsigned ind,
    VecDim & normal,
    Mat * normalGradDof
    )
{
    // number of shape fucntions
    const unsigned numFunctions = this->numFunctions();

    // compute non-unit tangents
    VecDim tang0, tang1;
    Mat * tang0GradDof = ( normalGradDof ) ? new Mat( dim, numFunctions*dof ) : NULL;
    Mat * tang1GradDof = ( normalGradDof ) ? new Mat( dim, numFunctions*dof ) : NULL;
    this->computeTangentsAtNode( conf, ind,
                                 tang0, tang1, tang0GradDof, tang1GradDof );

    // compute non-unit normal
    const VecDim z    = corlib::cross_prod( tang0, tang1 );
    const double zLen = z.norm( );

    // compute unit normal
    FTL_VERIFY( not corlib::fuzzyEqual( zLen, 0., 1.e-10 ) );
    normal = z / zLen;

    // compute gradient
    if ( normalGradDof ) {
        normalGradDof->resize( dim, numFunctions*dof );
        for ( unsigned l=0; l<numFunctions; ++l ) {
            for ( unsigned d=0; d<dof; ++d ) {
                const VecDim tang0DU_dl = ( *tang0GradDof ).col( l*dof+d );
                const VecDim tang1DU_dl = ( *tang1GradDof ).col( l*dof+d );

                // gradient of non-unit normal
                double zDU_dl[dim];
                SurfaceGeometer::crossprodVectorGradDof( &(tang0(0)), &(tang0DU_dl(0)),
                                                         &(tang1(0)), &(tang1DU_dl(0)),
                                                         zDU_dl );

                // gradient of unit normal-
                VecDim nDU_dl;
                SurfaceGeometer::unitVectorGradDof( zLen, zDU_dl, &(normal(0)), &(nDU_dl(0)) );
                // store latter
                ( *normalGradDof ).col( l*dof+d ) = nDU_dl;
            }
        }
    }

    // remove allocated matrices
    if ( tang0GradDof ) delete tang0GradDof;
    if ( tang1GradDof ) delete tang1GradDof;
}

//------------------------------------------------------------------------------
// Resize and (rows x NF)-matrix of ENTRYs and optionally set to zero
template< typename BELEMENT, typename MAT >
template< typename ENTRY >
void gshell::fem::ElementKLStatic<BELEMENT,MAT>::resizeMatrixNF_(
    Eigen::Matrix<ENTRY, Eigen::Dynamic, Eigen::Dynamic> & mat,
    unsigned rows,
    bool setToZero
    ) const
{
    const unsigned numFunctions = this->BasisElement::numFunctions( );
    mat.resize( rows, numFunctions );
    if ( setToZero )
        for ( unsigned a=0; a<rows; ++a )
            for ( unsigned k=0; k<numFunctions; ++k )
                mat( a, k ) = corlib::detail_::AccessorTraits< ENTRY >::zero( );
    return;
}



//------------------------------------------------------------------------------
// Write (rows x NF)-matrix of ENTRYs to stream
template< typename BELEMENT, typename MAT >
template< typename ENTRY >
std::ostream & gshell::fem::ElementKLStatic<BELEMENT,MAT>::writeMatrixNF_(
    std::ostream & os,
    const Eigen::Matrix<ENTRY, Eigen::Dynamic, Eigen::Dynamic> & mat,
    const std::string desc
    ) const
{
    os << desc;
    for ( unsigned d=0; d<mat.rows(); ++d ) {
        for ( unsigned l=0; l<mat.cols(); ++l ) {
            os << "<" << d << "," << l << "> " << mat( d, l ) << std::endl;
        }
    }
    return os;
}

//------------------------------------------------------------------------------
// Compute vectors and their derivatives w.r.t. nodal DOFs
template< typename BELEMENT, typename MAT >
void gshell::fem::ElementKLStatic<BELEMENT,MAT>::vectors_(
    const enum gshell::fem::config conf,
    const VecLDim & xi,
    const VecNF & phi,
    const MatLDimNF & dPhiDXi,
    const MatSDimNF & ddPhiDDXi,
    double aCur[dim][dim],
    double bCur[dim][localDim][dim],
    double pCur[dim],
    double pLGrad[localDim][dim]
    ) const
{
    const unsigned numFunctions = this->BasisElement::numFunctions( );

    // Collect nodal quantities
    MatDimNF xNF;  this->BasisElement::supportNodeCoordinates( xNF );
    if ( conf == gshell::fem::CURRENT ) {
        this->updateNodalQuantity_( xNF, 1.0, std::mem_fun( &Node::giveDisplacements ) );
    }
    else if ( conf == LASTCONVERGED ) {
        this->updateNodalQuantity_( xNF,  1.0, std::mem_fun( &Node::giveDisplacements ) );
        this->updateNodalQuantity_( xNF, -1.0, std::mem_fun( &Node::giveIncrement ) );
    }
    else {
        // this IMPLIES that the INITIAL DISPLACEMENTS are ZERO
        ;  // do nothing
    }

    // Compute base
    const MatLDimNF dPhiDXiT = dPhiDXi.transpose();
    const MatDimNF  xNFT     = xNF.transpose();
    SurfaceGeometer::covariantBase( numFunctions, &(dPhiDXiT(0,0)), &(xNFT(0,0)), aCur );

    // Local gradient of base vectors
    MatSDimNF ddPhiDDXiT = ddPhiDDXi.transpose();
    SurfaceGeometer::covariantBaseDeriv1( numFunctions, &(ddPhiDDXiT(0,0)), &(xNFT(0,0)),
                                          aCur, bCur );

    // \vec{p} = \vec{a}_1 \times \vec{a}_2
    la::crossProduct( pCur, aCur[ 0 ], aCur[ 1 ] );
    double pLen; la::norm2<dim>( pLen, pCur );
    SurfaceGeometer::crossprodVectorDeriv1( aCur[ 0 ], bCur[ 0 ], aCur[ 1 ], bCur[ 1 ],
                                            pLGrad );
    
    // done
    return;
}

//------------------------------------------------------------------------------
// Compute vectors and their derivatives w.r.t. nodal DOFs
template< typename BELEMENT, typename MAT >
void gshell::fem::ElementKLStatic<BELEMENT,MAT>::vectorsGradCur_(
    const unsigned & d,
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
    double pLGradDU_dl[localDim][dim]
    ) const
{
    // co-variant base vectors
    const double (&a0Cur)[dim] = aCur[ 0 ];
    const double (&a1Cur)[dim] = aCur[ 1 ];
    const double (&a2Cur)[dim] = aCur[ 2 ];

    // local gradient of co-variant base vectors \vec{a}_{\alpha,\beta}
    const double (&b0Cur)[localDim][dim] = bCur[ 0 ];
    const double (&b1Cur)[localDim][dim] = bCur[ 1 ];

    // 1st derivative of \vec{a}_\alpha w.r.t. nodal DOFs
    double (&a0DU_dl)[dim] = aDU_dl[ 0 ];
    la::zero<dim>( a0DU_dl );
    a0DU_dl[ d ] = dPhiDXi( 0, l );
    double (&a1DU_dl)[dim] = aDU_dl[ 1 ];
    la::zero<dim>( a1DU_dl );
    a1DU_dl[ d ] = dPhiDXi( 1, l );

    // 1st derivative of local gradient of co-variant base vectors \vec{a}_{\alpha,\beta}
    // w.r.t. nodal DOFs
    double (&b0DU_dl)[localDim][dim] = bDU_dl[ 0 ];
    double (&b1DU_dl)[localDim][dim] = bDU_dl[ 1 ];
    la::zero<localDim,dim>( b0DU_dl );
    la::zero<localDim,dim>( b1DU_dl );
    for ( unsigned be=0; be<localDim; ++be ) {
        b0DU_dl[ be ][ d ] = ddPhiDDXi( voigtForward[ 0 ][ be ], l );
        b1DU_dl[ be ][ d ] = ddPhiDDXi( voigtForward[ 1 ][ be ], l );
    }

    // \vec{p} = \vec{a}_1 \times \vec{a}_2
    // P^A{}_D{}^L = \frac{ \partial p^A }{ \partial u^D{}_L }
    SurfaceGeometer::crossprodVectorGradDof( a0Cur, a0DU_dl, a1Cur, a1DU_dl,
                                             pDU_dl );
    // \frac{ \partial p^A{}_{,\beta} }{ \partial u^D{}_L }
    SurfaceGeometer::crossprodVectorDeriv1GradDof( a0Cur, a0DU_dl, b0Cur, b0DU_dl,
                                                   a1Cur, a1DU_dl, b1Cur, b1DU_dl,
                                                   pLGradDU_dl );
    
    // s^A{}_{\beta} = p^A{}_{,\beta} / |p^A|
    double pLen;
    la::norm2<dim>( pLen, pCur );
    double sCur[localDim][dim];
    la::assign<localDim,dim>( sCur, pLGrad );
    la::divide<localDim,dim>( sCur, pLen );
        
    // DU^D{}_D{}^L = \frac{ \partial d^A }{ \partial u^D{}_L }
    SurfaceGeometer::unitVectorGradDof( pLen, pDU_dl, a2Cur, aDU_dl[ 2 ] );
    
    // CU^A{}_\beta{}_D{}^L = \frac{ \partial c^A{}_\beta }{ \partial u^D{}_L }
    SurfaceGeometer::unitVectorDeriv1GradDof( pLen, pDU_dl, pLGradDU_dl,
                                              a2Cur, aDU_dl[ 2 ], sCur, bDU_dl[ 2 ] );

    // done
    return;
}

//------------------------------------------------------------------------------
// Computes the Green-Lagrange strain resultants
// \f$E_{ij} = \alpha_{ij} + \xi^3\beta_{ij}\f$
// (components in reference curvi-linear system \f$\{\bar{\vec{g}}^i\}\f$)
template< typename BELEMENT, typename MAT >
void gshell::fem::ElementKLStatic<BELEMENT,MAT>::strainRes_(
    const double aRef[dim][dim],
    const double cRef[localDim][dim],
    const double aCur[dim][dim],
    const double cCur[localDim][dim],
    double alpha[dim][dim],
    double beta[dim][dim]
    ) const
{
    // membrane strain resultants al,be=be,al
    la::innerProduct<dim>(         alpha[ 0 ][ 0 ], aCur[ 0 ], aCur[ 0 ] );
    la::subtractInnerProduct<dim>( alpha[ 0 ][ 0 ], aRef[ 0 ], aRef[ 0 ] );
    la::innerProduct<dim>(         alpha[ 0 ][ 1 ], aCur[ 0 ], aCur[ 1 ] );
    la::subtractInnerProduct<dim>( alpha[ 0 ][ 1 ], aRef[ 0 ], aRef[ 1 ] );
    la::innerProduct<dim>(         alpha[ 1 ][ 0 ], aCur[ 1 ], aCur[ 0 ] );
    la::subtractInnerProduct<dim>( alpha[ 1 ][ 0 ], aRef[ 1 ], aRef[ 0 ] );
    la::innerProduct<dim>(         alpha[ 1 ][ 1 ], aCur[ 1 ], aCur[ 1 ] );
    la::subtractInnerProduct<dim>( alpha[ 1 ][ 1 ], aRef[ 1 ], aRef[ 1 ] );
    // out-of-plane shear al,3=3,al
    la::innerProduct<dim>(         alpha[ 0 ][ 2 ], aCur[ 0 ], aCur[ 2 ] );
    la::subtractInnerProduct<dim>( alpha[ 0 ][ 2 ], aRef[ 0 ], aRef[ 2 ] );
    la::innerProduct<dim>(         alpha[ 1 ][ 2 ], aCur[ 1 ], aCur[ 2 ] );
    la::subtractInnerProduct<dim>( alpha[ 1 ][ 2 ], aRef[ 1 ], aRef[ 2 ] );
    la::innerProduct<dim>(         alpha[ 2 ][ 0 ], aCur[ 0 ], aCur[ 2 ] );
    la::subtractInnerProduct<dim>( alpha[ 2 ][ 0 ], aRef[ 0 ], aRef[ 2 ] );
    la::innerProduct<dim>(         alpha[ 2 ][ 1 ], aCur[ 1 ], aCur[ 2 ] );
    la::subtractInnerProduct<dim>( alpha[ 2 ][ 1 ], aRef[ 1 ], aRef[ 2 ] );
    // remaining entry is zero 3,3
    alpha[ 2 ][ 2 ] = 0.0;
    // scale
    la::multiply<dim,dim>( alpha, 0.5 );

    // bending strain resultants al,be=be,al
    la::zero<dim,dim>( beta );
    for ( unsigned al = 0; al < localDim; ++al ) {
        for ( unsigned be = 0; be < localDim; ++be ) {
            double tmp_albe = 0.0;
            la::addInnerProduct<dim>(      tmp_albe, aCur[ al ], cCur[ be ] );
            la::addInnerProduct<dim>(      tmp_albe, cCur[ al ], aCur[ be ] );
            la::subtractInnerProduct<dim>( tmp_albe, aRef[ al ], cRef[ be ] );
            la::subtractInnerProduct<dim>( tmp_albe, cRef[ al ], aRef[ be ] );
            beta[ al ][ be ] = 0.5 * tmp_albe;
        }
    }

    // remaining entries are zero
    
    return;
}

//------------------------------------------------------------------------------
// Computes the Green-Lagrange strain resultants
// and their derivatives w.r.t. nodal DOFs
template< typename BELEMENT, typename MAT >
void gshell::fem::ElementKLStatic<BELEMENT,MAT>::strainResGradCur_(
    const unsigned & d,
    const unsigned & l,
    const double aCur[dim][dim],
    const double aDU_dl[dim][dim],
    const double cCur[localDim][dim],
    const double cDU_dl[localDim][dim],
    double etaU_dl[dim][dim],
    double zetaU_dl[dim][dim]
    ) const
{

    // al,be
    for ( unsigned al=0; al<localDim; ++al ) {
        for ( unsigned be=0; be<localDim; ++be ) {
            etaU_dl[ al ][ be ] = 0.5 * (
                aCur[ be ][ d ] * aDU_dl[ al ][ d ] +
                aCur[ al ][ d ] * aDU_dl[ be ][ d ]
                );
            zetaU_dl[ al ][ be ] =
                cCur[ be ][ d ] * aDU_dl[ al ][ d ] +
                cCur[ al ][ d ] * aDU_dl[ be ][ d ];
            la::addInnerProduct<dim>( zetaU_dl[ al ][ be ], aCur[ al ], cDU_dl[ be ] );
            la::addInnerProduct<dim>( zetaU_dl[ al ][ be ], aCur[ be ], cDU_dl[ al ] );
            zetaU_dl[ al ][ be ] *= 0.5;
        }
    }

    // al,3=3,al
    for ( unsigned al=0; al<localDim; ++al ) {
        double etaU_dlal2 = aCur[ 2 ][ d ] * aDU_dl[ al ][ d ];
        la::addInnerProduct<dim>( etaU_dlal2, aCur[ al ], aDU_dl[ 2 ] );
        etaU_dlal2 *= 0.5;
        etaU_dl[ al ][ 2 ] = etaU_dlal2;
        etaU_dl[ 2 ][ al ] = etaU_dlal2;
        zetaU_dl[ al ][ 2 ] = 0.0;
        zetaU_dl[ 2 ][ al ] = 0.0;
    }

    // 3,3
    etaU_dl[ 2 ][ 2 ] = 0.;
    zetaU_dl[ 2 ][ 2 ] = 0.;

    // done
    return;
}


//------------------------------------------------------------------------------
// Computes the elasticity 4-tensor \f$H^{ijkl}\f$
// (components in co-variant reference base \f$\{\bar{\vec{g}}_i\}\f$)
template< typename BELEMENT, typename MAT >
void gshell::fem::ElementKLStatic<BELEMENT,MAT>::elasticityTensorVK_(
    const double aRef[dim][dim],
    double hRef[dim][dim][dim][dim]
    ) const
{
    // material constants 
    const double E = this->giveYoungsModulus( );
    const double nu = this->givePoissonRatio( );
    const double lam = nu*E/(1.+nu)/(1.-2.*nu);
    const double mu = 0.5*E/(1.+nu);
    //const double fact = 2.*mu*lam/(2.*mu+lam);  // due to "elmination"
    const double fact = lam;  // due to "application"

    // shear correction factor
    const double kappa = 5./6.;
//    const double kappa = 1.;

    // metric coefficients [\bar{a}_{al be}]
    double aInv3x3[dim][dim];
    la::productNT<dim,dim>( aInv3x3, aRef, aRef );
    double aInv[localDim][localDim];
    la::subAssign<localDim,localDim,dim,dim>( aInv, aInv3x3 );
    // in-place inversion [\bar{a}^{al be}]
    la::inverse<localDim>( aInv );

    // blank elasticity tensor
    la::zero<dim,dim,dim,dim>( hRef );

    // H^{\alpha\beta\gamma\delta}
    for ( unsigned al=0; al<localDim; ++al ) {
        for ( unsigned be=0; be<localDim; ++be ) {
            for ( unsigned ga=0; ga<localDim; ++ga ) {
                for ( unsigned de=0; de<localDim; ++de ) {
                      hRef[ al ][ be ][ ga ][ de ] = 
                          mu * aInv[ al ][ ga ] * aInv[ be ][ de ] +
                          mu * aInv[ al ][ de ] * aInv[ be ][ ga ] +
                          fact * aInv[ al ][ be ] * aInv[ ga ][ de ];
                }
            }
        }
    }

    // H^{\alpha 3 \gamma 3}
    for ( unsigned al=0; al<localDim; ++al ) {
        for ( unsigned ga=0; ga<localDim; ++ga ) {
            const double hRef_alga = kappa * mu * aInv[ al ][ ga ];
            hRef[ al ][  2 ][ ga ][  2 ] = hRef_alga;
            hRef[ al ][  2 ][  2 ][ ga ] = hRef_alga;
            hRef[  2 ][ al ][ ga ][  2 ] = hRef_alga;
            hRef[  2 ][ al ][  2 ][ ga ] = hRef_alga;
        }
    }

    // done
    return;
}

//------------------------------------------------------------------------------
// Computes the 2nd Piola-Kirchhoff stress resultants
// \f$n^{ij}\f$ and \f$m^{ij}\f$
template< typename BELEMENT, typename MAT >
void gshell::fem::ElementKLStatic<BELEMENT,MAT>::stressResVK_(
    const double aRef[dim][dim],
    const double alphaCur[dim][dim],
    const double betaCur[dim][dim],
    double nCur[dim][dim],
    double mCur[dim][dim]
    ) const
{

    // elasticity tensor
    double hRef[dim][dim][dim][dim];  this->elasticityTensorVK_( aRef, hRef );

    // thickness
    const double & t = thickness_;
    // second moment of area
    const double s = t*t*t / 12.;

    // membrane and bending stress resultant
    for ( unsigned i=0; i<dim; ++i ) {
        for ( unsigned j=0; j<dim; ++j ) {
            la::innerProduct<dim,dim>( nCur[ i ][ j ], hRef[ i ][ j ], alphaCur );
            nCur[ i ][ j ] *= t;
            la::innerProduct<dim,dim>( mCur[ i ][ j ], hRef[ i ][ j ], betaCur );
            mCur[ i ][ j ] *= s;
        }
    }

    return;
}

//------------------------------------------------------------------------------
// Compute nodal forces due to body force
template< typename BELEMENT, typename MAT >
template< typename FUNC >
void gshell::fem::ElementKLStatic<BELEMENT,MAT>::bodyForceDetailed(
    const VecLDim & xi,
    const double & weight,
    FUNC func,      // physical coord => force density N/m^3 or N/m^2 (see below)
    const double & factor,
    const enum gshell::fem::config conf,
    const bool isThicknessIntegrated
    )
{
    // constants
    const unsigned numFunctions = this->BasisElement::numFunctions( );

    // nodal coordinates
    MatDimNF xNF;  this -> supportNodeCoordinates( xNF );
    if ( conf == gshell::fem::CURRENT ) {
        MatDimNF uNF;  this -> nodalDisplacements( uNF );
        xNF += uNF;
    }

    // shape functions
    VecNF phi; this -> sfun( xi, phi );
    // location
    VecDim x;
    const MatDimNF xNFT = xNF.transpose();
    SurfaceGeometer::interpolate( numFunctions,
                                  &(phi(0)), &(xNFT(0,0)),
                                  &(x(0)) );

    // shape function derivatives
    MatLDimNF dPhiDXi; this -> sfunGrad( xi, dPhiDXi );
    // tangents
    double tang[localDim][dim];
    const MatLDimNF dPhiDXiT = dPhiDXi.transpose();
    SurfaceGeometer::interpolate( numFunctions,
                                  &(dPhiDXiT(0,0)), &(xNFT(0,0)),
                                  tang[0] );
    SurfaceGeometer::interpolate( numFunctions,
                                  &(dPhiDXiT(0,1)), &(xNFT(0,0)),
                                  tang[1] );

    // normal
    double normal[dim];
    la::crossProduct( normal, tang[0], tang[1] );
    double metric;
    la::norm2<dim>( metric, normal );
    if ( not isThicknessIntegrated ) metric *= thickness_;

    // load density
    const VecDim forceDensityConst = factor * func( x );

    // pass each nodal force to the nodes
    for ( unsigned l=0; l<numFunctions; ++l ) {
        VecDof nodalForce;
        nodalForce.setZero( );
        for ( unsigned d=0; d<dim; ++d ) {
            // 1st derivatives of vectors
            double uDU_dl[dim];
            la::zero<dim>( uDU_dl );
            uDU_dl[ d ] = phi( l );

            // nodal force
            double tmp;
            la::innerProduct<dim>( tmp, &(forceDensityConst(0)), uDU_dl );
            nodalForce[ d ] = tmp * metric * weight;
        }
        
        // store nodal force
        this->BasisElement::giveSupportNodePtr( l )->addToForce( nodalForce );  // negative RHS

    }

    return;
}

//------------------------------------------------------------------------------
// Compute nodal forces due to body force
template< typename BELEMENT, typename MAT >
template< typename FUNC >
void gshell::fem::ElementKLStatic<BELEMENT,MAT>::bodyForceLin(
    const VecLDim & xi,
    const double & weight,
    FUNC func, 
    const double & factor
    )
{
    // constants
    const unsigned numFunctions = this->BasisElement::numFunctions( );

    // material point
    const VecDim xRef = this->giveCoordinate( xi );

    // force density
    const VecDim forceDensityLin   = factor * func( xRef );  // N/m^3

    // collect shape functions and their derivatives
    VecNF phi;              //this->BasisElement::sfun( xi, phi );
    MatLDimNF dPhiDXi;      //this->BasisElement::sfunGrad( xi, dPhiDXi );
    MatSDimNF ddPhiDDXi;    //this->BasisElement::sfunHess( xi, ddPhiDDXi );
    this->BasisElement::sfunGradHess( xi, phi, dPhiDXi, ddPhiDDXi );

    // get reference mid-vectors
    double aRef[dim][dim];       double bRef[dim][localDim][dim];
    double pRef[dim];            double pRefLGrad[localDim][dim];
    this->vectors_( gshell::fem::REFERENCE, xi, phi, dPhiDXi, ddPhiDDXi,
                    aRef, bRef, pRef, pRefLGrad );

    // get current mid-vectors
    double aCur[dim][dim];       double bCur[dim][localDim][dim];
    double pCur[dim];            double pLGrad[localDim][dim];
    this->vectors_( gshell::fem::CURRENT, xi, phi, dPhiDXi, ddPhiDDXi,
                    aCur, bCur, pCur, pLGrad );

    // convenience
    const double (&a0Cur)[dim] = aCur[ 0 ];
    const double (&a1Cur)[dim] = aCur[ 1 ];
    const double (&dCur)[dim]  = aCur[ 2 ];
    double pLen;  la::norm2<dim>( pLen, pCur );
    double zeroDim[dim];  la::zero<dim>( zeroDim );

    // surface metric
    double detA = la::determinant<dim>( aRef );

    // thickness
    const double & t = thickness_;
    // second moment of area
    const double s = t*t*t / 12.;

    // pass each nodal force to the nodes
    for ( unsigned l=0; l<numFunctions; ++l ) {
        VecDof nodalForce;
        nodalForce.clear( );
        for ( unsigned d=0; d<dim; ++d ) {
            // 1st derivatives of vectors
            double uDU_dl[dim];
            la::zero<dim>( uDU_dl );
            uDU_dl[ d ] = phi( l );

            // nodal force
            nodalForce[ d ] = 0.;

            // external stiffness matrix is need
            if ( forceDensityLin.norm( ) > 1.e-6 ) {

                // 1st derivatives of vectors
                double aDU_dl[dim][dim];      double bDU_dl[dim][localDim][dim];
                double pDU_dl[dim];           double pLGradDU_dl[localDim][dim];
                this->vectorsGradCur_( d, l, xi, phi, dPhiDXi, ddPhiDDXi,
                                       aRef, aCur, bCur, pCur, pLGrad,
                                       aDU_dl, bDU_dl, pDU_dl, pLGradDU_dl );

                // convenience
                const double (&a0DU_dl)[dim] = aDU_dl[ 0 ];
                const double (&a1DU_dl)[dim] = aDU_dl[ 1 ];
                const double (&dDU_dl)[dim]  = aDU_dl[ 2 ];

                // nodal force
                double tmp;
                la::innerProduct<dim>( tmp, &(forceDensityLin(0)), dDU_dl );
                nodalForce( d ) += s * tmp * detA * weight;

            }
        }
        
        // store nodal force
        this->BasisElement::giveSupportNodePtr( l )->addToForce( nodalForce );  // negative RHS

    }

    return;
}


//------------------------------------------------------------------------------
// Compute nodal forces due to body force
template< typename BELEMENT, typename MAT >
template< typename FUNC >
void gshell::fem::ElementKLStatic<BELEMENT,MAT>::bodyForceLinStiffness(
    const VecLDim & xi,
    const double & weight,
    FUNC func, 
    const double & factor,
    Eigen::MatrixXd & elemStiff
    )
{
    // constants
    const unsigned numFunctions = this->BasisElement::numFunctions( );
    const unsigned matSize = dof * numFunctions;
    FTL_VERIFY( elemStiff.rows() == matSize and
                elemStiff.cols() == matSize );

    // material point
    const VecDim xRef = this->giveCoordinate( xi );

    // force density
    const VecDim forceDensityLin   = factor * func( xRef );  // N/m^3 

    // collect shape functions and their derivatives
    VecNF phi;              //this->BasisElement::sfun( xi, phi );
    MatLDimNF dPhiDXi;      //this->BasisElement::sfunGrad( xi, dPhiDXi );
    MatSDimNF ddPhiDDXi;    //this->BasisElement::sfunHess( xi, ddPhiDDXi );
    this->BasisElement::sfunGradHess( xi, phi, dPhiDXi, ddPhiDDXi );

    // get reference mid-vectors
    double aRef[dim][dim];       double bRef[dim][localDim][dim];
    double pRef[dim];            double pRefLGrad[localDim][dim];
    this->vectors_( gshell::fem::REFERENCE, xi, phi, dPhiDXi, ddPhiDDXi,
                    aRef, bRef, pRef, pRefLGrad );

    // get current mid-vectors
    double aCur[dim][dim];       double bCur[dim][localDim][dim];
    double pCur[dim];            double pLGrad[localDim][dim];
    this->vectors_( gshell::fem::CURRENT, xi, phi, dPhiDXi, ddPhiDDXi,
                    aCur, bCur, pCur, pLGrad );

    // convenience
    const double (&a0Cur)[dim] = aCur[ 0 ];
    const double (&a1Cur)[dim] = aCur[ 1 ];
    const double (&dCur)[dim]  = aCur[ 2 ];
    double pLen;  la::norm2<dim>( pLen, pCur );
    double zeroDim[dim];  la::zero<dim>( zeroDim );

    // surface metric
    double detA = la::determinant<dim>( aRef );

    // thickness
    const double & t = thickness_;
    // second moment of area
    const double s = t*t*t / 12.;

    // pass each nodal force to the nodes
    for ( unsigned l=0; l<numFunctions; ++l ) {
        VecDof nodalForce;
        nodalForce.clear( );
        for ( unsigned d=0; d<dim; ++d ) {
            // 1st derivatives of vectors
            double uDU_dl[dim];
            la::zero<dim>( uDU_dl );
            uDU_dl[ d ] = phi( l );

            // external stiffness matrix
            {

                // 1st derivatives of vectors
                double aDU_dl[dim][dim];      double bDU_dl[dim][localDim][dim];
                double pDU_dl[dim];           double pLGradDU_dl[localDim][dim];
                this->vectorsGradCur_( d, l, xi, phi, dPhiDXi, ddPhiDDXi,
                                       aRef, aCur, bCur, pCur, pLGrad,
                                       aDU_dl, bDU_dl, pDU_dl, pLGradDU_dl );

                // convenience
                const double (&a0DU_dl)[dim] = aDU_dl[ 0 ];
                const double (&a1DU_dl)[dim] = aDU_dl[ 1 ];
                const double (&dDU_dl)[dim]  = aDU_dl[ 2 ];

                for ( unsigned k=0; k<numFunctions; ++k ) {
                    for ( unsigned c=0; c<dim; ++c ) {

                        // 1st derivatives of vectors
                        double aDU_ck[dim][dim];      double bDU_ck[dim][localDim][dim];
                        double pDU_ck[dim];           double pLGradDU_ck[localDim][dim];
                        this->vectorsGradCur_( c, k, xi, phi, dPhiDXi, ddPhiDDXi,
                                               aRef, aCur, bCur, pCur, pLGrad,
                                               aDU_ck, bDU_ck, pDU_ck, pLGradDU_ck );

                        // convenience
                        const double (&a0DU_ck)[dim] = aDU_ck[ 0 ];
                        const double (&a1DU_ck)[dim] = aDU_ck[ 1 ];
                        const double (&dDU_ck)[dim]  = aDU_ck[ 2 ];

                        double pDUDU_dlck[dim];
                        SurfaceGeometer::crossprodVectorHessDof( a0Cur, a0DU_dl, a0DU_ck, zeroDim,
                                                                 a1Cur, a1DU_dl, a1DU_ck, zeroDim,
                                                                 pDUDU_dlck );

                        double dDUDU_dlck[dim];
                        SurfaceGeometer::unitVectorHessDof( pLen, pDU_dl, pDU_ck, pDUDU_dlck,
                                                            dCur, dDU_dl, dDU_ck, dDUDU_dlck );

                        double kUU_dlck;
                        la::innerProduct<dim>( kUU_dlck, &(forceDensityLin(0)), dDUDU_dlck );
                        kUU_dlck *= s;

                        elemStiff( l*dof+d, k*dof+c ) -= kUU_dlck * detA * weight;  // positive LHS

                    }
                }
            }
        }
        
    }

    return;
}

//------------------------------------------------------------------------------
// Compute nodal forces due to body force
template< typename BELEMENT, typename MAT >
template< typename FUNC >
void gshell::fem::ElementKLStatic<BELEMENT,MAT>::externalPressure(
    const VecLDim & xi,
    const double & weight,
    FUNC func,
    const double & factor
    )
{
    // constants
    const unsigned numFunctions = this->BasisElement::numFunctions( );

    // material point
    const VecDim xRef = this->giveCoordinate( xi );

    // force density
    const double pressure = factor * func( xRef );  // N/m^2

    // collect shape functions and their derivatives
    VecNF phi;              //this->BasisElement::sfun( xi, phi );
    MatLDimNF dPhiDXi;      //this->BasisElement::sfunGrad( xi, dPhiDXi );
    MatSDimNF ddPhiDDXi;    //this->BasisElement::sfunHess( xi, ddPhiDDXi );
    this->BasisElement::sfunGradHess( xi, phi, dPhiDXi, ddPhiDDXi );

    // Compute reference base
    MatDimNF xNF;  this->BasisElement::supportNodeCoordinates( xNF );
    double aRef[dim][dim];
    const MatLDimNF dPhiDXiT = dPhiDXi.transpose();
    const MatDimNF xNFT = xNF.transpose();
    const double detA = SurfaceGeometer::covariantBase( numFunctions, &(dPhiDXiT(0,0)),
                                                        &(xNFT(0,0)), aRef );

    // Compute current base
    this->updateNodalQuantity_( xNF, 1.0, std::mem_fun( &Node::giveDisplacements ) );
    double aCur[dim][dim];
    const MatDimNF xNFTcur = xNF.transpose();
    SurfaceGeometer::covariantBase( numFunctions, &(dPhiDXiT(0,0)), &(xNFTcur(0,0)), aCur );

    // convenience
    const double (&a2Cur)[dim] = aCur[ 2 ];

    // current surface pressure vector
    double normalPressure[dim];  // N/m^2
    la::assign<dim>( normalPressure, a2Cur );
    la::multiply<dim>( normalPressure, pressure );
    

    // pass each nodal force to the nodes
    VecDof nodalForce;
    // double uDU_dl[dim];
    double tmp;
    for ( unsigned l=0; l<numFunctions; ++l ) {
        
        for ( unsigned d=0; d<dim; ++d ) {
            // 1st derivatives of vectors
            // const VecDim uDU_dl = phi( l ) * ublas::column( delta, d );
            // la::zero<dim>( uDU_dl );  uDU_dl[ d ] = phi( l );

            // nodal force
            //la::innerProduct<dim>( tmp, normalPressure, uDU_dl );
            tmp = normalPressure[ d ] * phi( l );
            nodalForce[ d ] = tmp * detA * weight;
        }

        // store nodal force
        this->BasisElement::giveSupportNodePtr( l )->addToForce( nodalForce );  // negative RHS
    }

    return;
}

//------------------------------------------------------------------------------
// Compute nodal forces due to body force
template< typename BELEMENT, typename MAT >
template< typename FUNC >
void gshell::fem::ElementKLStatic<BELEMENT,MAT>::externalPressureStiffness(
    const VecLDim & xi,
    const double & weight,
    FUNC func,
    const double & factor,
    Eigen::MatrixXd & elemStiff
    )
{
    // constants
    const unsigned numFunctions = this->BasisElement::numFunctions( );
    const unsigned matSize = dof * numFunctions;
    FTL_VERIFY( elemStiff.rows() == matSize and
                elemStiff.cols() == matSize );

    // material point
    const VecDim xRef = this->giveCoordinate( xi );

    // force density
    const double pressure = factor * func( xRef );  // N/m^2

    // collect shape functions and their derivatives
    VecNF phi;              //this->BasisElement::sfun( xi, phi );
    MatLDimNF dPhiDXi;      //this->BasisElement::sfunGrad( xi, dPhiDXi );
    MatSDimNF ddPhiDDXi;    //this->BasisElement::sfunHess( xi, ddPhiDDXi );
    this->BasisElement::sfunGradHess( xi, phi, dPhiDXi, ddPhiDDXi );

    // Compute reference base
    MatDimNF xNF;  this->BasisElement::supportNodeCoordinates( xNF );
    double aRef[dim][dim];
    const MatLDimNF dPhiDXiT = dPhiDXi.transpose();
    const MatDimNF xNFT = xNF.transpose();
    const double detA = SurfaceGeometer::covariantBase( numFunctions, &(dPhiDXiT(0,0)),
                                                        &(xNFT(0,0)), aRef );
    // Compute current base
    this->updateNodalQuantity_( xNF, 1.0, std::mem_fun( &Node::giveDisplacements ) );
    double aCur[dim][dim];
    
    xNFT = xNF.transpose();
    SurfaceGeometer::covariantBase( numFunctions, &(dPhiDXiT(0,0)), &(xNFT(0,0)), aCur );

    // convenience
    const double (&a2Cur)[dim] = aCur[ 2 ];

    // current surface pressure vector
    double normalPressure[dim];  // N/m^2
    la::assign<dim>( normalPressure, a2Cur );
    la::multiply<dim>( normalPressure, pressure );
    

    // // pass each nodal force to the nodes
    // VecDof nodalForce;
    // // double uDU_dl[dim];
    // double tmp;
    // for ( unsigned l=0; l<numFunctions; ++l ) {
        
    //     for ( unsigned d=0; d<dim; ++d ) {
    //         // 1st derivatives of vectors
    //         // const VecDim uDU_dl = phi( l ) * ublas::column( delta, d );
    //         // la::zero<dim>( uDU_dl );  uDU_dl[ d ] = phi( l );

    //         // nodal force
    //         //la::innerProduct<dim>( tmp, normalPressure, uDU_dl );
    //         tmp = normalPressure[ d ] * phi( l );
    //         nodalForce[ d ] = tmp * detA * weight;
    //     }

    //     // store nodal force
    //     this->BasisElement::giveSupportNodePtr( l )->addToForce( nodalForce );  // negative RHS
    // }

    // external stiffness matrix (non-symmetric)

    // convenience
    const double (&a0Cur)[dim] = aCur[ 0 ];
    const double (&a1Cur)[dim] = aCur[ 1 ];
    double pCur[dim];  la::crossProduct( pCur, a0Cur, a1Cur );
    double pLen;       la::norm2<dim>( pLen, pCur );

    for ( unsigned l=0; l<numFunctions; ++l ) {
        for ( unsigned d=0; d<dim; ++d ) {

            // 1st derivatives of vectors
            double uDU_dl[dim];
            la::zero<dim>( uDU_dl );  uDU_dl[ d ] = phi( l );
        
            FTL_VERIFY_DESCRIPTIVE( false, "WARNING : not verified\n" );

            for ( unsigned k=0; k<numFunctions; ++k ) {
                for ( unsigned c=0; c<dim; ++c ) {

                    double a0DU_ck[dim], a1DU_ck[dim];
                    la::zero<dim>( a0DU_ck );  a0DU_ck[ c ] = dPhiDXi( 0, k );
                    la::zero<dim>( a1DU_ck );  a1DU_ck[ c ] = dPhiDXi( 1, k );

                    // \vec{p} = \vec{a}_1 \times \vec{a}_2
                    // P^A{}_D{}^L = \frac{ \partial p^A }{ \partial u^D{}_L }
                    double pDU_ck[dim];
                    SurfaceGeometer::crossprodVectorGradDof( a0Cur, a0DU_ck, a1Cur, a1DU_ck,
                                                             pDU_ck );
                    // \vec{a}_3 = \vec{p}/|\vec{p}|
                    // (A_3)^A{}_D{}^L = \frac{ \partial a_3^A }{ \partial u^D{}_L }
                    double a2DU_ck[dim];
                    SurfaceGeometer::unitVectorGradDof( pLen, pDU_ck, a2Cur, a2DU_ck );

                    // add to stiffness
                    double kUU_dlck;
                    la::innerProduct<dim>( kUU_dlck, uDU_dl, a2DU_ck );
                    kUU_dlck *= pressure;
                    elemStiff( l*dof+d, k*dof+c ) -= kUU_dlck * detA * weight;  // positive LHS
                    
                }
            }

        }
    }

    return;
}

//------------------------------------------------------------------------------
// Evaluate internal force kernel at xi (store results in nodes)
template< typename BELEMENT, typename MAT >
void gshell::fem::ElementKLStatic<BELEMENT,MAT>::internalForceIntegrand_(
    const VecLDim & xi,
    const double & weight,
    const enum gshell::fem::config & conf,
    const double & factor
    )
{
    // Collect shape functions and their derivatives
    VecNF phi;              //this->BasisElement::sfun( xi, phi );
    MatLDimNF dPhiDXi;      //this->BasisElement::sfunGrad( xi, dPhiDXi );
    MatSDimNF ddPhiDDXi;    //this->BasisElement::sfunHess( xi, ddPhiDDXi );
    this->BasisElement::sfunGradHess( xi, phi, dPhiDXi, ddPhiDDXi );

    // get reference mid-vectors
    double aRef[dim][dim];       double bRef[dim][localDim][dim];
    double pRef[dim];            double pRefLGrad[localDim][dim];
    this->vectors_( gshell::fem::REFERENCE, xi, phi, dPhiDXi, ddPhiDDXi,
                    aRef, bRef, pRef, pRefLGrad );

    // get current mid-vectors
    double aCur[dim][dim];       double bCur[dim][localDim][dim];
    double pCur[dim];            double pLGrad[localDim][dim];
    this->vectors_( conf, xi, phi, dPhiDXi, ddPhiDDXi,
                    aCur, bCur, pCur, pLGrad );

    // Green-Lagrange strain resultants
    double alphaCur[dim][dim], betaCur[dim][dim];
    this->strainRes_( aRef, bRef[ 2 ], aCur, bCur[ 2 ],
                      alphaCur, betaCur );

    // 2nd PK stress resultants
    double nCur[dim][dim], mCur[dim][dim];
    this->stressResVK_( aRef, alphaCur, betaCur, nCur, mCur );

    // metric
    const double detA = la::determinant<dim>( aRef );

    // pass each nodal force to the nodes
    const unsigned numFunctions = this->BasisElement::numFunctions();
    for ( unsigned l=0; l<numFunctions; ++l ) {
        VecDof nodalForce;
        for ( unsigned d=0; d<dim; ++d ) {

            // 1st derivatives of vectors
            double aDU_dl[dim][dim];      double bDU_dl[dim][localDim][dim];
            double pDU_dl[dim];           double pLGradDU_dl[localDim][dim];
            this->vectorsGradCur_( d, l, xi, phi, dPhiDXi, ddPhiDDXi,
                                   aRef, aCur, bCur, pCur, pLGrad,
                                   aDU_dl, bDU_dl, pDU_dl, pLGradDU_dl );

            // 1st derivatives of strains w.r.t. nodal DOFs
            double alphaDU_dl[dim][dim], betaDU_dl[dim][dim];
            this->strainResGradCur_( d, l, 
                                     aCur, aDU_dl, bCur[ 2 ], bDU_dl[ 2 ],
                                     alphaDU_dl, betaDU_dl );

            // nodal force
            double tmp;
            la::innerProduct<dim,dim>( tmp, nCur, alphaDU_dl );
            la::addInnerProduct<dim,dim>( tmp, mCur, betaDU_dl );
            nodalForce[ d ] = tmp * detA * weight;

        }
        nodalForce *= -factor;
        this->BasisElement::giveSupportNodePtr( l )->addToForce( nodalForce );  // negative RHS
    }

    return;
}

//------------------------------------------------------------------------------
// Evaluate stiffness kernel at xi, weigh result and store it on elemStiff
template< typename BELEMENT, typename MAT >
void gshell::fem::ElementKLStatic<BELEMENT,MAT>::stiffnessIntegrand(
    const VecLDim & xi,
    const double & weight,
    Eigen::MatrixXd & elemStiff
    ) const
{
    // number of shape functions
    const unsigned numFunctions = this->BasisElement::numFunctions( );
    const unsigned matSize = dof * numFunctions;
    FTL_VERIFY( elemStiff.rows() == matSize and
                elemStiff.cols() == matSize );

    // Collect shape functions and their derivatives
    VecNF phi;              //this->BasisElement::sfun( xi, phi );
    MatLDimNF dPhiDXi;      //this->BasisElement::sfunGrad( xi, dPhiDXi );
    MatSDimNF ddPhiDDXi;    //this->BasisElement::sfunHess( xi, ddPhiDDXi );
    this->BasisElement::sfunGradHess( xi, phi, dPhiDXi, ddPhiDDXi );

    // get reference mid-vectors
    double aRef[dim][dim];       double bRef[dim][localDim][dim];
    double pRef[dim];            double pRefLGrad[localDim][dim];
    this->vectors_( gshell::fem::REFERENCE, xi, phi, dPhiDXi, ddPhiDDXi,
                    aRef, bRef, pRef, pRefLGrad );

    // get current mid-vectors
    double aCur[dim][dim];       double bCur[dim][localDim][dim];
    double pCur[dim];            double pLGrad[localDim][dim];
    this->vectors_( gshell::fem::CURRENT, xi, phi, dPhiDXi, ddPhiDDXi,
                    aCur, bCur, pCur, pLGrad );

    // metric
    const double detA = la::determinant<dim>( aRef );

    const double (&a0Cur)[dim] = aCur[ 0 ];
    const double (&a1Cur)[dim] = aCur[ 1 ];
    const double (&dCur)[dim]  = aCur[ 2 ];
    double pLen;  la::norm2<dim>( pLen, pCur );
    double sCur[dim][dim];
    la::assign<dim,dim>( sCur, pLGrad );
    la::divide<dim,dim>( sCur, pLen );
    const double (&b0Cur)[localDim][dim] = bCur[ 0 ];
    const double (&b1Cur)[localDim][dim] = bCur[ 1 ];
    const double (&cCur)[localDim][dim]  = bCur[ 2 ];

    // Green-Lagrange strain resultants
    double alphaCur[dim][dim], betaCur[dim][dim];
    this->strainRes_( aRef, bRef[ 2 ], aCur, cCur,
                      alphaCur, betaCur );

    // 2nd PK stress resultants
    double nCur[dim][dim], mCur[dim][dim];
    this->stressResVK_( aRef, alphaCur, betaCur, nCur, mCur );
    // material elasticity tensor
    double hRef[dim][dim][dim][dim];
    this->elasticityTensorVK_( aRef, hRef );

    // thickness
    const double & t = thickness_;
    // second moment of area
    const double s = t*t*t / 12.;

    // constants
    double zeroDim[dim];                la::zero<dim>( zeroDim );
    double zeroLDimDim[localDim][dim];  la::zero<localDim,dim>( zeroLDimDim );

#define GSHELL_KL_PRECOMPUTEVECTORSFORSTIFFNESS
#ifdef GSHELL_KL_PRECOMPUTEVECTORSFORSTIFFNESS
    //--------------------------------------------------------------------------
    // pre-compute mid-plane vectors and their various derivatives
    double (*aDU)[dim][dim]            = new double[dim*numFunctions][dim][dim];
    double (*bDU)[dim][localDim][dim]  = new double[dim*numFunctions][dim][localDim][dim];
    double (*pDU)[dim]                 = new double[dim*numFunctions][dim];
    double (*pLGradDU)[localDim][dim]  = new double[dim*numFunctions][localDim][dim];
    double (*etaU)[dim][dim]           = new double[dim*numFunctions][dim][dim];
    double (*zetaU)[dim][dim]          = new double[dim*numFunctions][dim][dim];

    for ( unsigned l=0; l<numFunctions; ++l ) {
        for ( unsigned d=0; d<dim; ++d ) {

            // combined index
            const unsigned dl = l*dim + d;

            // 1st derivatives of vectors
            this->vectorsGradCur_( d, l, xi, phi, dPhiDXi, ddPhiDDXi,
                                   aRef, aCur, bCur, pCur, pLGrad,
                                   aDU[ dl ], bDU[ dl ], pDU[ dl ], pLGradDU[ dl ] );

            // 1st derivatives of strains w.r.t. nodal DOFs
            this->strainResGradCur_( d, l,
                                     aCur, aDU[ dl ], cCur, bDU[ dl ][ 2 ],
                                     etaU[ dl ], zetaU[ dl ] );

        }
    }
#endif

    //--------------------------------------------------------------------------
    // stiffness matrix
    for ( unsigned l=0; l<numFunctions; ++l ) {
        for ( unsigned d=0; d<dim; ++d ) {

#ifdef GSHELL_KL_PRECOMPUTEVECTORSFORSTIFFNESS
            const unsigned dl = l*dim + d;

            // base vectors
            const double (&a0DU_dl)[dim] = aDU[ dl ][ 0 ];
            const double (&a1DU_dl)[dim] = aDU[ dl ][ 1 ];
            const double (&dDU_dl)[dim]  = aDU[ dl ][ 2 ];
            const double (&pDU_dl)[dim]  = pDU[ dl ];
            
            // local deriv of base vectors
            const double (&b0DU_dl)[localDim][dim]     = bDU[ dl ][ 0 ];
            const double (&b1DU_dl)[localDim][dim]     = bDU[ dl ][ 1 ];
            const double (&cDU_dl)[localDim][dim]      = bDU[ dl ][ 2 ];
            const double (&pLGradDU_dl)[localDim][dim] = pLGradDU[ dl ];

            // 1st derivatives of strains w.r.t. nodal DOFs
            const double (&etaU_dl)[dim][dim]  = etaU[ dl ];
            const double (&zetaU_dl)[dim][dim] = zetaU[ dl ];
#else
            // 1st derivatives of vectors
            double aDU_dl[dim][dim];     double bDU_dl[dim][localDim][dim];
            double pDU_dl[dim];          double pLGradDU_dl[localDim][dim];
            this->vectorsGradCur_( d, l, xi, phi, dPhiDXi, ddPhiDDXi,
                                   aRef, aCur, bCur, pCur, pLGrad,
                                   aDU_dl, bDU_dl, pDU_dl, pLGradDU_dl );
            
            // base vectors
            const double (&a0DU_dl)[dim] = aDU_dl[ 0 ];
            const double (&a1DU_dl)[dim] = aDU_dl[ 1 ];
            const double (&dDU_dl)[dim]  = aDU_dl[ 2 ];
            // local deriv of base vectors
            const double (&b0DU_dl)[localDim][dim] = bDU_dl[ 0 ];
            const double (&b1DU_dl)[localDim][dim] = bDU_dl[ 1 ];
            const double (&cDU_dl)[localDim][dim]  = bDU_dl[ 2 ];

            // 1st derivatives of strains w.r.t. nodal DOFs
            double etaU_dl[dim][dim], zetaU_dl[dim][dim];
            this->strainResGradCur_( d, l,
                                     aCur, aDU_dl, cCur, cDU_dl,
                                     etaU_dl, zetaU_dl );
#endif

            //------------------------------------------------------------------
            // "elastic and initial-displacement" constribution to stiffness matrix
//#define GSHELL_KL_SYMMETRICSTIFFNESS
#ifdef GSHELL_KL_SYMMETRICSTIFFNESS
            const unsigned kBegin = l;
#else
            const unsigned kBegin = 0;
#endif
            for ( unsigned k=kBegin; k<numFunctions; ++k ) {
                for ( unsigned c=0; c<dim; ++c ) {

#ifdef GSHELL_KL_PRECOMPUTEVECTORSFORSTIFFNESS
                    const unsigned ck = k*dim + c;

                    // base vectors
                    const double (&aDU_ck)[dim][dim] = aDU[ ck ];
                    const double (&pDU_ck)[dim]      = pDU[ ck ];
            
                    // local deriv of base vectors
                    const double (&bDU_ck)[dim][localDim][dim] = bDU[ ck ];
                    const double (&pLGradDU_ck)[localDim][dim] = pLGradDU[ ck ];

                    // 1st derivatives of strains w.r.t. nodal DOFs
                    const double (&etaU_ck)[dim][dim]  = etaU[ ck ];
                    const double (&zetaU_ck)[dim][dim] = zetaU[ ck ];
#else
                    // 1st derivatives of vectors
                    double aDU_ck[dim][dim];     double bDU_ck[dim][localDim][dim];
                    double pDU_ck[dim];          double pLGradDU_ck[localDim][dim];
                    this->vectorsGradCur_( c, k, xi, phi, dPhiDXi, ddPhiDDXi,
                                           aRef, aCur, bCur, pCur, pLGrad,
                                           aDU_ck, bDU_ck, pDU_ck, pLGradDU_ck );

                    // 1st derivatives of strains w.r.t. nodal DOFs
                    double etaU_ck[dim][dim], zetaU_ck[dim][dim];
                    this->strainResGradCur_( c, k,
                                             aCur, aDU_ck, cCur, bDU_ck[ 2 ],
                                             etaU_ck, zetaU_ck );
#endif

                    // "elastic and initial-displacement" stiffness components
                    double kUU_dlck = 0.0;
                    for ( unsigned i=0; i<dim; ++i ) {
                        for ( unsigned j=0; j<dim; ++j ) {
                            double tmp;
                            la::innerProduct<dim,dim>( tmp, hRef[ i ][ j ], etaU_ck );
                            kUU_dlck +=  etaU_dl[ i ][ j ] * t * tmp;
                            la::innerProduct<dim,dim>( tmp, hRef[ i ][ j ], zetaU_ck );
                            kUU_dlck += zetaU_dl[ i ][ j ] * s * tmp;
                        }
                    }
#ifdef GSHELL_KL_SYMMETRICSTIFFNESS
                    if ( k != l )  elemStiff( k*dof+c, l*dof+d ) += kUU_dlck * detA * weight;  // positive LHS
#endif
                    elemStiff( l*dof+d, k*dof+c ) += kUU_dlck * detA * weight;  // positive LHS

                    //----------------------------------------------------------
                    // "geometric" contribution to stiffness matrix
//#define GSHELL_KL_SKIP_GEOMETRICSTIFFNESS
#ifndef GSHELL_KL_SKIP_GEOMETRICSTIFFNESS
                    // base vectors
                    const double (&a0DU_ck)[dim] = aDU_ck[ 0 ];
                    const double (&a1DU_ck)[dim] = aDU_ck[ 1 ];
                    const double (&dDU_ck)[dim]  = aDU_ck[ 2 ];
                    const double (&b0DU_ck)[localDim][dim] = bDU_ck[ 0 ];
                    const double (&b1DU_ck)[localDim][dim] = bDU_ck[ 1 ];
                    const double (&cDU_ck)[localDim][dim]  = bDU_ck[ 2 ];

                    double pDUDU_dlck[dim];
                    SurfaceGeometer::crossprodVectorHessDof( a0Cur, a0DU_dl, a0DU_ck, zeroDim,
                                                             a1Cur, a1DU_dl, a1DU_ck, zeroDim,
                                                             pDUDU_dlck );

                    double pLGradDUDU_dlck[localDim][dim];
                    SurfaceGeometer::crossprodVectorDeriv1HessDof( a0Cur, a0DU_dl, a0DU_ck, zeroDim,
                                                                   b0Cur, b0DU_dl, b0DU_ck, zeroLDimDim,
                                                                   a1Cur, a1DU_dl, a1DU_ck, zeroDim,
                                                                   b1Cur, b1DU_dl, b1DU_ck, zeroLDimDim,
                                                                   pLGradDUDU_dlck );
                    
                    double dDUDU_dlck[dim];
                    SurfaceGeometer::unitVectorHessDof( pLen, pDU_dl, pDU_ck, pDUDU_dlck,
                                                        dCur, dDU_dl, dDU_ck, dDUDU_dlck );

                    double cDUDU_dlck[localDim][dim];
                    SurfaceGeometer::unitVectorDeriv1HessDof( pLen, pDU_dl, pDU_ck, pDUDU_dlck,
                                                              pLGradDU_dl, pLGradDU_ck, pLGradDUDU_dlck,
                                                              dCur, dDU_dl, dDU_ck, dDUDU_dlck,
                                                              sCur, cDUDU_dlck );
                    
                    // "geometric" stiffness components
                    double kgUU_dlck = 0.0;
                    // al,be
                    for ( unsigned al=0; al<localDim; ++al ) {
                        for ( unsigned be=0; be<localDim; ++be ) {
                            if ( d == c )
                                kgUU_dlck += 0.5 * nCur[ al ][ be ] * (
                                    dPhiDXi( be, l ) * dPhiDXi( al, k ) +
                                    dPhiDXi( al, l ) * dPhiDXi( be, k )
                                    );
                            double tmp =
                                dPhiDXi( al, l ) * cDU_ck[ be ][ d ] +
                                dPhiDXi( be, l ) * cDU_ck[ al ][ d ] +
                                cDU_dl[ be ][ c ] * dPhiDXi( al, k ) +
                                cDU_dl[ al ][ c ] * dPhiDXi( be, k );
                            la::addInnerProduct<dim>( tmp, aCur[ al ], cDUDU_dlck[ be ] );
                            la::addInnerProduct<dim>( tmp, aCur[ be ], cDUDU_dlck[ al ] );
                            kgUU_dlck += 0.5 * mCur[ al ][ be ] * tmp;
                        }
                    }
                    // al,3=3,al
                    for ( unsigned al=0; al<localDim; ++al ) {
                        double tmp =
                            dPhiDXi( al, l ) * dDU_ck[ d ] +
                            dDU_dl[ c ] * dPhiDXi( al, k );
                        la::addInnerProduct<dim>( tmp, aCur[ al ], dDUDU_dlck );
                        kgUU_dlck += 2.0 * 0.5 * nCur[ al ][ 2 ] * tmp;
                    }
                    // 3,3
                    ;  // do nothing
#ifdef GSHELL_KL_SYMMETRICSTIFFNESS
                    if ( k != l )  elemStiff( k*dof+c, l*dof+d ) += kgUU_dlck * detA * weight;  // positive LHS
#endif
                    elemStiff( l*dof+d, k*dof+c ) += kgUU_dlck * detA * weight;  // positive LHS
#endif
                }
            }
        }
    }

#ifdef GSHELL_KL_PRECOMPUTEVECTORSFORSTIFFNESS
    //--------------------------------------------------------------------------
    // remove pre-computed mid-plane vectors and their various derivatives
    delete [] aDU;
    delete [] bDU;
    delete [] pDU;
    delete [] pLGradDU;
    delete [] etaU;
    delete [] zetaU;
#endif

    // done
    return;
}

//------------------------------------------------------------------------------
// Compute strain resultants -- for output
template< typename BELEMENT, typename MAT >
corlib::eigenX::MatrixSd< 3, 3 >
gshell::fem::ElementKLStatic<BELEMENT,MAT>::giveStrainResultant(
    const enum plane planeMode, 
    const enum strains strainMode,
    const VecLDim & xi
    ) const
{
    // Collect shape functions and their derivatives
    VecNF phi;
    MatLDimNF dPhiDXi;
    MatSDimNF ddPhiDDXi;
    this->BasisElement::sfunGradHess( xi, phi, dPhiDXi, ddPhiDDXi );

    // get reference mid-vectors
    double aRef[dim][dim];       double bRef[dim][localDim][dim];
    double pRef[dim];            double pRefLGrad[localDim][dim];
    this->vectors_( gshell::fem::REFERENCE, xi, phi, dPhiDXi, ddPhiDDXi,
                    aRef, bRef, pRef, pRefLGrad );

    // get current mid-vectors
    double aCur[dim][dim];       double bCur[dim][localDim][dim];
    double pCur[dim];            double pLGrad[localDim][dim];
    this->vectors_( gshell::fem::CURRENT, xi, phi, dPhiDXi, ddPhiDDXi,
                    aCur, bCur, pCur, pLGrad );

    // Green-Lagrange strain resultants
    double alphaCur[dim][dim], betaCur[dim][dim];
    this->strainRes_( aRef, bRef[ 2 ], aCur, bCur[ 2 ],
                      alphaCur, betaCur );

    // convert to global frame
    double aInv[dim][dim];
    if      ( strainMode == LOCALGL ) {        // local Green-Lagrange
        la::identity<dim>( aInv );
    }
    else if ( strainMode == GREENLAGRANGE ) {  // global Green-Lagrange
        la::assign<dim,dim>( aInv, aRef );
        la::inverse<dim>( aInv );
    }
    else if ( strainMode == EULERALMANSI ) {   // global Euler-Almansi
        la::assign<dim,dim>( aInv, aCur );
        la::inverse<dim>( aInv );
    }
    else {                                     // error
        FTL_VERIFY_DESCRIPTIVE( false, "Unknown strain measure\n" );
    }
    double tmp[dim][dim];
    la::productNT<dim,dim,dim>( tmp, alphaCur, aInv );
    la::productNN<dim,dim,dim>( alphaCur, aInv, tmp );
    la::productNT<dim,dim,dim>( tmp, betaCur, aInv );
    la::productNN<dim,dim,dim>( betaCur, aInv, tmp );

    // return strain
    Mat3x3 strainResultant;
    if      ( planeMode == INPLANE ) {
        la::assign<dim,dim>( &(strainResultant(0,0)), alphaCur );
    }
    else if ( planeMode == OUTOFPLANE ) {
        la::assign<dim,dim>( &(strainResultant(0,0)), betaCur );
    }
    else {
        strainResultant.setZero( );
    }
    return strainResultant;
}

//------------------------------------------------------------------------------
//! Compute stress resultants -- for output
template< typename BELEMENT, typename MAT >
corlib::eigenX::MatrixSd< 3, 3 >
gshell::fem::ElementKLStatic<BELEMENT,MAT>::giveStressResultant(
    const enum plane planeMode, 
    const enum stresses stressMode,
    const VecLDim & xi
    ) const
{
    // Collect shape functions and their derivatives
    VecNF phi;
    MatLDimNF dPhiDXi;
    MatSDimNF ddPhiDDXi;
    this->BasisElement::sfunGradHess( xi, phi, dPhiDXi, ddPhiDDXi );

    // get reference mid-vectors
    double aRef[dim][dim];       double bRef[dim][localDim][dim];
    double pRef[dim];            double pRefLGrad[localDim][dim];
    this->vectors_( gshell::fem::REFERENCE, xi, phi, dPhiDXi, ddPhiDDXi,
                    aRef, bRef, pRef, pRefLGrad );

    // get mid-vectors
    double aCur[dim][dim];       double bCur[dim][localDim][dim];
    double pCur[dim];            double pLGrad[localDim][dim];
    this->vectors_( gshell::fem::CURRENT, xi, phi, dPhiDXi, ddPhiDDXi,
                    aCur, bCur, pCur, pLGrad );

    // Green-Lagrange strain resultants
    double alphaCur[dim][dim], betaCur[dim][dim];
    this->strainRes_( aRef, bRef[ 2 ], aCur, bCur[ 2 ],
                      alphaCur, betaCur );

    // 2nd PK stress resultants
    double nCur[dim][dim], mCur[dim][dim];
    this->stressResVK_( aRef, alphaCur, betaCur,
                        nCur, mCur );

    // convert to global frame
    double a[dim][dim];
    if      ( stressMode == LOCALPK2 ) {         // local 2nd Piola-Kirchhoff
        la::identity<dim>( a );
    }
    else if ( stressMode == PIOLAKIRCHHOFF2 ) {  // global 2nd Piola-Kirchhoff
        la::assign<dim,dim>( a, aRef );
    }
    else if ( stressMode == CAUCHY ) {           // global Cauchy (or true) stress
        la::assign<dim,dim>( a, aCur );
    }
    else {                                       // error
        FTL_VERIFY_DESCRIPTIVE( false, "Unknown stress measure\n" );
    }
    const double detA = la::determinant<dim>( a );
    double tmp[dim][dim];
    la::productNN<dim,dim,dim>( tmp, nCur, a );
    la::productTN<dim,dim,dim>( nCur, a, tmp );
    la::divide<dim,dim>( nCur, detA );
    la::productNN<dim,dim,dim>( tmp, mCur, a );
    la::productTN<dim,dim,dim>( mCur, a, tmp );
    la::divide<dim,dim>( mCur, detA );

    // return stress
    Mat3x3 stressResultant;
    if      ( planeMode == INPLANE ) {
        la::assign<dim,dim>( &(stressResultant(0,0)), nCur );
    }
    else if ( planeMode == OUTOFPLANE ) {
        la::assign<dim,dim>( &(stressResultant(0,0)), mCur );
    }
    else {
        stressResultant.setZero( );
    }
    return stressResultant;
}

