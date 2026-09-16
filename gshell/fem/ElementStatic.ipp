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
const unsigned gshell::fem::ElementStatic<BELEMENT,MAT>::
voigtForward[localDim][localDim] = {
    { BELEMENT::voigtForward[0][0], BELEMENT::voigtForward[0][1] },
    { BELEMENT::voigtForward[1][0], BELEMENT::voigtForward[1][1] }
};

//------------------------------------------------------------------------------
// Give general nodal quantity accessed by NodeFunctor op
template< typename BELEMENT, typename MAT >
template< typename OP >
void gshell::fem::ElementStatic<BELEMENT,MAT>::nodalQuantity_( MatDimNF & quanU,
                                                                MatDimNF & quanW,
                                                                OP op ) const
{
    const unsigned numNodes = this->numFunctions( );
    quanU.conservativeResize( dim, numNodes );
    quanW.conservativeResize( dim, numNodes );
    for ( unsigned j = 0; j < numNodes; j ++ ) {
        const VecDof qVec = op( this->BasisElement::giveSupportNodePtr( j ) );
        for ( unsigned d=0; d<dim; ++d )
            quanU( d, j ) = qVec[ d ];
        for ( unsigned d=0; d<localDim; ++d )
            quanW( d, j ) = qVec[ dim+d ];
        for ( unsigned d=localDim; d<dim; ++d )
            quanW( d, j ) = 0.0;
    }
    return;
}

//------------------------------------------------------------------------------
// Give dof indices of nodes
template< typename BELEMENT, typename MAT >
void gshell::fem::ElementStatic<BELEMENT,MAT>::
getDofIndices( std::vector<unsigned> & dofIndices ) const
{
    const unsigned numNodes = this->numFunctions( );
    const unsigned matSize = dof * numNodes;
    dofIndices.resize( matSize );
    std::vector<unsigned>::iterator iter = dofIndices.begin();
    for ( unsigned v = 0; v < numNodes; v ++ ) { 
        std::vector<unsigned> aux;
        this->BasisElement::giveSupportNodePtr( v ) -> copyDofArray( aux );
        iter = std::copy( aux.begin(), aux.end(), iter );
    }
    return;
}

//------------------------------------------------------------------------------
// Return nodal displacements
template< typename BELEMENT, typename MAT >
void gshell::fem::ElementStatic<BELEMENT,MAT>::nodalDisplacements(
    MatDimNF & u
    ) const
{
    MatDimNF w;  // dummy shears
    this->nodalDisplacementsAndShears( u, w );
}

//------------------------------------------------------------------------------
// Return nodal displacements
template< typename BELEMENT, typename MAT >
void gshell::fem::ElementStatic<BELEMENT,MAT>::nodalDisplacementsAndShears(
    MatDimNF & u,
    MatDimNF & w
    ) const
{
    this->nodalQuantity_( u, w, std::mem_fun( &Node::giveDisplacements ) );

    const unsigned nf = this->numFunctions( );
    FTL_VERIFY_DESCRIPTIVE( u.rows( ) == dim,
                            "u.rows( )=%d unequal dim=%d\n", u.rows( ), dim );
    FTL_VERIFY_DESCRIPTIVE( u.cols( ) == nf,
                            "u.cols( )=%d unequal nf=%d\n",  u.cols( ), nf  );
    FTL_VERIFY_DESCRIPTIVE( w.rows( ) == dim,
                            "w.rows( )=%d unequal dim=%d\n", w.rows( ), dim );
    FTL_VERIFY_DESCRIPTIVE( w.cols( ) == nf,
                            "w.cols( )=%d unequal nf=%d\n",  w.cols( ), nf  );

    return;
}

//------------------------------------------------------------------------------
// Return nodal increments
template< typename BELEMENT, typename MAT >
void gshell::fem::ElementStatic<BELEMENT,MAT>::nodalIncrements(
    MatDimNF & deltaU,
    MatDimNF & deltaW
    ) const
{
    this->nodalQuantity_( deltaU, deltaW, std::mem_fun( &Node::giveIncrement ) );
    return;
}

//------------------------------------------------------------------------------
// Interpolate current displacement at local point
template< typename BELEMENT, typename MAT >
corlib::eigenX::VectorSd< BELEMENT::dim > 
gshell::fem::ElementStatic<BELEMENT,MAT>::giveDisplacement(
    const VecLDim & xi
    ) const
{
    // check xi
    FTL_VERIFY( corlib::ShapeTraits< myShape >::isInside( xi ) );
    // shape functions at xi
    VecNF phi;          this->BasisElement::sfun( xi, phi );
    // nodal displacements (and shears)
    MatDimNF uNF;       this->nodalDisplacements( uNF );
    // interpolate
    const MatDimNF uNFT = uNF.transpose();
    VecDim u;           SurfaceGeometer::interpolate( this->numFunctions(), &(phi(0)),
                                                      &(uNFT(0,0)), &(u(0)) );
    return u;
}

//------------------------------------------------------------------------------
// Interpolate current displacement at local point
template< typename BELEMENT, typename MAT >
corlib::eigenX::VectorSd< BELEMENT::dim > 
gshell::fem::ElementStatic<BELEMENT,MAT>::giveDisplacementAtNode(
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
void gshell::fem::ElementStatic<BELEMENT,MAT>::computeTangentsAtNode(
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
void gshell::fem::ElementStatic<BELEMENT,MAT>::computeNormalAtNode(
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
                const VecDim tang1DU_dl = ( *tang1GradDof).col( l*dof+d );

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
// Transform convected nodal shear vector to nodal Cartesian shear vector
template< typename BELEMENT, typename MAT >
void gshell::fem::ElementStatic<BELEMENT,MAT>::nodalCartesianShearsCur_(
    MatDimNF & wNF
    ) const
{
    const unsigned numFunctions = this->numFunctions( );
    for ( unsigned k=0; k<numFunctions; ++k ) {
        const Node * node_k = this->BasisElement::giveSupportNodePtr( k );

        std::array< VecDim, localDim > tang_k;
        tang_k[ 0 ] = node_k->getTangent( gshell::fem::CURRENT, 0 );
        tang_k[ 1 ] = node_k->getTangent( gshell::fem::CURRENT, 1 );

        // Cartesify shear vector
        VecDim wCart; wCart.setZero();
        for ( unsigned d=0; d<localDim; ++d )
            wCart += wNF( d, k ) * tang_k[ d ];
        for ( unsigned d=0; d<dim; ++d )
            wNF( d, k ) = wCart[ d ];
    }
    return;
}

//------------------------------------------------------------------------------
// Determine coefficients (derivatives) of tangents at element nodes w.r.t. node #l
template< typename BELEMENT, typename MAT >
void gshell::fem::ElementStatic<BELEMENT,MAT>::tangentsGradCur_(
    const unsigned & l,
    MatLDimNF & tangCoeff_l // L_{ga I}^{#l}
    ) const
{
    const unsigned numFunctions = this->numFunctions( );
    tangCoeff_l.conservativeResize( localDim, numFunctions );
    Node * node_l = this->BasisElement::giveSupportNodePtr( l );
    for ( unsigned v=0; v<numFunctions; ++v ) {
        const Node * node_v = this->BasisElement::giveSupportNodePtr( v );
        for ( unsigned ga=0; ga<localDim; ++ga ) {
            tangCoeff_l( ga, v ) = node_v->getTangentCoefficient( gshell::fem::CURRENT, node_l, ga );
        }
    }
    return;
}

//------------------------------------------------------------------------------
// Resize a double-type matrix of size (rows x NF) and optionally set to zero
template< typename BELEMENT, typename MAT >
template< typename ENTRY >
void gshell::fem::ElementStatic<BELEMENT,MAT>::resizeMatrixNF_(
    Eigen::Matrix<ENTRY, Eigen::Dynamic, Eigen::Dynamic> & mat,
    unsigned rows,
    bool setToZero
    ) const
{
    const unsigned numFunctions = this->numFunctions( );
    mat.resize( rows, numFunctions );
    if ( setToZero )
        for ( unsigned a=0; a<rows; ++a )
            for ( unsigned k=0; k<numFunctions; ++k )
                mat( a, k ) = corlib::detail_::AccessorTraits< ENTRY >::zero( );
    return;
}


//------------------------------------------------------------------------------
// Write a double-type matrix of size (rows x NF) to stream
template< typename BELEMENT, typename MAT >
template< typename ENTRY >
std::ostream & gshell::fem::ElementStatic<BELEMENT,MAT>::writeMatrixNF_(
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
// Compute mid-surface vectors and their derivatives w.r.t. nodal DOFs
template< typename BELEMENT, typename MAT >
void gshell::fem::ElementStatic<BELEMENT,MAT>::vectors_(
    const enum gshell::fem::config conf,
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
    double cCur[localDim][dim]  // dLGrad
    ) const
{
    const unsigned numFunctions = this->numFunctions( );

    // Collect nodal quantities
    MatDimNF xNF;  this->BasisElement::supportNodeCoordinates( xNF );
    MatDimNF uNF, wNF;
    if ( conf == gshell::fem::CURRENT ) {
        this->nodalDisplacementsAndShears( uNF, wNF );
        xNF += uNF;
    }
    else {
        // this IMPLIES that the INITIAL DISPLACEMENTS and SHEARS are ZERO
        //this->resizeMatrixNF_( uNF, dim, true );
        this->resizeMatrixNF_( wNF, dim, true );
    }

    // Compute base
    const MatLDimNF dPhiDXiT = dPhiDXi.transpose();
    const MatDimNF xNFT = xNF.transpose();
    SurfaceGeometer::covariantBase( numFunctions, &(dPhiDXiT(0,0)), &(xNFT(0,0)), aCur );
    // co-variant base vectors
    const double (&a0Cur)[dim] = aCur[ 0 ];
    const double (&a1Cur)[dim] = aCur[ 1 ];
    //const double (&a2Cur)[dim] = aCur[ 2 ];
    // Local gradient of base vectors
    const MatSDimNF ddPhiDDXiT = ddPhiDDXi.transpose();
    SurfaceGeometer::covariantBaseDeriv1( numFunctions, &(ddPhiDDXiT(0,0)), &(xNFT(0,0)), aCur, bCur );
    const double (&b0Cur)[localDim][dim] = bCur[ 0 ];
    const double (&b1Cur)[localDim][dim] = bCur[ 1 ];

    // Transform convected nodal shear to nodal Cartesian shear vector
    if ( conf == gshell::fem::CURRENT )
        this->nodalCartesianShearsCur_( wNF );
    // Shear vector field
    const MatDimNF wNFT = wNF.transpose();
    SurfaceGeometer::interpolate( numFunctions, &(phi(0)), &(wNFT(0,0)), wCur );
    // Local derivative of shear vector field
    SurfaceGeometer::interpolateDeriv1( numFunctions, &(dPhiDXiT(0,0)), &(wNFT(0,0)), zCur );

    // \vec{p} = \vec{a}_1 \times \vec{a}_2
    la::crossProduct( pCur, a0Cur, a1Cur );
    SurfaceGeometer::crossprodVectorDeriv1( a0Cur, b0Cur, a1Cur, b1Cur, pLGrad );

    // q^A = p^A + w^A
    la::assign<dim>( qCur, pCur );
    la::add<dim>( qCur, wCur );
    double qLen;  la::norm2<dim>( qLen, qCur );
    la::assign<localDim,dim>( qLGrad, pLGrad );
    la::add<localDim,dim>( qLGrad, zCur );
    
    // Thickness director
    SurfaceGeometer::director( pCur, wCur, dCur );
    SurfaceGeometer::unitVectorDeriv1( qLen, qLGrad, dCur, cCur );
    // done
    return;
}

//------------------------------------------------------------------------------
// Compute 1st derivatives of mid-surface vectors w.r.t. nodal DOFs
template< typename BELEMENT, typename MAT >
void gshell::fem::ElementStatic<BELEMENT,MAT>::vectorsGradCur_(
    const unsigned & d,
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
    double cDW_dl[localDim][dim]
    ) const
{
    const unsigned numFunctions = this->numFunctions( );

    // co-variant base vectors
    const double (&a0Cur)[dim] = aCur[ 0 ];
    const double (&a1Cur)[dim] = aCur[ 1 ];
    //const double (&a2Cur)[dim] = aCur[ 2 ];

    // 1st derivative of \vec{a}_\alpha w.r.t. nodal DOFs
    double a0DU_dl[dim];  la::zero<dim>( a0DU_dl );  a0DU_dl[ d ] = dPhiDXi( 0, l );
    double a1DU_dl[dim];  la::zero<dim>( a1DU_dl );  a1DU_dl[ d ] = dPhiDXi( 1, l );

    // local gradient of co-variant base vectors \vec{a}_{\alpha,\beta}
    const double (&b0Cur)[localDim][dim] = bCur[ 0 ];
    const double (&b1Cur)[localDim][dim] = bCur[ 1 ];

    // 1st derivative of local gradient of co-variant base vectors \vec{a}_{\alpha,\beta}
    // w.r.t. nodal DOFs
    double b0DU_dl[localDim][dim];  la::zero<localDim,dim>( b0DU_dl );
    double b1DU_dl[localDim][dim];  la::zero<localDim,dim>( b1DU_dl );
    for ( unsigned be=0; be<localDim; ++be ) {
        b0DU_dl[ be ][ d ] = ddPhiDDXi( voigtForward[ 0 ][ be ], l );
        b1DU_dl[ be ][ d ] = ddPhiDDXi( voigtForward[ 1 ][ be ], l );
    }

    // nodal displacements and shears (convected components)
    MatDimNF wNF, uNF;  this->nodalDisplacementsAndShears( uNF, wNF );

    // compute tangents at configuration
    std::array< VecDim, localDim > tang_l;
    Node * node_l = this->BasisElement::giveSupportNodePtr( l );
    for ( unsigned ga=0; ga<localDim; ++ga ) {
        tang_l[ ga ] = node_l->getTangent( gshell::fem::CURRENT, ga );
    }

    // 1st derivative of shear vectors w.r.t. nodal displacements
    double wDU_dl[dim];  la::zero<dim>( wDU_dl );
    for ( unsigned v=0; v<numFunctions; ++v ) {
        for ( unsigned al=0; al<localDim; ++al ) {
            wDU_dl[ d ] += wNF( al, v ) * phi( v ) * tangCoeff_l( al, v );
        }
    }

    double zDU_dl[localDim][dim];  la::zero<localDim,dim>( zDU_dl );
    for ( unsigned v=0; v<numFunctions; ++v ) {
        for ( unsigned be=0; be<localDim; ++be ) {
            for ( unsigned al=0; al<localDim; ++al ) {
                zDU_dl[ be ][ d ] += wNF( al, v ) * dPhiDXi( be, v ) * tangCoeff_l( al, v );
            }
        }
    }

    // 1st derivative of shear vectors w.r.t. nodal shears
    double wDW_dl[dim];
    if ( d < localDim ) {
        la::assign<dim>(   wDW_dl, &(tang_l[ d ][ 0 ]) );
        la::multiply<dim>( wDW_dl, phi( l ) );
    }
    else {
        la::zero<dim>( wDW_dl );
    }
    
    double zDW_dl[localDim][dim];
    for ( unsigned be=0; be<localDim; ++be ) {
        if ( d < localDim ) {
            la::assign<dim>(   zDW_dl[ be ], &(tang_l[ d ][ 0 ]) );
            la::multiply<dim>( zDW_dl[ be ], dPhiDXi( be, l ) );
        }
        else {
            la::zero<localDim,dim>( zDW_dl );
        }
    }

    // \vec{p} = \vec{a}_1 \times \vec{a}_2
    // P^A{}_D{}^L = \frac{ \partial p^A }{ \partial u^D{}_L }
    SurfaceGeometer::crossprodVectorGradDof( a0Cur, a0DU_dl, a1Cur, a1DU_dl,
                                             pDU_dl );
    // Q^A{}_\beta{}_D{}^L = \frac{ \partial p^A{}_{,\beta} }{ \partial u^D{}_L }
    SurfaceGeometer::crossprodVectorDeriv1GradDof( a0Cur, a0DU_dl, b0Cur, b0DU_dl,
                                                   a1Cur, a1DU_dl, b1Cur, b1DU_dl,
                                                   pLGradDU_dl );

    // q^A = p^A + w^A
    double qLen;  la::norm2<dim>( qLen, qCur );
    la::assign<dim>( qDU_dl, pDU_dl );
    la::add<dim>(    qDU_dl, wDU_dl );
    la::assign<dim>( qDW_dl, wDW_dl );
    la::assign<localDim,dim>( qLGradDU_dl, pLGradDU_dl ); 
    la::add<localDim,dim>(    qLGradDU_dl, zDU_dl );
    la::assign<localDim,dim>( qLGradDW_dl, zDW_dl );

    // s^A{}_{\beta} = q^A{}_{,\beta} / |q^A|
    double sCur[localDim][dim];
    la::assign<localDim,dim>( sCur, qLGrad );
    la::divide<localDim,dim>( sCur, qLen );

    // DU^D{}_D{}^L = \frac{ \partial d^A }{ \partial u^D{}_L }
    // DW^A{}_D{}^L = \frac{ \partial d^A }{ \partial w^D{}_L }
    SurfaceGeometer::unitVectorGradDof( qLen, qDU_dl, dCur, dDU_dl );
    SurfaceGeometer::unitVectorGradDof( qLen, qDW_dl, dCur, dDW_dl );

    // CU^A{}_\beta{}_D{}^L = \frac{ \partial c^A{}_\beta }{ \partial u^D{}_L }
    // CW^A{}_\beta{}_D{}^L = \frac{ \partial c^A{}_\beta }{ \partial w^D{}_L }
    SurfaceGeometer::unitVectorDeriv1GradDof( qLen, qDU_dl, qLGradDU_dl, dCur, dDU_dl, sCur, cDU_dl );
    SurfaceGeometer::unitVectorDeriv1GradDof( qLen, qDW_dl, qLGradDW_dl, dCur, dDW_dl, sCur, cDW_dl );
    // done
    return;
}

//------------------------------------------------------------------------------
// Computes the Green-Lagrange strain resultants
template< typename BELEMENT, typename MAT >
void gshell::fem::ElementStatic<BELEMENT,MAT>::strainRes_(
    const double aRef[dim][dim],
    const double dRef[dim],
    const double cRef[localDim][dim],
    const double aCur[dim][dim],
    const double dCur[dim],
    const double cCur[localDim][dim],
    double alpha[dim][dim],
    double beta[dim][dim]
    ) const
{

    // membrane strain resultants
    la::zero<dim,dim>( alpha );
    // al,be=be,al
    la::innerProduct<dim>(         alpha[ 0 ][ 0 ], aCur[ 0 ], aCur[ 0 ] );
    la::subtractInnerProduct<dim>( alpha[ 0 ][ 0 ], aRef[ 0 ], aRef[ 0 ] );
    la::innerProduct<dim>(         alpha[ 0 ][ 1 ], aCur[ 0 ], aCur[ 1 ] );
    la::subtractInnerProduct<dim>( alpha[ 0 ][ 1 ], aRef[ 0 ], aRef[ 1 ] );
    la::innerProduct<dim>(         alpha[ 1 ][ 0 ], aCur[ 1 ], aCur[ 0 ] );
    la::subtractInnerProduct<dim>( alpha[ 1 ][ 0 ], aRef[ 1 ], aRef[ 0 ] );
    la::innerProduct<dim>(         alpha[ 1 ][ 1 ], aCur[ 1 ], aCur[ 1 ] );
    la::subtractInnerProduct<dim>( alpha[ 1 ][ 1 ], aRef[ 1 ], aRef[ 1 ] );
    // al,3=3,al
    la::innerProduct<dim>(         alpha[ 0 ][ 2 ], aCur[ 0 ], dCur );
    la::subtractInnerProduct<dim>( alpha[ 0 ][ 2 ], aRef[ 0 ], dRef );
    la::innerProduct<dim>(         alpha[ 1 ][ 2 ], aCur[ 1 ], dCur );
    la::subtractInnerProduct<dim>( alpha[ 1 ][ 2 ], aRef[ 1 ], dRef );
    la::innerProduct<dim>(         alpha[ 2 ][ 0 ], aCur[ 0 ], dCur );
    la::subtractInnerProduct<dim>( alpha[ 2 ][ 0 ], aRef[ 0 ], dRef );
    la::innerProduct<dim>(         alpha[ 2 ][ 1 ], aCur[ 1 ], dCur );
    la::subtractInnerProduct<dim>( alpha[ 2 ][ 1 ], aRef[ 1 ], dRef );
    // 3,3 is zero
    la::multiply<dim,dim>( alpha, 0.5 );

    // bending strain resultants
    la::zero<dim,dim>( beta );
    // remaining entries are zero
    for ( unsigned al = 0; al < localDim; ++al ) {
        for ( unsigned be = 0; be < localDim; ++be ) {
            double beta_albe = 0.;
            la::addInnerProduct<dim>(      beta_albe, aCur[ al ], cCur[ be ] );
            la::addInnerProduct<dim>(      beta_albe, cCur[ al ], aCur[ be ] );
            la::subtractInnerProduct<dim>( beta_albe, aRef[ al ], cRef[ be ] );
            la::subtractInnerProduct<dim>( beta_albe, cRef[ al ], aRef[ be ] );
            beta[ al ][ be ] = 0.5 * beta_albe;
        }
    }
    return;
}

//------------------------------------------------------------------------------
// Computes Green-Lagrange strain resultants and their derivatives w.r.t. nodal DOFs
template< typename BELEMENT, typename MAT >
void gshell::fem::ElementStatic<BELEMENT,MAT>::strainResGradCur_(
    const unsigned & d,
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
    double betaDW_dl[dim][dim]
    ) const
{
    la::zero<dim,dim>( alphaDU_dl );
    la::zero<dim,dim>( alphaDW_dl );
    la::zero<dim,dim>( betaDU_dl );
    la::zero<dim,dim>( betaDW_dl );

    // al,be
    for ( unsigned al=0; al<localDim; ++al ) {
        for ( unsigned be=0; be<localDim; ++be ) {
            alphaDU_dl[ al ][ be ] = 0.5 * (
                aCur[ be ][ d ] * dPhiDXi( al, l ) +
                aCur[ al ][ d ] * dPhiDXi( be, l )
                );

            betaDU_dl[ al ][ be ] = 
                cCur[ be ][ d ] * dPhiDXi( al, l ) +
                cCur[ al ][ d ] * dPhiDXi( be, l );
            la::addInnerProduct<dim>( betaDU_dl[ al ][ be ], aCur[ al ], cDU_dl[ be ] );
            la::addInnerProduct<dim>( betaDU_dl[ al ][ be ], aCur[ be ], cDU_dl[ al ] );
            betaDU_dl[ al ][ be ] *= 0.5;

            betaDW_dl[ al ][ be ] = 0.0;
            la::addInnerProduct<dim>( betaDW_dl[ al ][ be ], aCur[ al ], cDW_dl[ be ] );
            la::addInnerProduct<dim>( betaDW_dl[ al ][ be ], aCur[ be ], cDW_dl[ al ] );
            betaDW_dl[ al ][ be ] *= 0.5;
        }
    }
    // al,3=3,al
    for ( unsigned al=0; al<localDim; ++al ) {
        double alphaDU_dlal2 = dCur[ d ] * dPhiDXi( al, l );
        la::addInnerProduct<dim>( alphaDU_dlal2, aCur[ al ], dDU_dl );
        alphaDU_dlal2 *= 0.5;
        alphaDU_dl[ al ][ 2 ] = alphaDU_dlal2;
        alphaDU_dl[ 2 ][ al ] = alphaDU_dlal2;

        double alphaDW_dlal2;
        la::innerProduct<dim>( alphaDW_dlal2, aCur[ al ], dDW_dl );
        alphaDW_dlal2 *= 0.5;
        alphaDW_dl[ al ][ 2 ] = alphaDW_dlal2;
        alphaDW_dl[ 2 ][ al ] = alphaDW_dlal2;
#if 1
        // omitted --- zero anyway
        double betaDU_dlal2;
        la::innerProduct<dim>( betaDU_dlal2, dCur, cDU_dl[ al ] );
        betaDU_dlal2 *= 0.5;
        betaDU_dl[ al ][ 2 ] = betaDU_dlal2;
        betaDU_dl[ 2 ][ al ] = betaDU_dlal2;

        double betaDW_dlal2;
        la::innerProduct<dim>( betaDW_dlal2, dCur, cDW_dl[ al ] );
        betaDW_dlal2 *= 0.5;
        betaDW_dl[ al ][ 2 ] = betaDW_dlal2;
        betaDW_dl[ 2 ][ al ] = betaDW_dlal2;
#endif
    }
#if 1
    // 3,3 --- zero anyway
    la::innerProduct<dim>( alphaDU_dl[ 2 ][ 2 ], dCur, dDU_dl );
    la::innerProduct<dim>( alphaDW_dl[ 2 ][ 2 ], dCur, dDW_dl );
#endif
    // done
    return;
}


//------------------------------------------------------------------------------
// Computes the elasticity 4-tensor of Saint Venant-Kirchhoff material
template< typename BELEMENT, typename MAT >
void gshell::fem::ElementStatic<BELEMENT,MAT>::elasticityTensorVK_(
    const double aRef[dim][dim],
    double hRef[dim][dim][dim][dim]
    ) const
{
    // material constants 
    const double E = this->giveYoungsModulus( );
    const double nu = this->givePoissonRatio( );
    const double lam = nu*E/(1.+nu)/(1.-2.*nu);
    const double mu = .5*E/(1.+nu);
    //const double fact = 2.*mu*lam/(2.*mu+lam);  // due to "elmination"
    const double fact = lam;  // due to "application"

    // shear correction factor
    const double kappa = 5./6.;

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
// Computes 2nd Piola-Kirchhoff stress resultants of Venant-Kirchhoff material
template< typename BELEMENT, typename MAT >
void gshell::fem::ElementStatic<BELEMENT,MAT>::stressResVK_(
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
void gshell::fem::ElementStatic<BELEMENT,MAT>::bodyForce(
    const VecLDim & xi,
    const double & weight,
    FUNC func, 
    const double & factor
    )
{
    // constants
    const unsigned numFunctions = this->numFunctions( );
    double zeroDim[dim];  la::zero<dim>( zeroDim );

    // material point
    const VecDim xRef = this->giveCoordinate( xi );

    // force density
    const VecDim forceDensityConst = factor * func( xRef );

    // collect shape functions and their derivatives
    VecNF phi;              //this->BasisElement::sfun( xi, phi );
    MatLDimNF dPhiDXi;      //this->BasisElement::sfunGrad( xi, dPhiDXi );
    MatSDimNF ddPhiDDXi;    //this->BasisElement::sfunHess( xi, ddPhiDDXi );
    this->BasisElement::sfunGradHess( xi, phi, dPhiDXi, ddPhiDDXi );
    // get reference mid-vectors
    double aRef[dim][dim];        double bRef[dim][localDim][dim];
    double wRef[dim];             double zRef[localDim][dim];
    double pRef[dim];             double pRefLGrad[localDim][dim];
    double qRef[dim];             double qRefLGrad[localDim][dim];
    double dRef[dim];             double cRef[localDim][dim];
    this->vectors_( gshell::fem::REFERENCE, xi, phi, dPhiDXi, ddPhiDDXi,
                    aRef, bRef, wRef, zRef,
                    pRef, pRefLGrad, qRef, qRefLGrad, dRef, cRef );

    // get current mid-vectors
    double aCur[dim][dim];        double bCur[dim][localDim][dim];
    double wCur[dim];             double zCur[localDim][dim];
    double pCur[dim];             double pLGrad[localDim][dim];
    double qCur[dim];             double qLGrad[localDim][dim];
    double dCur[dim];             double cCur[localDim][dim];
    this->vectors_( gshell::fem::CURRENT, xi, phi, dPhiDXi, ddPhiDDXi,
                    aCur, bCur, wCur, zCur,
                    pCur, pLGrad, qCur, qLGrad, dCur, cCur );

    // convenience
//    const double (&a0Cur)[dim] = aCur[ 0 ];
//    const double (&a1Cur)[dim] = aCur[ 1 ];
    //const double (&a2Cur)[dim] = aCur[ 2 ];
    double qLen;  la::norm2<dim>( qLen, qCur );
    //const double (&b0Cur)[dim][dim] = bCur[ 0 ];
    //const double (&b1Cur)[dim][dim] = bCur[ 1 ];

    // surface metric
    const double detA = la::determinant<dim>( aRef );

    // thickness
    const double & t = thickness_;
    // second moment of area
//    const double s = t*t*t / 12.;

    // pass each nodal force to the nodes
    for ( unsigned l=0; l<numFunctions; ++l ) {
        VecDof nodalForce;
        nodalForce.setZero( );
        for ( unsigned d=0; d<dim; ++d ) {
            // 1st derivatives of vectors
            double uDU_dl[dim];  la::zero<dim>( uDU_dl );  uDU_dl[ d ] = phi( l );

            // 1st derivative of tangents
            MatLDimNF tangCoeff_l;
            this->tangentsGradCur_( l, tangCoeff_l );

            // nodal force
            double tmp;
            la::innerProduct<dim>( tmp, &(forceDensityConst(0)), uDU_dl );
            nodalForce( d ) = t * tmp * detA * weight;
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
void gshell::fem::ElementStatic<BELEMENT,MAT>::bodyForceLin(
    const VecLDim & xi,
    const double & weight,
    FUNC func, 
    const double & factor
    )
{
    // constants
    const unsigned numFunctions = this->numFunctions( );
    double zeroDim[dim];  la::zero<dim>( zeroDim );

    // material point
    const VecDim xRef = this->giveCoordinate( xi );

    // force density
    const VecDim forceDensityLin   = factor * func( xRef );

    // collect shape functions and their derivatives
    VecNF phi;              //this->BasisElement::sfun( xi, phi );
    MatLDimNF dPhiDXi;      //this->BasisElement::sfunGrad( xi, dPhiDXi );
    MatSDimNF ddPhiDDXi;    //this->BasisElement::sfunHess( xi, ddPhiDDXi );
    this->BasisElement::sfunGradHess( xi, phi, dPhiDXi, ddPhiDDXi );

    // get reference mid-vectors
    double aRef[dim][dim];        double bRef[dim][localDim][dim];
    double wRef[dim];             double zRef[localDim][dim];
    double pRef[dim];             double pRefLGrad[localDim][dim];
    double qRef[dim];             double qRefLGrad[localDim][dim];
    double dRef[dim];             double cRef[localDim][dim];
    this->vectors_( gshell::fem::REFERENCE, xi, phi, dPhiDXi, ddPhiDDXi,
                    aRef, bRef, wRef, zRef,
                    pRef, pRefLGrad, qRef, qRefLGrad, dRef, cRef );

    // get current mid-vectors
    double aCur[dim][dim];        double bCur[dim][localDim][dim];
    double wCur[dim];             double zCur[localDim][dim];
    double pCur[dim];             double pLGrad[localDim][dim];
    double qCur[dim];             double qLGrad[localDim][dim];
    double dCur[dim];             double cCur[localDim][dim];
    this->vectors_( gshell::fem::CURRENT, xi, phi, dPhiDXi, ddPhiDDXi,
                    aCur, bCur, wCur, zCur,
                    pCur, pLGrad, qCur, qLGrad, dCur, cCur );

    // convenience
    const double (&a0Cur)[dim] = aCur[ 0 ];
    const double (&a1Cur)[dim] = aCur[ 1 ];
    //const double (&a2Cur)[dim] = aCur[ 2 ];
    double qLen;  la::norm2<dim>( qLen, qCur );
    //const double (&b0Cur)[dim][dim] = bCur[ 0 ];
    //const double (&b1Cur)[dim][dim] = bCur[ 1 ];

    // surface metric
    const double detA = la::determinant<dim>( aRef );

    // thickness
    const double & t = thickness_;
    // second moment of area
    const double s = t*t*t / 12.;

    // pass each nodal force to the nodes
    for ( unsigned l=0; l<numFunctions; ++l ) {
        VecDof nodalForce;
        nodalForce.setZero( );
        for ( unsigned d=0; d<dim; ++d ) {
            // 1st derivatives of vectors
            double uDU_dl[dim];  la::zero<dim>( uDU_dl );  uDU_dl[ d ] = phi( l );

            // 1st derivative of tangents
            MatLDimNF tangCoeff_l;
            this->tangentsGradCur_( l, tangCoeff_l );

            // nodal force
            double tmp;
            nodalForce.setZero();

            // compute nodal forces
            {
                // 1st derivatives of vectors
                double pDU_dl[dim];           double pLGradDU_dl[localDim][dim];
                double qDU_dl[dim];           double qLGradDU_dl[localDim][dim];
                double qDW_dl[dim];           double qLGradDW_dl[localDim][dim];
                double dDU_dl[dim];           double cDU_dl[localDim][dim];
                double dDW_dl[dim];           double cDW_dl[localDim][dim];
                this->vectorsGradCur_( d, l, xi, phi, dPhiDXi, ddPhiDDXi, tangCoeff_l,
                                       aRef, aCur, bCur, wCur, zCur, 
                                       pCur, pLGrad, qCur, qLGrad, dCur, cCur,
                                       pDU_dl, pLGradDU_dl, qDU_dl, qDW_dl,
                                       qLGradDU_dl, qLGradDW_dl,
                                       dDU_dl, dDW_dl, cDU_dl, cDW_dl );

                // nodal force
                la::innerProduct<dim>( tmp, &(forceDensityLin(0)), dDU_dl );
                nodalForce( d ) += s * tmp * detA * weight;
                if ( d < localDim ) {
                    la::innerProduct<dim>( tmp, &(forceDensityLin(0)), dDW_dl );
                    nodalForce( dim + d ) = s * tmp * detA * weight;
                }

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
void gshell::fem::ElementStatic<BELEMENT,MAT>::bodyForceLinStiffness(
    const VecLDim & xi,
    const double & weight,
    FUNC func, 
    const double & factor,
    Eigen::MatrixXd & elemStiff
    )
{
    // constants
    const unsigned numFunctions = this->numFunctions( );
    const unsigned matSize = dof * numFunctions;
    FTL_VERIFY( elemStiff.rows() == matSize and elemStiff.cols() == matSize );
    double zeroDim[dim];  la::zero<dim>( zeroDim );

    // material point
    const VecDim xRef = this->giveCoordinate( xi );

    // force density
    const VecDim forceDensityLin   = factor * func( xRef );

    // collect shape functions and their derivatives
    VecNF phi;              //this->BasisElement::sfun( xi, phi );
    MatLDimNF dPhiDXi;      //this->BasisElement::sfunGrad( xi, dPhiDXi );
    MatSDimNF ddPhiDDXi;    //this->BasisElement::sfunHess( xi, ddPhiDDXi );
    this->BasisElement::sfunGradHess( xi, phi, dPhiDXi, ddPhiDDXi );

    // get reference mid-vectors
    double aRef[dim][dim];        double bRef[dim][localDim][dim];
    double wRef[dim];             double zRef[localDim][dim];
    double pRef[dim];             double pRefLGrad[localDim][dim];
    double qRef[dim];             double qRefLGrad[localDim][dim];
    double dRef[dim];             double cRef[localDim][dim];
    this->vectors_( gshell::fem::REFERENCE, xi, phi, dPhiDXi, ddPhiDDXi,
                    aRef, bRef, wRef, zRef,
                    pRef, pRefLGrad, qRef, qRefLGrad, dRef, cRef );

    // get current mid-vectors
    double aCur[dim][dim];        double bCur[dim][localDim][dim];
    double wCur[dim];             double zCur[localDim][dim];
    double pCur[dim];             double pLGrad[localDim][dim];
    double qCur[dim];             double qLGrad[localDim][dim];
    double dCur[dim];             double cCur[localDim][dim];
    this->vectors_( gshell::fem::CURRENT, xi, phi, dPhiDXi, ddPhiDDXi,
                    aCur, bCur, wCur, zCur,
                    pCur, pLGrad, qCur, qLGrad, dCur, cCur );

    // convenience
    const double (&a0Cur)[dim] = aCur[ 0 ];
    const double (&a1Cur)[dim] = aCur[ 1 ];
    //const double (&a2Cur)[dim] = aCur[ 2 ];
    double qLen;  la::norm2<dim>( qLen, qCur );
    //const double (&b0Cur)[dim][dim] = bCur[ 0 ];
    //const double (&b1Cur)[dim][dim] = bCur[ 1 ];

    // surface metric
    const double detA = la::determinant<dim>( aRef );

    // thickness
    const double & t = thickness_;
    // second moment of area
    const double s = t*t*t / 12.;

    // pass each nodal force to the nodes
    for ( unsigned l=0; l<numFunctions; ++l ) {
        VecDof nodalForce;
        nodalForce.setZero( );
        for ( unsigned d=0; d<dim; ++d ) {
            // 1st derivatives of vectors
            double uDU_dl[dim];  la::zero<dim>( uDU_dl );  uDU_dl[ d ] = phi( l );

            // 1st derivative of tangents
            MatLDimNF tangCoeff_l;
            this->tangentsGradCur_( l, tangCoeff_l );

            // // nodal force
            // double tmp;
            // la::innerProduct<dim>( tmp, &(forceDensityConst(0)), uDU_dl );
            // nodalForce( d ) = t * tmp * detA * weight;

            // external stiffness matrix
            {

                // 1st derivatives of vectors
                double pDU_dl[dim];           double pLGradDU_dl[localDim][dim];
                double qDU_dl[dim];           double qLGradDU_dl[localDim][dim];
                double qDW_dl[dim];           double qLGradDW_dl[localDim][dim];
                double dDU_dl[dim];           double cDU_dl[localDim][dim];
                double dDW_dl[dim];           double cDW_dl[localDim][dim];
                this->vectorsGradCur_( d, l, xi, phi, dPhiDXi, ddPhiDDXi, tangCoeff_l,
                                       aRef, aCur, bCur, wCur, zCur, 
                                       pCur, pLGrad, qCur, qLGrad, dCur, cCur,
                                       pDU_dl, pLGradDU_dl, qDU_dl, qDW_dl,
                                       qLGradDU_dl, qLGradDW_dl,
                                       dDU_dl, dDW_dl, cDU_dl, cDW_dl );

                // // nodal force
                // la::innerProduct<dim>( tmp, &(forceDensityLin(0)), dDU_dl );
                // nodalForce( d ) += s * tmp * detA * weight;
                // if ( d < localDim ) {
                //     la::innerProduct<dim>( tmp, &(forceDensityLin(0)), dDW_dl );
                //     nodalForce( dim + d ) = s * tmp * detA * weight;
                // }

                double a0DU_dl[dim];  la::zero<dim>( a0DU_dl );
                a0DU_dl[ d ] = dPhiDXi( 0, l );
                double a1DU_dl[dim];  la::zero<dim>( a1DU_dl );
                a1DU_dl[ d ] = dPhiDXi( 1, l );

                for ( unsigned k=0; k<numFunctions; ++k ) {
                    for ( unsigned c=0; c<dim; ++c ) {

                        // 1st derivative of tangents
                        MatLDimNF tangCoeff_k;
                        this->tangentsGradCur_( k, tangCoeff_k );

                        // 1st derivatives of vectors
                        double pDU_ck[dim];           double pLGradDU_ck[localDim][dim];
                        double qDU_ck[dim];           double qLGradDU_ck[localDim][dim];
                        double qDW_ck[dim];           double qLGradDW_ck[localDim][dim];
                        double dDU_ck[dim];           double cDU_ck[localDim][dim];
                        double dDW_ck[dim];           double cDW_ck[localDim][dim];
                        this->vectorsGradCur_( c, k, xi, phi, dPhiDXi, ddPhiDDXi, tangCoeff_k,
                                               aRef, aCur, bCur, wCur, zCur,
                                               pCur, pLGrad, qCur, qLGrad, dCur, cCur,
                                               pDU_ck, pLGradDU_ck,
                                               qDU_ck, qDW_ck, qLGradDU_ck, qLGradDW_ck,
                                               dDU_ck, dDW_ck, cDU_ck, cDW_ck );

                        double a0DU_ck[dim];  la::zero<dim>( a0DU_ck );
                        a0DU_ck[ c ] = dPhiDXi( 0, k );
                        double a1DU_ck[dim];  la::zero<dim>( a1DU_ck );
                        a1DU_ck[ c ] = dPhiDXi( 1, k );
                        
                        double pDUDU_dlck[dim];
                        SurfaceGeometer::crossprodVectorHessDof( a0Cur, a0DU_dl, a0DU_ck, zeroDim,
                                                                 a1Cur, a1DU_dl, a1DU_ck, zeroDim,
                                                                 pDUDU_dlck );

                        double qDUDU_dlck[dim];  la::assign<dim>( qDUDU_dlck, pDUDU_dlck );
                        double qDUDW_dlck[dim];  la::zero<dim>( qDUDW_dlck );
                        if ( c < localDim )  qDUDW_dlck[ d ] = phi( k ) * tangCoeff_l( c, k );
                        double qDWDU_dlck[dim];  la::zero<dim>( qDWDU_dlck );
                        if ( d < localDim )  qDWDU_dlck[ c ] = phi( l ) * tangCoeff_k( d, l );
                        double qDWDW_dlck[dim];  la::zero<dim>( qDWDW_dlck );
                        
                        double dDUDU_dlck[dim], dDUDW_dlck[dim], dDWDU_dlck[dim], dDWDW_dlck[dim];
                        SurfaceGeometer::unitVectorHessDof( qLen, qDU_dl, qDU_ck, qDUDU_dlck,
                                                            dCur, dDU_dl, dDU_ck, dDUDU_dlck );
                        SurfaceGeometer::unitVectorHessDof( qLen, qDU_dl, qDW_ck, qDUDW_dlck,
                                                            dCur, dDU_dl, dDW_ck, dDUDW_dlck );
                        SurfaceGeometer::unitVectorHessDof( qLen, qDW_dl, qDU_ck, qDWDU_dlck,
                                                            dCur, dDW_dl, dDU_ck, dDWDU_dlck );
                        SurfaceGeometer::unitVectorHessDof( qLen, qDW_dl, qDW_ck, qDWDW_dlck,
                                                            dCur, dDW_dl, dDW_ck, dDWDW_dlck );

                        double kUU_dlck;
                        la::innerProduct<dim>( kUU_dlck, &(forceDensityLin(0)), dDUDU_dlck );
                        kUU_dlck *= s;
                        double kUW_dlck;
                        la::innerProduct<dim>( kUW_dlck, &(forceDensityLin(0)), dDUDW_dlck );
                        kUW_dlck *= s;
                        double kWU_dlck;
                        la::innerProduct<dim>( kWU_dlck, &(forceDensityLin(0)), dDWDU_dlck );
                        kWU_dlck *= s;
                        double kWW_dlck;
                        la::innerProduct<dim>( kWW_dlck, &(forceDensityLin(0)), dDWDW_dlck );
                        kWW_dlck *= s;

                        elemStiff( l*dof    +d, k*dof    +c ) -= kUU_dlck * detA * weight;  // positive LHS
                        if ( c < localDim )
                            elemStiff( l*dof    +d, k*dof+dim+c ) -= kUW_dlck * detA * weight;
                        if ( d < localDim )
                            elemStiff( l*dof+dim+d, k*dof    +c ) -= kWU_dlck * detA * weight;
                        if ( ( d < localDim ) and ( c < localDim ) )
                            elemStiff( l*dof+dim+d, k*dof+dim+c ) -= kWW_dlck * detA * weight;
                    }
                }
            }
        }
        
        // // store nodal force
        // this->BasisElement::giveSupportNodePtr( l )->addToForce( nodalForce );  // negative RHS

    }

    return;
}

//------------------------------------------------------------------------------
// Evaluate internal force kernel at xi (store results in nodes)
template< typename BELEMENT, typename MAT >
void gshell::fem::ElementStatic<BELEMENT,MAT>::internalForceIntegrand(
    const VecLDim & xi,
    const double & weight
    ) 
{
    // Collect shape functions and their derivatives
    VecNF phi;              //this->BasisElement::sfun( xi, phi );
    MatLDimNF dPhiDXi;      //this->BasisElement::sfunGrad( xi, dPhiDXi );
    MatSDimNF ddPhiDDXi;    //this->BasisElement::sfunHess( xi, ddPhiDDXi );
    this->BasisElement::sfunGradHess( xi, phi, dPhiDXi, ddPhiDDXi );

    // get reference mid-vectors
    double aRef[dim][dim];        double bRef[dim][localDim][dim];
    double wRef[dim];             double zRef[localDim][dim];
    double pRef[dim];             double pRefLGrad[localDim][dim];
    double qRef[dim];             double qRefLGrad[localDim][dim];
    double dRef[dim];             double cRef[localDim][dim];
    this->vectors_( gshell::fem::REFERENCE, xi, phi, dPhiDXi, ddPhiDDXi,
                    aRef, bRef, wRef, zRef,
                    pRef, pRefLGrad, qRef, qRefLGrad, dRef, cRef );

    // get current mid-vectors
    double aCur[dim][dim];        double bCur[dim][localDim][dim];
    double wCur[dim];             double zCur[localDim][dim];
    double pCur[dim];             double pLGrad[localDim][dim];
    double qCur[dim];             double qLGrad[localDim][dim];
    double dCur[dim];             double cCur[localDim][dim];
    this->vectors_( gshell::fem::CURRENT, xi, phi, dPhiDXi, ddPhiDDXi,
                    aCur, bCur, wCur, zCur,
                    pCur, pLGrad, qCur, qLGrad, dCur, cCur );

    // Green-Lagrange strain resultants
    double alphaCur[dim][dim], betaCur[dim][dim];
    this->strainRes_( aRef, dRef, cRef, aCur, dCur, cCur,
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

            // 1st derivative of tangents
            MatLDimNF tangCoeff_l;
            this->tangentsGradCur_( l, tangCoeff_l );

            // 1st derivatives of vectors
            double pDU_dl[dim];                double pLGradDU_dl[localDim][dim];
            double qDU_dl[dim];                double qLGradDU_dl[localDim][dim];
            double qDW_dl[dim];                double qLGradDW_dl[localDim][dim];
            double dDU_dl[dim];                double cDU_dl[localDim][dim];
            double dDW_dl[dim];                double cDW_dl[localDim][dim];
            this->vectorsGradCur_( d, l, xi, phi, dPhiDXi, ddPhiDDXi, tangCoeff_l,
                                   aRef, aCur, bCur, wCur, zCur, 
                                   pCur, pLGrad, qCur, qLGrad, dCur, cCur,
                                   pDU_dl, pLGradDU_dl, qDU_dl, qDW_dl, qLGradDU_dl, qLGradDW_dl,
                                   dDU_dl, dDW_dl, cDU_dl, cDW_dl );

            // 1st derivatives of strains w.r.t. nodal DOFs
            double alphaDU_dl[dim][dim], alphaDW_dl[dim][dim], betaDU_dl[dim][dim], betaDW_dl[dim][dim];
            this->strainResGradCur_( d, l, 
                                     aCur, dCur, dDU_dl, dDW_dl, cCur, cDU_dl, cDW_dl, dPhiDXi,
                                     alphaDU_dl, alphaDW_dl, betaDU_dl, betaDW_dl );

            // nodal force
            double tmp;
            la::innerProduct<dim,dim>( tmp, nCur, alphaDU_dl );
            la::addInnerProduct<dim,dim>( tmp, mCur, betaDU_dl );
            nodalForce( d ) = tmp * detA * weight;
            if ( d < localDim ) {
                la::innerProduct<dim,dim>( tmp, nCur, alphaDW_dl );
                la::addInnerProduct<dim,dim>( tmp, mCur,  betaDW_dl );
                nodalForce( dim + d ) = tmp * detA * weight;
            }

        }
        nodalForce *= -1.0;
        this->BasisElement::giveSupportNodePtr( l )->addToForce( nodalForce );  // negative RHS
    }

    return;
}



//------------------------------------------------------------------------------
// Evaluate stiffness kernel at xi, weigh result and store it on #elemStiff
template< typename BELEMENT, typename MAT >
void gshell::fem::ElementStatic<BELEMENT,MAT>::stiffnessIntegrand(
    const VecLDim & xi,
    const double & weight,
    Eigen::MatrixXd & elemStiff
    ) const
{
    // number of shape functions
    const unsigned numFunctions = this->numFunctions( );
    const unsigned matSize = dof * numFunctions;
    assert( elemStiff.rows() == matSize and
            elemStiff.cols() == matSize );

    // Collect shape functions and their derivatives
    VecNF phi;              //this->BasisElement::sfun( xi, phi );
    MatLDimNF dPhiDXi;      //this->BasisElement::sfunGrad( xi, dPhiDXi );
    MatSDimNF ddPhiDDXi;    //this->BasisElement::sfunHess( xi, ddPhiDDXi );
    this->BasisElement::sfunGradHess( xi, phi, dPhiDXi, ddPhiDDXi );

    // get reference mid-vectors
    double aRef[dim][dim];        double bRef[dim][localDim][dim];
    double wRef[dim];             double zRef[localDim][dim];
    double pRef[dim];             double pRefLGrad[localDim][dim];
    double qRef[dim];             double qRefLGrad[localDim][dim];
    double dRef[dim];             double cRef[localDim][dim];
    this->vectors_( gshell::fem::REFERENCE, xi, phi, dPhiDXi, ddPhiDDXi,
                    aRef, bRef, wRef, zRef,
                    pRef, pRefLGrad, qRef, qRefLGrad, dRef, cRef );

    // get current mid-vectors
    double aCur[dim][dim];        double bCur[dim][localDim][dim];
    double wCur[dim];             double zCur[localDim][dim];
    double pCur[dim];             double pLGrad[localDim][dim];
    double qCur[dim];             double qLGrad[localDim][dim];
    double dCur[dim];             double cCur[localDim][dim];
    this->vectors_( gshell::fem::CURRENT, xi, phi, dPhiDXi, ddPhiDDXi,
                    aCur, bCur, wCur, zCur,
                    pCur, pLGrad, qCur, qLGrad, dCur, cCur );

    const double detA = la::determinant<dim>( aRef );
    const double (&a0Cur)[dim] = aCur[ 0 ];
    const double (&a1Cur)[dim] = aCur[ 1 ];
    //const double (&a2Cur)[dim] = aCur[ 2 ];
    double qLen;  la::norm2<dim>( qLen, qCur );
    double sCur[localDim][dim];
    la::assign<localDim,dim>( sCur, qLGrad );
    la::divide<localDim,dim>( sCur, qLen );
    const double (&b0Cur)[localDim][dim] = bCur[ 0 ]; 
    const double (&b1Cur)[localDim][dim] = bCur[ 1 ];

    // Green-Lagrange strain resultants
    double alphaCur[dim][dim], betaCur[dim][dim];
    this->strainRes_( aRef, dRef, cRef, aCur, dCur, cCur,
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

    //--------------------------------------------------------------------------
    // pre-compute mid-plane vectors and their various derivatives
    double (*a0DU)[dim]               = new double[dim*numFunctions][dim];
    double (*a1DU)[dim]               = new double[dim*numFunctions][dim]; 
    double (*b0DU)[localDim][dim]     = new double[dim*numFunctions][localDim][dim]; 
    double (*b1DU)[localDim][dim]     = new double[dim*numFunctions][localDim][dim];
    double (*pDU)[dim]                = new double[dim*numFunctions][dim];  
    double (*pLGradDU)[localDim][dim] = new double[dim*numFunctions][localDim][dim];
    double (*qDU)[dim]                = new double[dim*numFunctions][dim];  
    double (*qDW)[dim]                = new double[dim*numFunctions][dim];  
    double (*qLGradDU)[localDim][dim] = new double[dim*numFunctions][localDim][dim];
    double (*qLGradDW)[localDim][dim] = new double[dim*numFunctions][localDim][dim];
    double (*dDU)[dim]                = new double[dim*numFunctions][dim];  
    double (*dDW)[dim]                = new double[dim*numFunctions][dim];  
    double (*cDU)[localDim][dim]      = new double[dim*numFunctions][localDim][dim];  
    double (*cDW)[localDim][dim]      = new double[dim*numFunctions][localDim][dim];  
    double (*alphaDU)[dim][dim]       = new double[dim*numFunctions][dim][dim]; 
    double (*alphaDW)[dim][dim]       = new double[dim*numFunctions][dim][dim]; 
    double (*betaDU)[dim][dim]        = new double[dim*numFunctions][dim][dim];
    double (*betaDW)[dim][dim]        = new double[dim*numFunctions][dim][dim];

    for ( unsigned l=0; l<numFunctions; ++l ) {
        for ( unsigned d=0; d<dim; ++d ) {
            const unsigned dl = l*dim + d;

            // 1st derivatives of vectors
            la::zero<dim>( a0DU[ dl ] );  a0DU[ dl ][ d ] = dPhiDXi( 0, l );
            la::zero<dim>( a1DU[ dl ] );  a1DU[ dl ][ d ] = dPhiDXi( 1, l );

            la::zero<localDim,dim>( b0DU[ dl ] );
            la::zero<localDim,dim>( b1DU[ dl ] );
            for ( unsigned be=0; be<localDim; ++be ) {
                b0DU[ dl ][ be ][ d ] = ddPhiDDXi( voigtForward[ 0 ][ be ], l );
                b1DU[ dl ][ be ][ d ] = ddPhiDDXi( voigtForward[ 1 ][ be ], l );
            }

            // 1st derivative of tangents
            MatLDimNF tangCoeff_l;
            this->tangentsGradCur_( l, tangCoeff_l );

            // 1st derivatives of vectors
            this->vectorsGradCur_( d, l, xi, phi, dPhiDXi, ddPhiDDXi, tangCoeff_l,
                                   aRef, aCur, bCur, wCur, zCur, 
                                   pCur, pLGrad, qCur, qLGrad, dCur, cCur,
                                   pDU[ dl ], pLGradDU[ dl ],
                                   qDU[ dl ], qDW[ dl ],
                                   qLGradDU[ dl ], qLGradDW[ dl ],
                                   dDU[ dl ], dDW[ dl ],
                                   cDU[ dl ], cDW[ dl ] );

            // 1st derivatives of strains w.r.t. nodal DOFs
            this->strainResGradCur_( d, l,
                                     aCur, dCur, dDU[ dl ], dDW[ dl ],
                                     cCur, cDU[ dl ], cDW[ dl ], dPhiDXi,
                                     alphaDU[ dl ], alphaDW[ dl ],
                                     betaDU[ dl ], betaDW[ dl ] );

        }
    }

    //--------------------------------------------------------------------------
    // stiffness matrix
    for ( unsigned l=0; l<numFunctions; ++l ) {
        for ( unsigned d=0; d<dim; ++d ) {
            const unsigned dl = l*dim + d;

            // 1st derivative of tangents
            MatLDimNF tangCoeff_l;
            this->tangentsGradCur_( l, tangCoeff_l );

            // base vectors
            const double (&a0DU_dl)[dim]               = a0DU[ dl ];
            const double (&a1DU_dl)[dim]               = a1DU[ dl ];
            const double (&b0DU_dl)[localDim][dim]     = b0DU[ dl ];
            const double (&b1DU_dl)[localDim][dim]     = b1DU[ dl ];

            // 1st derivatives of vectors
            //const double (&pDU_dl)[dim]      = pDU[ dl ];
            //const double (&pLGradDU_dl)[localDim][dim] = pLGradDU[ dl ];
            const double (&qDU_dl)[dim]                = qDU[ dl ];
            const double (&qDW_dl)[dim]                = qDW[ dl ];
            const double (&qLGradDU_dl)[localDim][dim] = qLGradDU[ dl ];
            const double (&qLGradDW_dl)[localDim][dim] = qLGradDW[ dl ];
            const double (&dDU_dl)[dim]                = dDU[ dl ];
            const double (&dDW_dl)[dim]                = dDW[ dl ];
            const double (&cDU_dl)[localDim][dim]      = cDU[ dl ];
            const double (&cDW_dl)[localDim][dim]      = cDW[ dl ];

            // 1st derivatives of strains w.r.t. nodal DOFs
            const double (&alphaDU_dl)[dim][dim]       = alphaDU[ dl ];
            const double (&alphaDW_dl)[dim][dim]       = alphaDW[ dl ];
            const double (&betaDU_dl)[dim][dim]        = betaDU[ dl ];
            const double (&betaDW_dl)[dim][dim]        = betaDW[ dl ];

            //------------------------------------------------------------------
            // "elastic and initial-displacement" constribution to stiffness matrix
            for ( unsigned k=0; k<numFunctions; ++k ) {
                for ( unsigned c=0; c<dim; ++c ) {
                    const unsigned ck = k*dim + c;

                    // 1st derivatives of vectors
                    //const double (&pDU_ck)[dim]      = pDU( c, k );
                    //const double (&pLGradDU_ck)[localDim][dim] = pLGradDU( c, k );
                    const double (&qDU_ck)[dim]                = qDU[ ck ];
                    const double (&qDW_ck)[dim]                = qDW[ ck ];
                    const double (&qLGradDU_ck)[localDim][dim] = qLGradDU[ ck ];
                    const double (&qLGradDW_ck)[localDim][dim] = qLGradDW[ ck ];
                    const double (&dDU_ck)[dim]                = dDU[ ck ];
                    const double (&dDW_ck)[dim]                = dDW[ ck ];
                    const double (&cDU_ck)[localDim][dim]      = cDU[ ck ];
                    const double (&cDW_ck)[localDim][dim]      = cDW[ ck ];

                    // 1st derivatives of strains w.r.t. nodal DOFs
                    const double (&alphaDU_ck)[dim][dim]       = alphaDU[ ck ];
                    const double (&alphaDW_ck)[dim][dim]       = alphaDW[ ck ];
                    const double (&betaDU_ck)[dim][dim]        = betaDU[ ck ];
                    const double (&betaDW_ck)[dim][dim]        = betaDW[ ck ];

                    // "elastic and initial-displacement" stiffness components
                    double kUU_dlck = 0.0;
                    double kUW_dlck = 0.0;
                    double kWU_dlck = 0.0;
                    double kWW_dlck = 0.0;
                    for ( unsigned i=0; i<dim; ++i ) {
                        for ( unsigned j=0; j<dim; ++j ) {
                            double tmp;
                            la::innerProduct<dim,dim>( tmp, hRef[ i ][ j ], alphaDU_ck );
                            tmp *= alphaDU_dl[ i ][ j ] * t;
                            kUU_dlck += tmp;
                            la::innerProduct<dim,dim>( tmp, hRef[ i ][ j ], alphaDW_ck );
                            tmp *= alphaDU_dl[ i ][ j ] * t;
                            kUW_dlck += tmp;
                            la::innerProduct<dim,dim>( tmp, hRef[ i ][ j ], alphaDU_ck );
                            tmp *= alphaDW_dl[ i ][ j ] * t;
                            kWU_dlck += tmp;
                            la::innerProduct<dim,dim>( tmp, hRef[ i ][ j ], alphaDW_ck );
                            tmp *= alphaDW_dl[ i ][ j ] * t;
                            kWW_dlck += tmp;
                            la::innerProduct<dim,dim>( tmp, hRef[ i ][ j ],  betaDU_ck );
                            tmp *= betaDU_dl[ i ][ j ] * s;
                            kUU_dlck += tmp;
                            la::innerProduct<dim,dim>( tmp, hRef[ i ][ j ],  betaDW_ck );
                            tmp *= betaDU_dl[ i ][ j ] * s;
                            kUW_dlck += tmp;
                            la::innerProduct<dim,dim>( tmp, hRef[ i ][ j ],  betaDU_ck );
                            tmp *= betaDW_dl[ i ][ j ] * s;
                            kWU_dlck += tmp;
                            la::innerProduct<dim,dim>( tmp, hRef[ i ][ j ],  betaDW_ck );
                            tmp *= betaDW_dl[ i ][ j ] * s;
                            kWW_dlck += tmp;
                        }
                    }

                    elemStiff( l*dof    +d, k*dof    +c ) += kUU_dlck * detA * weight;  // positive LHS
                    if ( c < localDim )
                        elemStiff( l*dof    +d, k*dof+dim+c ) += kUW_dlck * detA * weight;
                    if ( d < localDim )
                        elemStiff( l*dof+dim+d, k*dof    +c ) += kWU_dlck * detA * weight;
                    if ( ( d < localDim ) and ( c < localDim ) )
                        elemStiff( l*dof+dim+d, k*dof+dim+c ) += kWW_dlck * detA * weight;

                    //----------------------------------------------------------
                    // "geometric" contribution to stiffness matrix
//#define GSHELL_SKIP_GEOMETRICSTIFFNESS
#ifndef GSHELL_SKIP_GEOMETRICSTIFFNESS

                    // base vectors
                    const double (&a0DU_ck)[dim] = a0DU[ ck ];
                    const double (&a1DU_ck)[dim] = a1DU[ ck ];
                    double pDUDU_dlck[dim];
                    SurfaceGeometer::crossprodVectorHessDof( a0Cur, a0DU_dl, a0DU_ck, zeroDim,
                                                             a1Cur, a1DU_dl, a1DU_ck, zeroDim,
                                                             pDUDU_dlck );

                    const double (&b0DU_ck)[localDim][dim] = b0DU[ ck ];
                    const double (&b1DU_ck)[localDim][dim] = b1DU[ ck ];                    
                    double pLGradDUDU_dlck[localDim][dim];
                    SurfaceGeometer::crossprodVectorDeriv1HessDof( a0Cur, a0DU_dl, a0DU_ck, zeroDim,
                                                                   b0Cur, b0DU_dl, b0DU_ck, zeroLDimDim,
                                                                   a1Cur, a1DU_dl, a1DU_ck, zeroDim,
                                                                   b1Cur, b1DU_dl, b1DU_ck, zeroLDimDim,
                                                                   pLGradDUDU_dlck );

                    // 1st derivative of tangents
                    MatLDimNF tangCoeff_k;
                    this->tangentsGradCur_( k, tangCoeff_k );
                    
                    double qDUDU_dlck[dim];  la::assign<dim>( qDUDU_dlck, pDUDU_dlck );
                    double qDUDW_dlck[dim];  la::zero<dim>( qDUDW_dlck );
                    if ( c < localDim ) // gamma=c
                        qDUDW_dlck[ d ] = phi( k ) * tangCoeff_l( c, k );
                    double qDWDU_dlck[dim];  la::zero<dim>( qDWDU_dlck );
                    if ( d < localDim )  // delta=d
                        qDWDU_dlck[ c ] = phi( l ) * tangCoeff_k( d, l );
                    double qDWDW_dlck[dim];  la::zero<dim>( qDWDW_dlck );
                    
                    double dDUDU_dlck[dim], dDUDW_dlck[dim], dDWDU_dlck[dim], dDWDW_dlck[dim];
                    SurfaceGeometer::unitVectorHessDof( qLen, qDU_dl, qDU_ck, qDUDU_dlck,
                                                        dCur, dDU_dl, dDU_ck, dDUDU_dlck );
                    SurfaceGeometer::unitVectorHessDof( qLen, qDU_dl, qDW_ck, qDUDW_dlck,
                                                        dCur, dDU_dl, dDW_ck, dDUDW_dlck );
                    SurfaceGeometer::unitVectorHessDof( qLen, qDW_dl, qDU_ck, qDWDU_dlck,
                                                        dCur, dDW_dl, dDU_ck, dDWDU_dlck );
                    SurfaceGeometer::unitVectorHessDof( qLen, qDW_dl, qDW_ck, qDWDW_dlck,
                                                        dCur, dDW_dl, dDW_ck, dDWDW_dlck );

                    double qLGradDUDU_dlck[localDim][dim];  la::assign<localDim,dim>( qLGradDUDU_dlck, pLGradDUDU_dlck );
                    double qLGradDUDW_dlck[localDim][dim];  la::zero<localDim,dim>( qLGradDUDW_dlck );
                    if ( c < localDim )
                        for ( unsigned be=0; be<localDim; ++be )
                            qLGradDUDW_dlck[ be ][ d ] = dPhiDXi( be, k ) * tangCoeff_l( c, k );
                    double qLGradDWDU_dlck[localDim][dim];  la::zero<localDim,dim>( qLGradDWDU_dlck );
                    if ( d < localDim )
                        for ( unsigned be=0; be<localDim; ++be )
                            qLGradDWDU_dlck[ be ][ c ] = dPhiDXi( be, l ) * tangCoeff_k( d, l );
                    double qLGradDWDW_dlck[localDim][dim];  la::zero<localDim,dim>( qLGradDWDW_dlck );

                    double cDUDU_dlck[localDim][dim], cDUDW_dlck[localDim][dim];
                    double cDWDU_dlck[localDim][dim], cDWDW_dlck[localDim][dim];
                    SurfaceGeometer::unitVectorDeriv1HessDof( qLen, qDU_dl, qDU_ck, qDUDU_dlck,
                                                              qLGradDU_dl, qLGradDU_ck, qLGradDUDU_dlck,
                                                              dCur, dDU_dl, dDU_ck, dDUDU_dlck,
                                                              sCur, cDUDU_dlck );
                    SurfaceGeometer::unitVectorDeriv1HessDof( qLen, qDU_dl, qDW_ck, qDUDW_dlck,
                                                              qLGradDU_dl, qLGradDW_ck, qLGradDUDW_dlck,
                                                              dCur, dDU_dl, dDW_ck, dDUDW_dlck,
                                                              sCur, cDUDW_dlck );
                    SurfaceGeometer::unitVectorDeriv1HessDof( qLen, qDW_dl, qDU_ck, qDWDU_dlck,
                                                              qLGradDW_dl, qLGradDU_ck, qLGradDWDU_dlck,
                                                              dCur, dDW_dl, dDU_ck, dDWDU_dlck,
                                                              sCur, cDWDU_dlck );
                    SurfaceGeometer::unitVectorDeriv1HessDof( qLen, qDW_dl, qDW_ck, qDWDW_dlck,
                                                              qLGradDW_dl, qLGradDW_ck, qLGradDWDW_dlck,
                                                              dCur, dDW_dl, dDW_ck, dDWDW_dlck,
                                                              sCur, cDWDW_dlck );

                    // "geometric" stiffness components
                    double kgUU_dlck = 0.0;
                    double kgUW_dlck = 0.0;
                    double kgWU_dlck = 0.0;
                    double kgWW_dlck = 0.0;
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

                            tmp =
                                dPhiDXi( al, l ) * cDW_ck[ be ][ d ] +
                                dPhiDXi( be, l ) * cDW_ck[ al ][ d ];
                            la::addInnerProduct<dim>( tmp, aCur[ al ], cDUDW_dlck[ be ] );
                            la::addInnerProduct<dim>( tmp, aCur[ be ], cDUDW_dlck[ al ] );
                            kgUW_dlck += 0.5 * mCur[ al ][ be ] * tmp;

                            tmp =
                                cDW_dl[ be ][ c ] * dPhiDXi( al, k ) +
                                cDW_dl[ al ][ c ] * dPhiDXi( be, k );
                            la::addInnerProduct<dim>( tmp, aCur[ al ], cDWDU_dlck[ be ] );
                            la::addInnerProduct<dim>( tmp, aCur[ be ], cDWDU_dlck[ al ] );
                            kgWU_dlck += 0.5 * mCur[ al ][ be ] * tmp;

                            tmp = 0.;
                            la::addInnerProduct<dim>( tmp, aCur[ al ], cDWDW_dlck[ be ] );
                            la::addInnerProduct<dim>( tmp, aCur[ be ], cDWDW_dlck[ al ] );
                            kgWW_dlck += 0.5 * mCur[ al ][ be ] * tmp;
                        }
                    }
                    // al,3=3,al
                    for ( unsigned al=0; al<localDim; ++al ) {
                        double tmp =
                            dPhiDXi( al, l ) * dDU_ck[ d ] +
                            dDU_dl[ c ] * dPhiDXi( al, k );
                        la::addInnerProduct<dim>( tmp, aCur[ al ], dDUDU_dlck );
                        kgUU_dlck += 2.0 * 0.5 * nCur[ al ][ 2 ] * tmp;

                        tmp = dPhiDXi( al, l ) * dDW_ck[ d ];
                        la::addInnerProduct<dim>( tmp, aCur[ al ], dDUDW_dlck );
                        kgUW_dlck += 2.0 * 0.5 * nCur[ al ][ 2 ] * tmp;

                        tmp = dDW_dl[ c ] * dPhiDXi( al, k );
                        la::addInnerProduct<dim>( tmp, aCur[ al ], dDWDU_dlck );
                        kgWU_dlck += 2.0 * 0.5 * nCur[ al ][ 2 ] * tmp;
                        
                        tmp = 0.;
                        la::addInnerProduct<dim>( tmp, aCur[ al ], dDWDW_dlck );
                        kgWW_dlck += 2.0 * 0.5 * nCur[ al ][ 2 ] * tmp;
                    }
                    // 3,3
                    ;  // do nothing

                    elemStiff( l*dof    +d, k*dof    +c ) += kgUU_dlck * detA * weight;  // positive LHS
                    if ( c < localDim )
                        elemStiff( l*dof    +d, k*dof+dim+c ) += kgUW_dlck * detA * weight;
                    if ( d < localDim )
                        elemStiff( l*dof+dim+d, k*dof    +c ) += kgWU_dlck * detA * weight;
                    if ( ( d < localDim ) and ( c < localDim ) )
                        elemStiff( l*dof+dim+d, k*dof+dim+c ) += kgWW_dlck * detA * weight;
#endif

                }
            }

        }
    }

    // deallocate pre-computed mid-plane vectors and their various derivatives
    delete [] a0DU;
    delete [] a1DU;
    delete [] b0DU;
    delete [] b1DU;
    delete [] pDU;
    delete [] pLGradDU;
    delete [] qDU;
    delete [] qDW;
    delete [] qLGradDU;
    delete [] qLGradDW;
    delete [] dDU;
    delete [] dDW;
    delete [] cDU;
    delete [] cDW;
    delete [] alphaDU;
    delete [] alphaDW;
    delete [] betaDU;
    delete [] betaDW;

    // done
    return;
}

//------------------------------------------------------------------------------
// Compute strain resultants -- for output
template< typename BELEMENT, typename MAT >
corlib::eigenX::MatrixSd<3, 3>
gshell::fem::ElementStatic<BELEMENT,MAT>::giveStrainResultant(
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
    double aRef[dim][dim];        double bRef[dim][localDim][dim];
    double wRef[dim];             double zRef[localDim][dim];
    double pRef[dim];             double pRefLGrad[localDim][dim];
    double qRef[dim];             double qRefLGrad[localDim][dim];
    double dRef[dim];             double cRef[localDim][dim];
    this->vectors_( gshell::fem::REFERENCE, xi, phi, dPhiDXi, ddPhiDDXi,
                    aRef, bRef, wRef, zRef,
                    pRef, pRefLGrad, qRef, qRefLGrad, dRef, cRef );

    // get current mid-vectors
    double aCur[dim][dim];        double bCur[dim][localDim][dim];
    double wCur[dim];             double zCur[localDim][dim];
    double pCur[dim];             double pLGrad[localDim][dim];
    double qCur[dim];             double qLGrad[localDim][dim];
    double dCur[dim];             double cCur[localDim][dim];
    this->vectors_( gshell::fem::CURRENT, xi, phi, dPhiDXi, ddPhiDDXi,
                    aCur, bCur, wCur, zCur,
                    pCur, pLGrad, qCur, qLGrad, dCur, cCur );

    // Green-Lagrange strain resultants
    double alphaCur[dim][dim], betaCur[dim][dim];
    this->strainRes_( aRef, dRef, cRef, aCur, dCur, cCur,
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
// Compute stress resultants -- for output
template< typename BELEMENT, typename MAT >
corlib::eigenX::MatrixSd<3, 3>
gshell::fem::ElementStatic<BELEMENT,MAT>::giveStressResultant(
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
    double aRef[dim][dim];        double bRef[dim][localDim][dim];
    double wRef[dim];             double zRef[localDim][dim];
    double pRef[dim];             double pRefLGrad[localDim][dim];
    double qRef[dim];             double qRefLGrad[localDim][dim];
    double dRef[dim];             double cRef[localDim][dim];
    this->vectors_( gshell::fem::REFERENCE, xi, phi, dPhiDXi, ddPhiDDXi,
                    aRef, bRef, wRef, zRef,
                    pRef, pRefLGrad, qRef, qRefLGrad, dRef, cRef );

    // get mid-vectors
    double aCur[dim][dim];        double bCur[dim][localDim][dim];
    double wCur[dim];             double zCur[localDim][dim];
    double pCur[dim];             double pLGrad[localDim][dim];
    double qCur[dim];             double qLGrad[localDim][dim];
    double dCur[dim];             double cCur[localDim][dim];
    this->vectors_( gshell::fem::CURRENT, xi, phi, dPhiDXi, ddPhiDDXi,
                    aCur, bCur, wCur, zCur,
                    pCur, pLGrad, qCur, qLGrad, dCur, cCur );

    // Green-Lagrange strain resultants
    double alphaCur[dim][dim], betaCur[dim][dim];
    this->strainRes_( aRef, dRef, cRef, aCur, dCur, cCur,
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

