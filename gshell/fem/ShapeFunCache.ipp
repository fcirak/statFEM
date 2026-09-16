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
template< corlib::shape SHAPE >
const unsigned gshell::fem::ShapeFunCache<SHAPE>::
voigtForward[localDim][localDim] = { { 0, 1 }, { 1, 2 } };

//------------------------------------------------------------------------------
// Populate the shape functions and their (parametric) derivatives at given 
// quadrature rule
template< corlib::shape SHAPE >
template< typename SFUN, typename QUADRATURE >
void gshell::fem::ShapeFunCache<SHAPE>::populateEvaluatedShapeFunctions(
    SFUN * shapeFun,
    const QUADRATURE & quadrature
    )
{
    typedef QUADRATURE                       Quadrature;
    typedef typename Quadrature::QuadIter    QuadIter;

    const QuadIter qBegin = quadrature.begin();
    const QuadIter qEnd   = quadrature.end();

    VecNF     phi;
    MatLDimNF dPhiDXi;
    MatSDimNF ddPhiDDXi;

    for ( QuadIter qIter = qBegin; qIter != qEnd; ++qIter ) {
        const VecLDim xi = qIter -> second;
        phi.setZero();  dPhiDXi.setZero();  ddPhiDDXi.setZero();
        shapeFun -> evaluateGradHess( xi, phi, dPhiDXi, ddPhiDDXi );
        evaluatedFun_.insert( std::make_pair( xi, ShapeFunGradHess_( phi, dPhiDXi, ddPhiDDXi ) ) );
    }

    // add center (this going to be needed for output of element-wise data to VTU)
    const VecLDim xiCenter = corlib::ShapeTraits<myShape>::centroid();
    const MapLCoordToShapeFunConstIter_ xiEvalFun = evaluatedFun_.find( xiCenter );
    if ( xiEvalFun == evaluatedFun_.end() ) {  // only add if not already there
        phi.setZero();  dPhiDXi.setZero();  ddPhiDDXi.setZero();
        shapeFun -> evaluateGradHess( xiCenter, phi, dPhiDXi, ddPhiDDXi );
        evaluatedFun_.insert( std::make_pair( xiCenter, ShapeFunGradHess_( phi, dPhiDXi, ddPhiDDXi ) ) );
    }
}

//------------------------------------------------------------------------------
// Populate limit surface coefficients at vertices
template< corlib::shape SHAPE >
template< typename SFUN >
void gshell::fem::ShapeFunCache<SHAPE>::populateLimitCoefficientsAtVertices(
    SFUN * shapeFun
    )
{
    const unsigned numFunctions = this -> numFunctions();
    FTL_VERIFY( numFunctions > 0 );

    for ( unsigned v = 0; v < numVertices; ++v ) {
        VecNF & coeff = limitCoeffs_[ v ];
        //coeff.resize( localDim, numFunctions );
        coeff.setZero();
        shapeFun -> evaluate( v, coeff );
    }
}

//------------------------------------------------------------------------------
// Populate tangent coefficients at vertices
template< corlib::shape SHAPE >
template< typename SFUN >
void gshell::fem::ShapeFunCache<SHAPE>::populateTangentCoefficientsAtVertices(
    SFUN * shapeFun
    )
{
    const unsigned numFunctions = this -> numFunctions();
    FTL_VERIFY( numFunctions > 0 );

    for ( unsigned v = 0; v < numVertices; ++v ) {
        MatLDimNF & coeff = tangCoeffs_[ v ];
        //coeff.resize( localDim, numFunctions );
        coeff.setZero();
        shapeFun -> evaluateGradient( v, coeff );
    }
}


//------------------------------------------------------------------------------
// Evaluates the shape function at the given coordinate xi
template< corlib::shape SHAPE >
void gshell::fem::ShapeFunCache<SHAPE>::evaluate(
    const VecLDim & xi, 
    VecNF & phi
    ) const
{
    const unsigned numFunctions = this -> numFunctions();

    phi.resize( numFunctions );
 
    const MapLCoordToShapeFunConstIter_ xiEvalFun = evaluatedFun_.find( xi );
    FTL_VERIFY_DESCRIPTIVE( xiEvalFun != evaluatedFun_.end(),
                            "Could not find evaluated shape functions"
                            " at given co-ordinate\n" );
    phi = std::get< 0 >( xiEvalFun -> second );

    return;
}

//------------------------------------------------------------------------------
// Compute the gradient w.r.t xi of the shape functions at the given xi
template< corlib::shape SHAPE >
void gshell::fem::ShapeFunCache<SHAPE>::evaluateGradient(
    const VecLDim & xi, 
    MatLDimNF & dPhiDXi
    ) const
{
    const unsigned numFunctions = this -> numFunctions();

    dPhiDXi.resize( localDim, numFunctions );

    const MapLCoordToShapeFunConstIter_ xiEvalFun = evaluatedFun_.find( xi );
    FTL_VERIFY_DESCRIPTIVE( xiEvalFun != evaluatedFun_.end(),
                            "Could not find evaluated 1st derivatives of shape functions"
                            " at given co-ordinate\n" );
    dPhiDXi = std::get< 1 >( xiEvalFun -> second );

    return;
}

//------------------------------------------------------------------------------
// Compute the 2nd derivative w.r.t xi of the shape functions at the given xi:
template< corlib::shape SHAPE >
void gshell::fem::ShapeFunCache<SHAPE>::evaluateHessian(
    const VecLDim & xi, 
    MatSDimNF & ddPhiDDXi
    ) const
{
    const unsigned numFunctions = this -> numFunctions();

    ddPhiDDXi.resize( sDim, numFunctions );

    const MapLCoordToShapeFunConstIter_ xiEvalFun = evaluatedFun_.find( xi );
    FTL_VERIFY_DESCRIPTIVE( xiEvalFun != evaluatedFun_.end(),
                            "Could not find evaluated 2nd derivatives of shape functions"
                            " at given co-ordinate\n" );
    ddPhiDDXi = std::get< 2 >( xiEvalFun -> second );

    return;
}

//------------------------------------------------------------------------------
// Return value, 1st and 2nd derivatives of shape functions at xi
template< corlib::shape SHAPE >
void gshell::fem::ShapeFunCache<SHAPE>::evaluateGradHess(
    const VecLDim & xi,
    VecNF & phi,
    MatLDimNF & dPhiDXi,
    MatSDimNF & ddPhiDDXi
    ) const
{
    const unsigned numFunctions = this -> numFunctions();
    phi.resize( numFunctions );
    dPhiDXi.resize( localDim, numFunctions );
    ddPhiDDXi.resize( sDim, numFunctions );

    const MapLCoordToShapeFunConstIter_ xiEvalFun = evaluatedFun_.find( xi );
    FTL_VERIFY_DESCRIPTIVE( xiEvalFun != evaluatedFun_.end(),
                            "Could not find evaluated shape functions"
                            " (and their 1st and 2nd derivatives)"
                            " at given co-ordinate\n" );
    if ( xiEvalFun != evaluatedFun_.end() ) {
        FTL_VERIFY_DESCRIPTIVE( corlib::fuzzyEqual( (xiEvalFun -> first)[ 0 ], xi[ 0 ]),
                                "%g != %g\n", (xiEvalFun -> first)[ 0 ], xi[ 0 ] );
        FTL_VERIFY_DESCRIPTIVE( corlib::fuzzyEqual( (xiEvalFun -> first)[ 1 ], xi[ 1 ]),
                                "%g != %g\n", (xiEvalFun -> first)[ 1 ], xi[ 1 ] );
        phi       = std::get< 0 >( xiEvalFun -> second );
        dPhiDXi   = std::get< 1 >( xiEvalFun -> second );
        ddPhiDDXi = std::get< 2 >( xiEvalFun -> second );
    }

    return;
}

