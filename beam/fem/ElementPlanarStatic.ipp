// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//
//                 Computational Structural Mechanics Lab
//                         University of Cambridge
//
// This software is copyrighted and all rights are retained by the University of
// Cambridge. All use, disclosure, and/or reproduction of any part not expressly
// authorised by the University of Cambridge is prohibited. (C) 2010
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

//! @file   ElementPlanarStatic.ipp

#include <corlib/linalg.hpp>
#include <functional>

//------------------------------------------------------------------------------
/** Evaluates the shape function at the given local/parametric coordinate
 *  \param[in]  xi   Local cooordinate at which the functions is evaluated:
 *                   \f$ \theta^1 \f$
 *  \param[out] phi  Result of all shape functions:
 *                   \f$ \varphi^j(\theta^1) \f$
 */
template< typename NODE, typename SFUN, typename MAT >
void beam::fem::ElementPlanarStatic<NODE,SFUN,MAT>::sFun_( const VecLDim& xi, 
                                                           VecNF& phi ) const
{
    shapeFun_.evaluate( xi, phi );
    return;
}

//------------------------------------------------------------------------------
/** Compute the parametric gradient of the shape functions at the given point:
 * \param[in]  xi       Local coordinate at which the gradient is evaluated:
 *                      \f$ \theta^1 \f$
 * \param[out] dphiDxi  gradient:
 *                      \f$ (d \varphi^j / d \theta^1)(\theta^1) \f$
 */
template< typename NODE, typename SFUN, typename MAT >
void beam::fem::ElementPlanarStatic<NODE,SFUN,MAT>::sFunGrad_( const VecLDim& xi, 
                                                               MatLDimNF& dphiDxi ) const 
{
    shapeFun_.evaluateGradient( xi, dphiDxi );
    return;
}

//------------------------------------------------------------------------------
/** Compute the 2nd parametric derivative (Hessian) of the shape functions at the given point
 * \param[in]  xi        Local coordinate at which the gradient is evaluated:
 *                       \f$ \theta^1 \f$
 * \param[out] ddphiDDxi Hessian:
 *                       \f$ (d^2 \varphi^j / (d \theta^1)^2 )(\theta^1) \f$
 */
template< typename NODE, typename SFUN, typename MAT >
void beam::fem::ElementPlanarStatic<NODE,SFUN,MAT>::sFunGradGrad_( const VecLDim& xi, 
                                                                   MatVecNFLDimLDim & ddphiDDxi ) const
{
    shapeFun_.evaluateHessian( xi, ddphiDDxi );
    return;
}

//------------------------------------------------------------------------------
/** Compute resultant elasticity resultant, ie 'EA' and 'EI'
 *  \param[in]  baseRef       Basis vectors in reference configuration
 *  \param[out] elastModMemb  membrane elasticity modulus
 *  \param[out] elastModBend  bending elasticity modulus
 */
template< typename NODE, typename SFUN, typename MAT>
void beam::fem::ElementPlanarStatic<NODE,SFUN,MAT>::elasticityRes( const Mat3x3& baseRef,
                                                                   double& elastModMemb,
                                                                   double& elastModBend ) const
{
    // _simplified_ contra-variant metric 11-coefficient \f$G^{11}\f$
    // is the inverse of the squared Euclidian length of the axial
    // base vector \f$A^0\f$
    const double bas0RefLenSq = baseRef.col(0).squaredNorm( );

    // Young's modulus
    const double youngMod = material_->estimatedElasticityConst();

    /// - axial elasticity: \f$ \bar{EA} = \frac{EA}{ |A^0|^4 } \f$
    elastModMemb = youngMod * area_ / ( bas0RefLenSq * bas0RefLenSq );
    /// - bending elasticity: \f$ \bar{EI} = \frac{EI}{ |A^0|^4 } \f$
    elastModBend = youngMod * secMomArea_  / ( bas0RefLenSq * bas0RefLenSq );

    // quit
    return;
}

//------------------------------------------------------------------------------
/** Accessor to the indices of the nodal degrees of freedom
 *  \tparam    OUT  Output iterator
 *  \param[in] iter Iterator which provides access to some dof storage
 */
template< typename NODE, typename SFUN, typename MAT >
void beam::fem::ElementPlanarStatic<NODE,SFUN,MAT>::
getDofIndices( std::vector<unsigned> & dofIndices ) const
{
    dofIndices.resize( matSize );
    std::vector<unsigned>::iterator iter = dofIndices.begin();

    for ( unsigned n=0; n<numNodesSN; ++n ) { 
        std::vector<unsigned> aux;
        nodes_[ n ]->copyDofArray( aux );
        iter = std::copy( aux.begin(), aux.end(), iter );
    }
    return;
}

//------------------------------------------------------------------------------
/** Evaluate a quantity at given global coordinate. The global co-ordinate is
 * translated to local co-ordinate by looking at the line formed by
 * the left and right vertex.
 */
template< typename NODE, typename SFUN, typename MAT >
template< typename OP > 
typename OP::ReturnType 
beam::fem::ElementPlanarStatic<NODE,SFUN,MAT>::
interpolateNodalQuantity( const OP& op, const VecDim& coorCur ) const
{
    // get local co-ordinate of #coorCur with respect 
    typename ElementBasic_::NodeConstIterator nodeIter
        = ElementBasic_::nodesBegin();

    const VecDim leftVertexCur  = 
        (*nodeIter) -> giveCoordinates() +
        (*nodeIter) -> giveDisplacements();
    ++nodeIter;
    const VecDim rightVertexCur = 
        (*nodeIter) -> giveCoordinates() +
        (*nodeIter) -> giveDisplacements();
    const VecDim edge = rightVertexCur - leftVertexCur;
    const VecDim dist = coorCur - leftVertexCur;
    double coorLoc = 
        dist.dot( edge ) / 
        edge.dot( edge );
    
    // snap to [0, 1]
    if ( coorLoc < 0.0 ) coorLoc = 0.0;
    if ( coorLoc > 1.0 ) coorLoc = 1.0;
 
    VecLDim xi; xi[ 0 ] = coorLoc;

    // shape function
    VecNF phi;  this -> sFun_( xi, phi );
        
    typename OP::ReturnType result = OP::zero();
    for ( unsigned n = 0; n < numNodesSN; ++n ) {
        result += phi[ n ] * op( nodes_[ n ] );
    }
    return result;
}

//------------------------------------------------------------------------------
/** Return interpolated reference coordinate at evaluation point
 */
template< typename NODE, typename SFUN, typename MAT >
typename beam::fem::ElementPlanarStatic<NODE,SFUN,MAT>::VecDim
beam::fem::ElementPlanarStatic<NODE,SFUN,MAT>::giveReferenceCoordinate( const VecLDim& xi ) const
{
    // nodal coordinates
    MatDimNF coorRef;  this -> nodalCoordinatesRef( coorRef );
    // shape function derivatives
    VecNF phi;  this->sFun_( xi, phi );

    // interpolate geometry
    return coorRef * phi;
}

//------------------------------------------------------------------------------
/** Return interpolated displacement at evaluation point
 */
template< typename NODE, typename SFUN, typename MAT >
typename beam::fem::ElementPlanarStatic<NODE,SFUN,MAT>::VecDim
beam::fem::ElementPlanarStatic<NODE,SFUN,MAT>::giveDisplacement( const VecLDim& xi ) const
{
    // nodal coordinates
    MatDimNF dispN;  this->nodalDisplacements( dispN );
    // shape function derivatives
    VecNF phi;  this->sFun_( xi, phi );

    // interpolate displacements
    return dispN * phi;
}

//------------------------------------------------------------------------------
/** Return unit normal to current configuration
 */
template< typename NODE, typename SFUN, typename MAT >
typename beam::fem::ElementPlanarStatic<NODE,SFUN,MAT>::VecDim
beam::fem::ElementPlanarStatic<NODE,SFUN,MAT>::giveUnitNormalCur( const VecLDim& xi ) const
{
    Mat3x3 base;
    this -> baseCur_( xi, base );
    return - ( base.col( 1 ) ).head( dim );
}

//------------------------------------------------------------------------------
/** Return scaled normal to current configuration whose length carried an
 *  area-scale to the reference configuration
 */
template< typename NODE, typename SFUN, typename MAT >
typename beam::fem::ElementPlanarStatic<NODE,SFUN,MAT>::VecDim
beam::fem::ElementPlanarStatic<NODE,SFUN,MAT>::giveScaledNormal( const VecLDim& xi ) const
{
    Mat3x3 baseRef, baseCur;
    this -> baseRef_( xi, baseRef );
    this -> baseCur_( xi, baseCur );
    const VecDim unitNormalCur = - ( baseCur.col( 1 ) ).head( dim );
    const double scale = ( baseCur.col( 0 ).norm( ) ) /
        ( area_ * ( baseRef.col( 0 ).norm( ) ) );
    return ( scale * unitNormalCur );
}

//------------------------------------------------------------------------------
/** Collect nodal co-ordinates (reference config) an put them into a matrix
 *  \param[out] coorRef   Matrix of nodal reference co-ordinates
 */
template< typename NODE, typename SFUN, typename MAT >
void beam::fem::ElementPlanarStatic<NODE,SFUN,MAT>::nodalCoordinatesRef( MatDimNF& coorRef ) const
{
    // collect nodal vectors based on operation
    for ( unsigned n = 0; n < numNodesSN; ++n ) {
        coorRef.col( n ) = nodes_[ n ] -> giveCoordinates();
    }

    // done
    return;
}

//------------------------------------------------------------------------------
/** Collect nodal co-ordinates (current config) an put them into a matrix
 *  \param[out] coorRef   Matrix of nodal deformed co-ordinates
 */
template< typename NODE, typename SFUN, typename MAT >
void beam::fem::ElementPlanarStatic<NODE,SFUN,MAT>::nodalCoordinatesCur( MatDimNF& coorCur ) const
{
    // collect nodal vectors based on operation
    for ( unsigned n = 0; n < numNodesSN; ++n ) {
        coorCur.col( n ) = 
            ( nodes_[ n ] -> giveCoordinates() ) + 
            ( nodes_[ n ] -> giveDisplacements() );
    }

    // done
    return;
}

//------------------------------------------------------------------------------
/** Accessor to a general nodal quantitiy accessed by the pass operator
 *  \tparam      OP         Node functor which access the desired nodal quantity
 *  \param[out]  quantity   Matrix (dim x number of vertices) containing the result
 *  \param[in]   op         Node functor for access of nodal quantity
 */
template< typename NODE, typename SFUN, typename MAT >
template< typename OP >
void beam::fem::ElementPlanarStatic<NODE,SFUN,MAT>::nodalQuantity_( MatDofNF& quantity, 
                                                                    OP op ) const 
{
    // collect nodal vectors based on operation
    for ( unsigned n = 0; n < numNodesSN; ++n ) {
        quantity.col( n ) = op( nodes_[ n ] );
    }
    return;
}

//------------------------------------------------------------------------------
/** Collect nodal displacements and put them into a matrix
 *  \param[out] disp  Matrix with the nodal displacements as columns:
 *                    \f$ [ u^i{}_K ] \f$
 */
template< typename NODE, typename SFUN, typename MAT >
void beam::fem::ElementPlanarStatic<NODE,SFUN,MAT>::nodalDisplacements( MatDofNF& disp ) const
{
    this->nodalQuantity_( disp, std::mem_fun( &NODE::giveDisplacements ) );
    return;
}

//------------------------------------------------------------------------------
/** Collect nodal displacements and put them into a matrix
 *
 *  \param[out]  dispIncr  Matrix with the nodal displacements as columns:
 *                         \f$ [ u^i{}_K ] \f$
 */
template< typename NODE, typename SFUN, typename MAT >
void beam::fem::ElementPlanarStatic<NODE,SFUN,MAT>::nodalIncrements( MatDofNF& dispIncr ) const
{
    this->nodalQuantity_( dispIncr, std::mem_fun( &NODE::giveIncrement ) );
    return;
}

//------------------------------------------------------------------------------
/** Co-variant base vectors at centroid in reference configuration evaluated
 *  at given local co-ordinate
 *  \param[in]  xi    local/parametric evaluation co-ordinate
 *                    \f$ \theta^1 \f$
 *  \param[out] base  co-variant reference base vectors stored in (Jacobian) matrix;
 *                    the base vectors are given in global Cartesian components
 *                    (ie same co-ordinate system in which co-ordinates \f$ X^I \f$
 *                    of NODEs are given) and stored as column vectors in
 *                    resulting matrix: \f$ [ (A^I)_\alpha ] \f$
 */
template< typename NODE, typename SFUN, typename MAT >
void beam::fem::ElementPlanarStatic<NODE,SFUN,MAT>::baseRef_( const VecLDim& xi,
                                                              Mat3x3& base ) const
{
    // nodal coordinates
    MatDimNF coorRef;  this->nodalCoordinatesRef( coorRef );
    // shape function derivatives
    MatLDimNF dphiDxi; this->sFunGrad_( xi, dphiDxi );

    // blank
    base.setZero();

    /// - axial director, tangent to centerline: \f$ A^0 = d X / d \theta \f$
    base.block( 0, 0, dim, 1 ) = coorRef * dphiDxi.transpose( );

    /// - bi-normal, normal to plane: \f$ A^2 = (0,0,1) \f$
    base( 2, 2 ) = 1.0;

    /// - unit normal to centerline: \f$ A^1 = (A^2 \times A^0) / |A^2 \times A^0| \f$
    base.col( 1 ) = 
        ( base.col( 2 ).cross( base.col( 0 ) ) ) / 
        ( base.col( 0 ).norm( ) );

    // quit
    return;
}

//------------------------------------------------------------------------------
/** Co-variant base vectors at centroid in current configuration evaluated
 *  at given local co-ordinate
 *  \param[in]  xi    local/parametric evaluation co-ordinate
 *                    \f$ \theta^1 \f$
 *  \param[out] base  co-variant current base vectors stored in (Jacobian) matrix;
 *                    the base vectors are given in global Cartesian components
 *                    (ie same co-ordinate system in which co-ordinates \f$ x^i = X^I + u^i\f$
 *                    of NODEs are given, subject to co-inciding Cart. bases 
 *                    \f$ e_i = e_I \f$) and stored as column vectors in
 *                    resulting matrix: \f$ [ (a^i)_\alpha ] \f$
 */
template< typename NODE, typename SFUN, typename MAT >
void beam::fem::ElementPlanarStatic<NODE,SFUN,MAT>::baseCur_( const VecLDim & xi,
                                                              Mat3x3& base ) const
{
    // nodal reference coordinates
    MatDimNF coorCur;  this->nodalCoordinatesCur( coorCur );
    // shape function derivatives
    MatLDimNF dphiDxi;  this->sFunGrad_( xi, dphiDxi );  

    // blank
    base.setZero();

    /// - axial director, tangent to centerline: \f$ a^0 = d x / d \theta \f$
    base.block( 0, 0, dim, 1 ) = coorCur * dphiDxi.transpose( ); 

    /// - bi-normal, normal to plane: \f$ a^2 = (0,0,1) \f$
    base( 2, 2 ) = 1.0;

    /// - unit normal to centerline: \f$ a^1 = (a^2 \times a^0) / |a^2 \times a^0| \f$
    base.col( 1 ) = 
        ( base.col( 2 ).cross( base.col( 0 ) ) ) / 
        base.col( 0 ).norm( );

    // quit
    return;
}

//------------------------------------------------------------------------------
/** Local/Parametric gradient of co-variant base vectors at centroid 
 *  in reference configuration evaluated at given local co-ordinate
 *  \param[in]  xi           local/parametric evaluation co-ordinate:
 *                           \f$ \theta^1 \f$
 *  \param[in]  base         co-variant reference base vectors stored in (Jacobian) matrix
 *                           \f$ [ (A^I)_\alpha ] \f$
 *  \param[out] baseGradLoc  local/parametric gradient \f$ \dot{(\bullet)} = d (\bullet) / d \theta \f$
 *                           of co-variant reference base vectors
 *                           with respect to local co-ordinate \f$ \theta=\theta^1 \f$;
 *                           the base vectors are given in global Cartesian components
 *                           and stored as column vectors in
 *                           resulting matrix: \f$ [ (\dot{A}^I)_\alpha ] \f$
 */
template< typename NODE, typename SFUN, typename MAT >
void beam::fem::ElementPlanarStatic<NODE,SFUN,MAT>::baseRefGradLoc_( const VecLDim& xi,
                                                                     const Mat3x3& base,
                                                                     Mat3x3& baseGradLoc ) const
{
    // get nodal coordinates
    MatDimNF coorRef;  this->nodalCoordinatesRef( coorRef );
    // shape function derivatives
    MatVecNFLDimLDim ddphi; this -> sFunGradGrad_( xi, ddphi );
    //MatLDimNF ddphi;  this->sFunGradGrad_( xi, ddphi );  

    // blank
    baseGradLoc.setZero();

    /// - axial director:  \f$ \dot{A}^0 = d^2 X / d \theta^2 \f$
    const VecDim tmp = coorRef * ddphi(0,0); 
    for ( unsigned d = 0; d < dim; d ++ ) baseGradLoc( d, 0 ) = tmp(d);

    /// - bi-normal: \f$ \dot{A}^2 = 0 \f$
    baseGradLoc( 2, 2 ) = 0.0;

    // length of A_0
    const double base0Length = base.col( 0 ).norm( );

    /** - normal to centerline: 
     *   \f[ 
     *       \dot{A}^1 = \frac{ \dot{A}^2 \times A^0 + A^2 \times \dot{A}^0 }{ |A^2 \times A^0 |}
     *                   + \frac{ d (1/|A^2 \times A^0|) }{ d \theta } A^2 \times A^0
     *                 = \frac{ A^2 \times \dot{A}^0 }{ |A^0| }
     *                   + \frac{ \dot{A}^0 \cdot A^0 }{ |A^0|^2 } A^1
     *   \f]
     */
    baseGradLoc.col( 1 ) = 
        base.col( 2 ).cross( baseGradLoc.col( 0 ) ) / base0Length
        - baseGradLoc.col( 0 ).dot( base.col( 0 ) ) / (base0Length * base0Length) * ( base.col( 1 ) );

    // quit
    return;
}

//------------------------------------------------------------------------------
/** Local/Parametric gradient of co-variant base vectors at centroid 
 *  in current configuration evaluated at given local co-ordinate
 *  \param[in]  xi           local/parametric evaluation co-ordinate:
 *                           \f$ \theta^1 \f$
 *  \param[in]  base         co-variant current base vectors stored in (Jacobian) matrix
 *                           \f$ [ (a^i)_\alpha ] \f$
 *  \param[out] baseGradLoc  local/parametric gradient \f$ \dot{(\bullet)} = d (\bullet) / d \theta \f$
 *                           of co-variant current base vectors
 *                           with respect to local co-ordinate \f$ \theta=\theta^1 \f$;
 *                           the base vectors are given in global Cartesian components
 *                           and stored as column vectors in
 *                           resulting matrix: \f$ [ (\dot{a}^i)_\alpha ] \f$
 */
template< typename NODE, typename SFUN, typename MAT >
void beam::fem::ElementPlanarStatic<NODE,SFUN,MAT>::baseCurGradLoc_( const VecLDim& xi,
                                                                     const Mat3x3& base,
                                                                     Mat3x3& baseGradLoc ) const
{
    // nodal coordinates
    MatDimNF coorCur;  this->nodalCoordinatesCur( coorCur );
    // shape function derivatives
    MatVecNFLDimLDim ddphi; this -> sFunGradGrad_( xi, ddphi );

    // blank
    baseGradLoc.setZero();

    /// - axial director: \f$ \dot{a}^0 = d^2 x / d \theta^2 \f$
    const VecDim tmp = coorCur * ddphi(0,0); 
    for ( unsigned d = 0; d < dim; d ++ ) baseGradLoc( d, 0 ) = tmp(d);

    /// - bi-normal: \f$ \dot{a}^2 = 0
    baseGradLoc( 2, 2 ) = 0.0;

    // length of base0
    const double base0length = base.col( 0 ).norm( );

    /// - normal to centerline (see baseRefGradLoc_ above)
    baseGradLoc.col( 1 ) = 
        base.col( 2 ).cross( baseGradLoc.col( 0 ) ) / base0length 
        - baseGradLoc.col( 0 ).dot( base.col( 0 ) ) / (base0length * base0length) * ( base.col( 1 ) );

    // quit
    return;
}

//------------------------------------------------------------------------------
/** Gradient of co-variant current base w.r.t. to current nodal co-ordinates
 *  \param[in]  xi           local/parametric evaluation co-ordinate:
 *                           \f$ \theta^1 \f$
 *  \param[in]  base         co-variant current base vectors stored in (Jacobian) matrix
 *                           \f$ [ (a^i)_\alpha ] \f$
 *  \param[out] baseGrad     gradient of co-variant current base vectors
 *                           with respect to current nodal co-ordinates \f$ (x^j)_K \f$;
 *                           and stored as column vectors in
 *                           resulting matrix: \f$ [ \ell^i{}_\alpha{}_j{}^K ] 
 *                           = [ \partial ( a^i{}_\alpha ) / \partial x^j{}_K ] \f$
 */
template< typename NODE, typename SFUN, typename MAT >
void beam::fem::ElementPlanarStatic<NODE,SFUN,MAT>::baseCurGradCur_( const VecLDim& xi,
                                                                     const Mat3x3& base,
                                                                     Mat3x3xDimxNF_& baseGrad ) const
{
    // shape function derivatives
    MatLDimNF dphiDxi;  this->sFunGrad_( xi, dphiDxi );  

    // blank
    for ( unsigned n=0; n<numNodesSN; ++n )
        for ( unsigned d=0; d<dim; ++d )
            baseGrad( d, n ).setZero();

    /** - gradient of axial director
     *    \f[
     *          \frac{\partial a^0_i}{\partial x^j_K} = \delta_{ij} \dot{\phi}^K
     *    \f]
     */
    for ( unsigned n=0; n<numNodesSN; ++n ) {
        for ( unsigned d=0; d<dim; ++d ) {
            baseGrad( d, n )( d, 0 ) = dphiDxi( 0, n );
        }
    }

    // gradient of thickness director
    const Vec3 bas0 = base.col( 0 );
    const Vec3 bas1 = base.col( 1 );
    const Vec3 bas2 = base.col( 2 );

    const double bas0Len     = bas0.norm( );
    const double bas0LenPow2 = bas0Len * bas0Len;

    /** - gradient of normal vector (in-plane)
     *    \f[
     *         \frac{\partial a^1}{\partial x^j_K} = 
     *               \frac{1}{|a^0|}   (a^2 \times \frac{\partial a^0}{\partial x^j_K} ) -
     *               \frac{1}{|a^0|^2} (a^0 \cdot  \frac{\partial a^0}{\partial x^j_K} ) a^1
     *    \f]
     */
    for ( unsigned n=0; n<numNodesSN; ++n ) {
        for ( unsigned d=0; d<dim; ++d ) {
            // \partial a^0 / \partial x
            const Vec3 bas0Grad = baseGrad( d, n ).col( 0 );
            baseGrad( d, n ).col( 1 ) = ( bas2.cross( bas0Grad ) ) / bas0Len;
            const double fact = bas0.dot( bas0Grad ) / bas0LenPow2;
            baseGrad( d, n ).col( 1 ) -= fact * bas1;
        }
    }

    // gradient of bi-normal
    // is zero
}

//------------------------------------------------------------------------------
/** Gradient of co-variant current base local gradient w.r.t. current nodal co-ordinates
 *  \param[in]  xi           local/parametric evaluation co-ordinate:
 *                           \f$ \theta^1 \f$
 *  \param[in]  base         co-variant current base vectors stored in (Jacobian) matrix
 *                           \f$ [ (a^i)_\alpha ] \f$
 *  \param[in]  baseGradLoc  co-variant current base vectors local gradient
 *                           \f$ [ (\dot{a}^i)_\alpha ] \f$
 *  \param[in]  baseGradCur  co-variant current base vectors current gradient
 *                           \f$ [ \ell^i{}_\alpha{}_j{}^K ] \f$
 *  \param[out] baseGradLocGradCur  gradient of co-variant current base vector local gradient
 *                           with respect to current nodal co-ordinates \f$ (x^j)_K \f$;
 *                           and stored as column vectors in
 *                           resulting matrix: \f$ [ \dot{\ell}^i{}_\alpha{}_j{}^K ] 
 *                           = [ \partial ( \dot{a}^i{}_\alpha ) / \partial ( x^j{}_K ) ] \f$
 */
template< typename NODE, typename SFUN, typename MAT >
void beam::fem::ElementPlanarStatic<NODE,SFUN,MAT>::baseCurGradLocGradCur_( const VecLDim & xi,
                                                                            const Mat3x3& base,
                                                                            const Mat3x3& baseGradLoc,
                                                                            const Mat3x3xDimxNF_& baseGradCur,
                                                                            Mat3x3xDimxNF_& baseGradLocGradCur ) const
{
    // shape function derivatives
    //MatLDimNF ddphiDxi;  this->sFunGradGrad_( xi, ddphiDxi );  
    MatVecNFLDimLDim ddphiDxi; this -> sFunGradGrad_( xi, ddphiDxi );

    // blank
    for ( unsigned n=0; n<numNodesSN; ++n )
        for ( unsigned d=0; d<dim; ++d )
            baseGradLocGradCur( d, n ).setZero();

    /** - gradient of axial director
     *    \f[
     *          \frac{\partial \dot{a}^0_i}{\partial x^j_K} = \delta_{ij} \ddot{\phi}^K
     *    \f]
     */
    for ( unsigned n=0; n<numNodesSN; ++n ) {
        for ( unsigned d=0; d<dim; ++d ) {
            baseGradLocGradCur( d, n )( d, 0 ) = ddphiDxi( 0, 0)( n );
        }
    }

    // gradient of thickness director
    const Vec3 bas0 = base.col( 0 );
    const Vec3 bas1 = base.col( 1 );
    const Vec3 bas2 = base.col( 2 );

    const double bas0Len     = bas0.norm( );
    const double bas0LenPow2 = bas0Len * bas0Len;
    const double bas0LenPow3 = bas0Len * bas0LenPow2;
    const double bas0LenPow4 = bas0LenPow2 * bas0LenPow2;

    const Vec3 bas0GradLoc = baseGradLoc.col( 0 );
    const Vec3 bas1GradLoc = baseGradLoc.col( 1 );
    const Vec3 bas2GradLoc = baseGradLoc.col( 2 );

    /** - gradient of the normal vector
     *    \f[
     *        \frac{\partial \dot{a}^1}{\partial x^j_K} =
     *             \frac{1}{|a^0|} (a^2 \times \partial_{x^j_K} \dot{a}^0)
     *           + \frac{ a^0 \cdot \partial_{x^j_K} a^0}{ |a^0|^3 } (\dot{a}^0 \times a^2)
     *           + \left( 2 \frac{(a^0\cdot \partial_{x^j_K} a^0)(a^0 \cdot \dot{a}^0)}{ |a^0|^4 }
     *              - \frac{\partial_{x_j^K} a^0 \cdot \dot{a}^0}{|a^0|^2}
     *              - \frac{a^0 \cdot \partial_{x_j^K} \dot{a}^0}{|a^0|^2} \right) a^1
     *           -  \frac{a^0 \cdot \dot{a}^0}{|a^0|^2} \partial_{x_j^K} a^1
     *    \f]
     */
    for ( unsigned n=0; n<numNodesSN; ++n ) {
        for ( unsigned d=0; d<dim; ++d ) {
            const Vec3 bas0GradCur = baseGradCur( d, n ).col( 0 );
            const Vec3 bas1GradCur = baseGradCur( d, n ).col( 1 );
            const Vec3 bas0GradLocGradCur = baseGradLocGradCur( d, n ).col( 0 );

            const double fact0  = bas0.dot( bas0GradCur );
            const double fact01 = fact0 / bas0LenPow3;
            const double fact02 = 2.0 * fact0 * bas0.dot( bas0GradLoc ) / bas0LenPow4;
            const double fact11 = bas0GradCur.dot( bas0GradLoc ) / bas0LenPow2;
            const double fact12 = bas0.dot( bas0GradLocGradCur ) / bas0LenPow2;
            const double fact2  = bas0.dot( bas0GradLoc ) / bas0LenPow2;

            baseGradLocGradCur( d, n ).col( 1 ) = 
                bas2.cross( bas0GradLocGradCur ) / bas0Len +
                fact01 * bas0GradLoc.cross( bas2 ) +
                ( fact02 - fact11 - fact12 ) * bas1 -
                fact2 * bas1GradCur;
        }
    }

    // gradient of bi-normal
    // is zero
}

//------------------------------------------------------------------------------
/** Strain resultants: membrane and bending contributions
 *  \param[in]  xi           Local coordinate at which kinematics are computed
 *                           \f$ \theta^1 \f$
 *  \param[in]  baseRef      Co-variant reference base vectors
 *                           \f$ [ (A^I)_\alpha ] \f$
 *  \param[in]  baseCur      Co-variant current base vectors
 *                           \f$ [ (a^i)_\alpha ] \f$
 *  \param[in]  baseRefGradLoc  Co-variant reference base vectors local gradient
 *                           \f$ [ (\dot{A}^I)_\alpha ] \f$
 *  \param[in]  baseCurGradLoc  Co-variant current base vectors local gradient
 *                           \f$ [ (\dot{a}^i)_\alpha ] \f$
 *  \param[out] strainMemb   Membrane strain resultant (axial strain)
 *                           \f$ \alpha \f$
 *  \param[out] strainBend   Bending strain resultant (curvature)
 *                           \f$ \beta \f$
 */
template< typename NODE, typename SFUN, typename MAT>
void beam::fem::ElementPlanarStatic<NODE,SFUN,MAT>::strainRes_( const VecLDim& xi, 
                                                                const Mat3x3& baseRef,
                                                                const Mat3x3& baseCur,
                                                                const Mat3x3& baseRefGradLoc,
                                                                const Mat3x3& baseCurGradLoc,
                                                                double& strainMemb,
                                                                double& strainBend ) const
{
    /** - membrane strain
     *    \f[
     *        \alpha = \frac{1}{2} ( a^0 \cdot a^0 - A^0 \cdot A^0 )
     *    \f]
     */
    strainMemb = 0.5 * (
        baseCur.col( 0 ).squaredNorm( ) -
        baseRef.col( 0 ).squaredNorm( )
        );

    /** - bending strain:
     *    \f[
     *        \beta = - ( a^1 \cdot \dot{a}^0 - A^1 \cdot \dot{A}^0 )
     *    \f]
     */
    strainBend = - (
        baseCur.col( 1 ).dot( baseCurGradLoc.col( 0 ) ) -
        baseRef.col( 1 ).dot( baseRefGradLoc.col( 0 ) )
        );

    // bye
    return;
}

//------------------------------------------------------------------------------
/** Gradient of strain resultants w.r.t. current nodal co-ordinates
 *
 *  \param[in]  xi           Local coordinate at which kinematics are computed
 *                           \f$ \theta^1 \f$
 *  \param[in]  baseRef      Co-variant reference base vectors
 *                           \f$ [ (A^I)_\alpha ] \f$
 *  \param[in]  baseCur      Co-variant current base vectors
 *                           \f$ [ (a^i)_\alpha ] \f$
 *  \param[in]  baseCurGradLoc  Co-variant current base vectors local gradient
 *                           \f$ [ (\dot{a}^i)_\alpha ] \f$
 *  \param[in]  baseCurGradCur  Gradient of co-variant base vector w.r.t. 
 *                           current nodal co-ordinates
 *                           \f$ [ \ell^i{}_\alpha{}_j{}^K ] \f$
 *  \param[in]  baseCurGradLocGradCur  Gradient of co-variant base vector local gradient w.r.t. 
 *                           current nodal co-ordinates
 *                           \f$ [ \dot{\ell}^i{}_\alpha{}_j{}^K ] \f$
 *  \param[out] strainMembGradCur   Membrane strain resultant gradient
 *                           \f$ [ \partial \alpha / \partial x^j{}_K ] \f$
 *  \param[out] strainBendGradCur   Bending strain resultant gradient
 *                           \f$ [ \partial \beta / \partial x^j{}_K ] \f$
 */
template< typename NODE, typename SFUN, typename MAT>
void beam::fem::ElementPlanarStatic<NODE,SFUN,MAT>::strainResGradCur_( const VecLDim& xi,
                                                                       const Mat3x3& baseRef,
                                                                       const Mat3x3& baseCur,
                                                                       const Mat3x3& baseCurGradLoc,
                                                                       const Mat3x3xDimxNF_& baseCurGradCur,
                                                                       const Mat3x3xDimxNF_& baseCurGradLocGradCur,
                                                                       MatDimNF& strainMembGradCur,
                                                                       MatDimNF& strainBendGradCur ) const
{
    /** - membrane strain gradient:
     *  \f[ 
     *       \frac{\partial \alpha}{\partial x^j_K} = 
     *                       a^0 \cdot \frac{\partial a^0}{\partial x^K_j}
     *  \f]
     */
    for ( unsigned d=0; d<dim; ++d ) {
        for ( unsigned n=0; n<numNodesSN; ++n ) {
            strainMembGradCur( d, n ) = (
                ( baseCur.col( 0 ) ).dot( baseCurGradCur( d, n ).col( 0 ) )
                );
        }
    }
    

    /** - bending strain gradient:
     *    \f[
     *        \frac{\partial \beta}{\partial x^j_K} =
     *                - \left( \dot{a}^0 \cdot \frac{\partial a^1}{\partial x^j_K} +
     *                         a^1 \cdot \frac{\partial \dot{a}^0}{\partial x^j_K}
     *                  \right)
     *    \f]
     */
    for ( unsigned d=0; d<dim; ++d ) {
        for ( unsigned n=0; n<numNodesSN; ++n ) {
            strainBendGradCur( d, n ) = - (
                ( baseCurGradLoc.col( 0 ) ).dot( baseCurGradCur( d, n ).col( 1 ) ) +
                ( baseCur.col( 1 ) ).dot( baseCurGradLocGradCur( d, n ).col( 0 ) )
                );
        }
    }

    // bye
    return;
}

//------------------------------------------------------------------------------
/** Hessian of strain resultants w.r.t. current nodal co-ordinates
 *
 *  \param[in]  xi           Local coordinate at which kinematics are computed
 *                           \f$ \theta^1 \f$
 *  \param[in]  baseCur      Co-variant current base vectors
 *                           \f$ [ (a^i)_\alpha ] \f$
 *  \param[in]  baseRef      Co-variant reference base vectors
 *                           \f$ [ (A^I)_\alpha ] \f$
 *  \param[in]  baseCurGradLoc  Co-variant current base vectors local gradient
 *                           \f$ [ (\dot{a}^i)_\alpha ] \f$
 *  \param[in]  baseCurGradCur  Gradient of co-variant base vector w.r.t. 
 *                           current nodal co-ordinates
 *                           \f$ [ \ell^i{}_\alpha{}_j{}^K ] \f$
 *  \param[in]  baseCurGradLocGradCur  Gradient of co-variant base vector local gradient w.r.t. 
 *                           current nodal co-ordinates
 *                           \f$ [ \dot{\ell}^i{}_\alpha{}_j{}^K ] \f$
 *  \param[out] strainMembGradGradCur   Membrane strain resultant Hessian
 *                           \f$ [ \partial^2 \alpha / \partial x^i{}_K \partial x^j{}_L ] \f$
 *  \param[out] strainBendGradGradCur   Bending strain resultant Hessian
 *                           \f$ [ \partial^2 \beta / \partial x^i{}_K \partial x^i{}_L ] \f$
 */
template< typename NODE, typename SFUN, typename MAT>
void beam::fem::ElementPlanarStatic<NODE,SFUN,MAT>::strainResGradGradCur_( const VecLDim& xi,
                                                                           const Mat3x3& baseRef,
                                                                           const Mat3x3& baseCur,
                                                                           const Mat3x3& baseCurGradLoc,
                                                                           const Mat3x3xDimxNF_& baseCurGradCur,
                                                                           const Mat3x3xDimxNF_& baseCurGradLocGradCur,
                                                                           ElemMat_& strainMembGradGradCur,
                                                                           ElemMat_& strainBendGradGradCur ) const
{
    /// - membrane strain: \f$ \frac{\partial^2 \alpha}{\partial x^i_K x^j_L} = \dots \f$
    // membrane strain Hessian
    // [ \partial^2 \alpha / \partial (x^i)_K \partial (x^j)_L]
    //     = [ \ell^a{}_1{}_i{}^K ] . [ \ell^a{}_1{}_j{}^L ]
    for ( unsigned iDof=0; iDof<dof; ++iDof ) {
        for ( unsigned kNod=0; kNod<numNodesSN; ++kNod ) {
            const unsigned ikDofNod = kNod*dof + iDof;

            for ( unsigned jDof=0; jDof<dof; ++jDof ) {
                for ( unsigned lNod=0; lNod<numNodesSN; ++lNod ) {
                    const unsigned jlDofNod = lNod*dof + jDof;
                    
                    strainMembGradGradCur( ikDofNod, jlDofNod ) = (
                        ( baseCurGradCur( iDof, kNod ).col( 0 ) ).dot( baseCurGradCur( jDof, lNod ).col( 0 ) )
                        );

                }
            }

        }
    }

    /// - bending strain: \f$ \frac{\partial^2 \beta}{\partial x^i_K x^j_L} = \dots \f$
    // bending strain Hessian
    const Vec3 bas0 = baseCur.col( 0 );
    const Vec3 bas1 = baseCur.col( 1 );
    const Vec3 bas2 = baseCur.col( 2 );

    const double bas0Len = bas0.norm( );
    const double bas0LenPow2 = bas0Len * bas0Len;
    const double bas0LenPow3 = bas0Len * bas0LenPow2;

    const Vec3 bas0GradLoc = baseCurGradLoc.col( 0 );

    for ( unsigned iDof=0; iDof<dof; ++iDof ) {
        for ( unsigned kNod=0; kNod<numNodesSN; ++kNod ) {
            const unsigned ikDofNod = kNod*dof + iDof;

            const Vec3 bas0GradCur_ik = baseCurGradCur( iDof, kNod ).col( 0 );
            const Vec3 bas1GradCur_ik = baseCurGradCur( iDof, kNod ).col( 1 );
            const Vec3 bas2xBas0GradCur_ik = bas2.cross( bas0GradCur_ik );

            const Vec3 bas0GradLocGradCur_ik = baseCurGradLocGradCur( iDof, kNod ).col( 0 );

            for ( unsigned jDof=0; jDof<dof; ++jDof ) {
                for ( unsigned lNod=0; lNod<numNodesSN; ++lNod ) {
                    const unsigned jlDofNod = lNod*dof + jDof;

                    const Vec3 bas0GradCur_jl = baseCurGradCur( jDof, lNod ).col( 0 );
                    const Vec3 bas1GradCur_jl = baseCurGradCur( jDof, lNod ).col( 1 );

                    const Vec3 bas0GradLocGradCur_jl = baseCurGradLocGradCur( jDof, lNod ).col( 0 );

                    const double fact0_jl = -bas0.dot( bas0GradCur_jl ) / bas0LenPow3;
                    const double fact1_ik = bas0.dot( bas0GradCur_ik ) / bas0LenPow2;
                    const double fact2_ikjl = bas0GradCur_ik.dot( bas0GradCur_jl ) / bas0LenPow2;
                    const double fact3_ikjl = 2.0 * bas0Len * fact0_jl * fact1_ik;

                    strainBendGradGradCur( ikDofNod, jlDofNod ) = - (
                        bas1GradCur_ik.dot( bas0GradLocGradCur_jl ) +
                        bas1GradCur_jl.dot( bas0GradLocGradCur_ik ) +
                        fact0_jl * bas0GradLoc.dot( bas2xBas0GradCur_ik ) -
                        fact1_ik * bas0GradLoc.dot( bas1GradCur_jl ) +
                        ( -fact2_ikjl + fact3_ikjl ) * bas0GradLoc.dot( bas1 )
                        );

                }
            }

        }
    }

    // bye
    return;
}


//------------------------------------------------------------------------------
/** Stress resultants
 *
 *  \param[in]  xi           Local coordinate at which kinematics are computed
 *                           \f$ \theta^1 \f$
 *  \param[in]  baseRef      Co-variant reference base vectors
 *                           \f$ [ (A^I)_\alpha ] \f$
 *  \param[in]  strainMemb   Membrane strain resultant (axial strain)
 *                           \f$ \alpha \f$
 *  \param[in]  strainBend   Bending strain resultant (curvature)
 *                           \f$ \beta \f$
 *  \param[out] stressMemb   Membrane stress resultant (axial/normal force)
 *                           \f$ n \f$
 *  \param[out] stressBend   Bending stress resultant (torque)
 *                           \f$ m \f$
 */
template< typename NODE, typename SFUN, typename MAT>
void beam::fem::ElementPlanarStatic<NODE,SFUN,MAT>::stressRes_( const VecLDim& xi,
                                                                const Mat3x3& baseRef,
                                                                const double& strainMemb,
                                                                const double& strainBend,
                                                                double& stressMemb,
                                                                double& stressBend ) const
{
    // elastic moduli
    double elastModMemb, elastModBend;
    this->elasticityRes( baseRef, elastModMemb, elastModBend );

    /// - membrane stress resultant \f$ n = \bar{EA} \alpha \f$
    stressMemb = elastModMemb * strainMemb;

    /// - bending stress resultant  \f$ m = \bar{EI} \beta \f$
    stressBend = elastModBend * strainBend;

    // bye
    return;
}

//------------------------------------------------------------------------------
/** Give Jacobian determinant
 *
 *  \param[in]  xi              Local co-ordinate \f$[\theta^\alpha]\f$
 *  \param[in]  currentConfig   If true, axial metric W.R.T. current configuration,
 *                                 otherwise W.R.T. reference configuration
 */
template< typename NODE, typename SFUN, typename MAT >
double beam::fem::ElementPlanarStatic<NODE,SFUN,MAT>::jac_( const VecLDim& xi,
                                                            const bool currentConfig ) const
{
    // nodal coordinates
    MatDimNF coorRef;
    if ( currentConfig )  this->nodalCoordinatesCur( coorRef );
    else                  this->nodalCoordinatesRef( coorRef );
    // shape function derivatives
    MatLDimNF dphiDxi;  this->sFunGrad_( xi, dphiDxi );
    // tangential vector = a^0
    const VecDim base0 = ( coorRef * dphiDxi.transpose( ) ).col( 0 );
    // return its Euclidean norm
    // which is the _simplified_ Jacobian
    return base0.norm( );
}

//------------------------------------------------------------------------------
/** Computation of the internal strain energy. Evaluation of kernel function at
 *  given local coordinate, multiplication with a weight and addition on the
 *  energy scalar: energy =  \f$ \int W(F(\xi)) \det G dX \f$. The result
 *  is stored on member #energy_.
 *  \f[
 *        E = \frac{1}{2} \int m \alpha + n \beta 
 *  \f]
 *  \param[in] xi      Local coordinate at which this function is evaluated
 *  \param[in] weight  Corresponding quadrature weight
 */
template<typename NODE, typename SFUN, typename MAT>
void beam::fem::ElementPlanarStatic<NODE,SFUN,MAT>::strainEnergyIntegrand( const VecLDim& locCoor,
                                                                           const double& weight ) 
{
    // get kinematic quantities
    Mat3x3 baseRef;        this -> baseRef_(        locCoor, baseRef );
    Mat3x3 baseCur;        this -> baseCur_(        locCoor, baseCur );
    Mat3x3 baseRefGradLoc; this -> baseRefGradLoc_( locCoor, baseRef, baseRefGradLoc );
    Mat3x3 baseCurGradLoc; this -> baseCurGradLoc_( locCoor, baseCur, baseCurGradLoc );

    // strain resultants and gradient
    double strainMemb, strainBend;
    this -> strainRes_( locCoor, baseRef, baseCur, baseRefGradLoc, baseCurGradLoc,
                        strainMemb, strainBend );

    // stress resultants
    double stressMemb, stressBend;
    this -> stressRes_( locCoor, strainMemb, strainBend, stressMemb, stressBend );

    // pass back the strain energy density
    const double det = baseRef.col( 0 ).norm( );
    (this->energy_) += 0.5 * ( strainMemb * stressMemb + 
                               strainBend * stressBend ) * det * weight;

    return;
}

//------------------------------------------------------------------------------
/** Computation of internal force due to non-equilibrium stress. Evaluates the kernel
 *  function at given quadrature point xi, multiplies the result with the weight
 *  and passes it to the nodes.
 *  \f[
 *       F[k d + i] = \int m \frac{\partial \alpha}{\partial x^i_K} +
 *                         n \frac{\partial \beta }{\partial x^i_K}
 *  \f]
 *
 *  \param[in] xi      Local coordinate at which this function is evaluated
 *  \param[in] weight  Corresponding quadrature weight
 */
template<typename NODE, typename SFUN, typename MAT >
void beam::fem::ElementPlanarStatic<NODE,SFUN,MAT>::internalForceIntegrand( const VecLDim& locCoor, 
                                                                            const double& weight ) 
{
    // get kinematic quantities
    Mat3x3 baseRef;
    this->baseRef_( locCoor, baseRef );
    Mat3x3 baseCur;
    this->baseCur_( locCoor, baseCur );
    Mat3x3 baseRefGradLoc;
    this->baseRefGradLoc_( locCoor, baseRef, baseRefGradLoc );
    Mat3x3 baseCurGradLoc;
    this->baseCurGradLoc_( locCoor, baseCur, baseCurGradLoc );
    Mat3x3xDimxNF_ baseCurGradCur;
    this->baseCurGradCur_( locCoor, baseCur, baseCurGradCur );
    Mat3x3xDimxNF_ baseCurGradLocGradCur;
    this->baseCurGradLocGradCur_( locCoor, baseCur, baseCurGradLoc, baseCurGradCur,
                                  baseCurGradLocGradCur );

    // strain resultants and gradient
    double strainMemb, strainBend;
    this->strainRes_( locCoor, baseRef, baseCur, baseRefGradLoc, baseCurGradLoc,
                      strainMemb, strainBend );
    MatDimNF strainMembGradCur, strainBendGradCur;
    this->strainResGradCur_( locCoor, baseRef, baseCur,
                             baseCurGradLoc, baseCurGradCur, baseCurGradLocGradCur,
                             strainMembGradCur, strainBendGradCur );

    // stress resultants
    double stressMemb, stressBend;
    this->stressRes_( locCoor, baseRef, strainMemb, strainBend,
                      stressMemb, stressBend );

    // compute matrix of nodal forces 
    const double det = baseRef.col( 0 ).norm( );
    const MatDofNF nodalForces = (
        stressMemb * strainMembGradCur.block( 0, 0, dof, numNodesSN ) +
        stressBend * strainBendGradCur.block( 0, 0, dof, numNodesSN )
        ) * det * weight;

    // pass each nodal force to the nodes
    for ( unsigned n=0; n<numNodesSN; ++n ) {
        nodes_[ n ]->addToForce( - nodalForces.col( n ) );
    }

    return;
}

//------------------------------------------------------------------------------
/** Compute nodal forces due to an applied body force. This function is given
 *  the coordinates and weight of a quadrature point. Moreover, a function object
 *  is passed with the operator of type VecDof = func(VecDim), where the 
 *  argument refers to the reference coordinate of the integration point.
 *
 *  \tparam FUNC                 Type of function object
 *  \param[in] xi                Local coordinate of quadrature point
 *  \param[in] weight            Weight of quadrature point
 *  \param[in] func              Function defining a load field along beam axis
 *  \param[in] factor            Scalar multiplier
 *  \param[in] currentConfig     If true, a current force field is considered
 *                                  (Evaluation occurs in current configuration of beam axis,
 *                                  and the load field is subject to current beam length too)
 *  \param[in] isAreaIntegrated  If true, the force field is a line load (N/m)
 *                                  otherwise it's a density (N/m^3)
 */
template< typename NODE, typename SFUN, typename MAT >
template< typename FUNC >
void beam::fem::ElementPlanarStatic<NODE,SFUN,MAT>::bodyForceDetailed( const VecLDim& xi,
                                                                       const double& weight,
                                                                       FUNC func,
                                                                       const double& factor,
                                                                       const bool currentConfig,
                                                                       const bool isAreaIntegrated )
{
    // get global coordinate of xi
    MatDimNF verts;
    if ( currentConfig ) this -> nodalCoordinatesCur( verts );
    else                 this -> nodalCoordinatesRef( verts );
    VecNF phi;  this -> sFun_( xi, phi );
    const VecDim coor = verts * phi;

    // compute value of body force density function [N/m^3] or readily integrated line load [N/m]
    const VecDof forceDens = func( coor );

    // get jacobian of this element
    double det = this->jac_( xi, currentConfig );
    
    // if not readily integrated over cross section area
    if ( not isAreaIntegrated )
        det *= area_;

    // Compute the nodal forces 
    for ( unsigned m=0; m<numNodesSN; ++m ) {
        const VecDof force = phi( m ) * forceDens * det * weight * factor;
        nodes_[ m ]->addToForce( force );
    }
    return;
}

//------------------------------------------------------------------------------
/** Compute nodal forces due to an applied body force. This function is given
 *  the coordinates and weight of a quadrature point. Moreover, a function object
 *  is passed with the operator of type VecDof = func(VecDim), where the 
 *  argument refers to the reference coordinate of the integration point.
 *
 *  \tparam FUNC       Type of function object
 *  \param[in] xi      Local coordinate of quadrature point
 *  \param[in] weight  Weight of quadrature point
 *  \param[in] func    Function defining a material force density field [N/m^3] along beam element length
 *  \param[in] factor  Scalar multiplier
 */
template< typename NODE, typename SFUN, typename MAT >
template< typename FUNC >
void beam::fem::ElementPlanarStatic<NODE,SFUN,MAT>::bodyForce( const VecLDim& xi,
                                                               const double& weight,
                                                               FUNC func, 
                                                               const double& factor )
{
    this->bodyForceDetailed<FUNC>( xi, weight, func, factor, false, false );
    return;
}

//------------------------------------------------------------------------------
/** Compute equivalent nodal forces due to concentrated forces and torques
 *
 *  \param[in]  xi       Local co-ordinate on beam centroid at which
 *                       concentrated loads are applied. This is usually
 *                       0.0 or 1.0, i.e. left or right end of topological beam
 *                       element.
 *  \param[in]  forces   Concentrated force vector: F_x, F_y, F_z
 *  \param[in]  torques  Concentrated moment vector: M_xx, M_yy, M_zz
 *  \param[in]  factor   (Global) load factor
 */
template< typename NODE, typename SFUN, typename MAT >
void beam::fem::ElementPlanarStatic<NODE,SFUN,MAT>::concentratedLoad( const VecLDim& coorLoc,
                                                                      const Vec3& forces,
                                                                      const Vec3& torques,
                                                                      const double& factor )
{
    // get kinematic quantities

    // shape function
    VecNF shpFct;  this->sFun_( coorLoc, shpFct );

    // store derivative of centroidal with respect to nodal co-ordinates
    // HINT: Could be done more efficiently
    eigenX::MatrixSd< 3, numNodesSN > coorGradCur;
    coorGradCur.setZero();
    for ( unsigned n=0; n<numNodesSN; ++n ) {
        for ( unsigned d=0; d<dim; ++d ) {
            coorGradCur( d, n ) = shpFct( n );
        }
    }

    //Mat3x3 baseRef;
    //this->baseRef_( coorLoc, baseRef );
    Mat3x3 baseCur;
    this->baseCur_( coorLoc, baseCur );
    Mat3x3xDimxNF_ baseCurGradCur;
    this->baseCurGradCur_( coorLoc, baseCur, baseCurGradCur );

    // convenience
    const Vec3 bas0 = baseCur.col( 0 );
    const Vec3 bas1 = baseCur.col( 1 );
    const Vec3 bas2 = baseCur.col( 2 );

    const double bas0Len = bas0.norm( );
    const double bas0LenPow2 = bas0Len * bas0Len;
    const double bas0LenPow3 = bas0LenPow2 * bas0Len;

    // Compute the nodal forces 
    for ( unsigned kNod=0; kNod<numNodesSN; ++kNod ) {
        VecDof nodalForce;
        for ( unsigned iDof=0; iDof<dof; ++iDof ) {
            nodalForce( iDof ) = (
                // force in global x-dir OR force in global y-dir
                forces( iDof ) * coorGradCur( iDof, kNod ) +
                // torque around z-axis
                torques( 2 ) / bas0Len * bas1.dot( baseCurGradCur( iDof , kNod ).col( 0 ) )
                );
        }

        // store onto nodal force vectors
        nodes_[ kNod ]->addToForce( factor * nodalForce );
    }

    // compute tangent
    ElemMat_ tang;
    tang.setZero();
    for ( unsigned iDof=0; iDof<dof; ++iDof ) {
        for ( unsigned kNod=0; kNod<numNodesSN; ++kNod ) {
            const unsigned ikDofNod = kNod*dof + iDof;

            const Vec3 bas0GradCur_ik = baseCurGradCur( iDof, kNod ).col( 0 );
            const Vec3 bas1GradCur_ik = baseCurGradCur( iDof, kNod ).col( 1 );

            for ( unsigned jDof=0; jDof<dof; ++jDof ) {
                for ( unsigned lNod=0; lNod<numNodesSN; ++lNod ) {
                    const unsigned jlDofNod = lNod*dof + jDof;

                    const Vec3 bas0GradCur_jl = baseCurGradCur( jDof, lNod ).col( 0 );
                    const Vec3 bas1GradCur_jl = baseCurGradCur( jDof, lNod ).col( 1 );

                    const double fact1_jl = bas0.dot( bas0GradCur_jl ) / bas0LenPow3;

                    tang( ikDofNod, jlDofNod ) = torques( 2 ) * (
                        bas1GradCur_jl.dot( bas0GradCur_ik ) / bas0Len -
                        fact1_jl * bas1.dot( bas0GradCur_ik )
                        );

                }
            }

        }
    }
    // store on element stiffness matrix
    // WARNING: This is likely to fail for dynamic simulations
    elemStiff_ -= factor * tang;

    // done
    return;
}

//------------------------------------------------------------------------------
/** The integral kernel for the computation of the element stiffness matrix. 
 *  Evaluates the kernel function at the local coordinate xi, multiplies the 
 *  result with the weight and stores the result in a matrix.
 *  \f[
 *      K[K d + i, L d + j] =  \int
 *        \frac{\partial \alpha}{\partial x_K^i} \bar{EA} 
 *        \frac{\partial \alpha}{\partial x_L^j}          +
 *        \frac{\partial \beta }{\partial x_K^i} \bar{EI} 
 *        \frac{\partial \beta }{\partial x_L^j}          +
 *        m  \frac{\partial^2 \alpha}{\partial x_K^i \partial x_L^j} +
 *        n  \frac{\partial^2 \beta}{\partial x_K^i \partial x_L^j} 
 *  \f]
 *
 *  \param[in] xi     Local coordinate at which this function is evaluated
 *  \param[in] weight Corresponding quadrature weight
 */
template< typename NODE, typename SFUN, typename MAT >
void beam::fem::ElementPlanarStatic<NODE,SFUN,MAT>::stiffnessIntegrand( const VecLDim & coorLoc,
                                                                        const double  & weight )
{
    // get kinematic quantities
    Mat3x3 baseRef;
    this->baseRef_( coorLoc, baseRef );
    Mat3x3 baseCur;
    this->baseCur_( coorLoc, baseCur );
    Mat3x3 baseRefGradLoc;
    this->baseRefGradLoc_( coorLoc, baseRef, baseRefGradLoc );
    Mat3x3 baseCurGradLoc;
    this->baseCurGradLoc_( coorLoc, baseCur, baseCurGradLoc );
    Mat3x3xDimxNF_ baseCurGradCur;
    this->baseCurGradCur_( coorLoc, baseCur, baseCurGradCur );
    Mat3x3xDimxNF_ baseCurGradLocGradCur;
    this->baseCurGradLocGradCur_( coorLoc, baseCur, baseCurGradLoc, baseCurGradCur,
                                  baseCurGradLocGradCur );

    // strain resultants and gradients
    double strainMemb, strainBend;
    this->strainRes_( coorLoc, baseRef, baseCur, baseRefGradLoc, baseCurGradLoc,
                      strainMemb, strainBend );
    MatDimNF strainMembGradCur, strainBendGradCur;
    this->strainResGradCur_( coorLoc, baseRef, baseCur,
                             baseCurGradLoc, baseCurGradCur, baseCurGradLocGradCur,
                             strainMembGradCur, strainBendGradCur );
    ElemMat_ strainMembGradGradCur, strainBendGradGradCur;
    this->strainResGradGradCur_( coorLoc, baseRef, baseCur,
                                 baseCurGradLoc, baseCurGradCur, baseCurGradLocGradCur,
                                 strainMembGradGradCur, strainBendGradGradCur );

    // stress resultants
    double stressMemb, stressBend;
    this->stressRes_( coorLoc, baseRef,
                      strainMemb, strainBend,
                      stressMemb, stressBend );

    // elastic moduli
    double elastModMemb, elastModBend;
    this->elasticityRes( baseRef, elastModMemb, elastModBend );

    // jacobian determinant
    const double det = ( baseRef.col( 0 ).norm( ) );

    for ( unsigned iDof=0; iDof<dof; ++iDof ) {
        for ( unsigned kNod=0; kNod<numNodesSN; ++kNod ) {
            const unsigned ikDofNod = kNod*dof + iDof;

            for ( unsigned jDof=0; jDof<dof; ++jDof ) {
                for ( unsigned lNod=0; lNod<numNodesSN; ++lNod ) {
                    const unsigned jlDofNod = lNod*dof + jDof;
                    
                    elemStiff_( ikDofNod, jlDofNod ) += (
                        // so-called elastic and initial displacement
                        strainMembGradCur( iDof, kNod ) * elastModMemb * strainMembGradCur( jDof, lNod ) +
                        strainBendGradCur( iDof, kNod ) * elastModBend * strainBendGradCur( jDof, lNod ) +
                        // so-called geometric
                        stressMemb * strainMembGradGradCur( ikDofNod, jlDofNod ) +
                        stressBend * strainBendGradGradCur( ikDofNod, jlDofNod )
                        ) * det * weight;

                }
            }
        }
    }

    return;
}

//------------------------------------------------------------------------------
/** integrate volume of beam element in reference configuration
 *
 *  \param[in] xi     Local coordinate at which this function is evaluated
 *  \param[in] weight Corresponding quadrature weight
 */
template< typename NODE, typename SFUN, typename MAT >
void beam::fem::ElementPlanarStatic<NODE,SFUN,MAT>::volumeIntegrand( const VecLDim & xi,
                                                                     const double & weight )
{
    // copy [0,1] -> [0,1]
    const VecLDim coorLoc =  xi;

    // jacobian determinant
    const double det = this->jac_( coorLoc, false );

    // add contribution of quadrature point to volume
    volume_ += area_ * det * weight;
    return;
}

//------------------------------------------------------------------------------
/** Computation of internal force due to non-equilibrium stress. Evaluates the kernel
 *  function at given quadrature point xi, multiplies the result with the weight
 *  and passes it to the nodes.
 *
 *  \param[in]  locCoor      Local coordinate at which this function is evaluated
 *  \param[out] stressMemb   Membrane stress resultant (axial/normal force)
 *                           \f$ n \f$
 *  \param[out] stressBend   Bending stress resultant (torque)
 *                           \f$ m \f$
 */
template< typename NODE, typename SFUN, typename MAT >
void beam::fem::ElementPlanarStatic<NODE,SFUN,MAT>::giveStressResultants( const VecLDim & locCoor,
                                                                          double & stressMemb,
                                                                          double & stressBend ) const
{
    // get kinematic quantities
    Mat3x3 baseRef;
    this->baseRef_( locCoor, baseRef );
    Mat3x3 baseCur;
    this->baseCur_( locCoor, baseCur );
    Mat3x3 baseRefGradLoc;
    this->baseRefGradLoc_( locCoor, baseRef, baseRefGradLoc );
    Mat3x3 baseCurGradLoc;
    this->baseCurGradLoc_( locCoor, baseCur, baseCurGradLoc );
    Mat3x3xDimxNF_ baseCurGradCur;
    this->baseCurGradCur_( locCoor, baseCur, baseCurGradCur );
    Mat3x3xDimxNF_ baseCurGradLocGradCur;
    this->baseCurGradLocGradCur_( locCoor, baseCur, baseCurGradLoc, baseCurGradCur,
                                  baseCurGradLocGradCur );

    // strain resultants and gradient
    double strainMemb, strainBend;
    this->strainRes_( locCoor, baseRef, baseCur, baseRefGradLoc, baseCurGradLoc,
                      strainMemb, strainBend );
    MatDimNF strainMembGradCur, strainBendGradCur;
    this->strainResGradCur_( locCoor, baseRef, baseCur,
                             baseCurGradLoc, baseCurGradCur, baseCurGradLocGradCur,
                             strainMembGradCur, strainBendGradCur );

    // stress resultants
    this->stressRes_( locCoor, baseRef, strainMemb, strainBend,
                      stressMemb, stressBend );
    return;
}

