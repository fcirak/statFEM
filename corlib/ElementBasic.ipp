// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//
//                 Computational Structural Mechanics Lab
//                         University of Cambridge
//
// This software is copyrighted and all rights are retained by the University of
// Cambridge. All use, disclosure, and/or reproduction of any part not expressly
// authorised by the University of Cambridge is prohibited. (C) 2010
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

//! @file ElementBasic.ipp

#include <algorithm>
#include <Eigen/Core>
#include <corlib/misc.hpp>
#include <corlib/eigenX.hpp>

//------------------------------------------------------------------------------
//! Get nodal coordinates
template<typename NODE, typename SFUN>
void corlib::ElementBasic<NODE,SFUN>::nodalCoordinates_( MatDimNN & X ) const
{
    for ( unsigned j = 0; j < numNodes; j ++ ) {
        X.col( j ) = nodes_[j] -> giveCoordinates();
    }
    return;
}

//------------------------------------------------------------------------------
/** Read node IDs from stream and set node pointers from container
 *  \param  inp     Input stream
 *  \param  nodes   Container with node pointers
 */
template<typename NODE, typename SFUN>
std::istream & corlib::ElementBasic<NODE,SFUN>:: 
readSelf( std::istream & inp, const std::vector<Node*> & nodes ) 
{
    // read node IDs and set pointers
    for ( unsigned v = 0; v < numNodes; v ++ ) {
        unsigned vId;
        inp >> vId; 
        nodes_[v] = nodes[vId];
        if ( inp.peek() == ',' ) inp.ignore();
    }
    return inp;
}

//------------------------------------------------------------------------------
/** Apply a given functor to all nodes by means of the for_each algorithm
 *  (const version)
 *  \tparam OP       Type of node functor
 *  \param[in] op    Functor applicable to NODE*
 *  \return          Given functor op
 */
template<typename NODE, typename SFUN>
template<typename OP> 
OP corlib::ElementBasic<NODE,SFUN>::iterateOverNodes(OP op) const
{
    return std::for_each( nodes_.begin(), nodes_.end(), op );
}

//------------------------------------------------------------------------------
/** Apply a given functor to all nodes by means of the for_each algorithm
 *  (non-const version)
 *  \tparam OP       Type of node functor
 *  \param[in] op    Functor applicable to NODE*
 *  \return          Given functor op
 */
template<typename NODE, typename SFUN>
template<typename OP> 
OP corlib::ElementBasic<NODE,SFUN>::iterateOverNodes(OP op) 
{
    return std::for_each( nodes_.begin(), nodes_.end(), op );
}

//------------------------------------------------------------------------------
/** Evaluates the shape function at the given coordinate xi
 *  \param[in]  xi   Local cooordinate at which the functions is evaluated
 *  \param[out] phi  Result of all shape functions at xi 
 */
template<typename NODE, typename SFUN>
void corlib::ElementBasic<NODE,SFUN>::sfun( const VecLDim & xi, 
                                            VecNN & phi ) const
{
    ShapeFun shapeFun;
    shapeFun.evaluate( xi, phi );
}

//------------------------------------------------------------------------------
/** Compute the gradient w.r.t. xi of the shape functions at the given xi:
 * \f$ \frac{d \varphi}{d \xi}|_ {\xi} \f$
 * \param[in]  xi       Local coordinate at which the gradient is evaluated
 * \param[out] dphiDxi  dphiDxi[i,K] \f$ = d \varphi^K / d \xi_i \f$
 */
template<typename NODE, typename SFUN>
void corlib::ElementBasic<NODE,SFUN>::sfunGrad( const VecLDim & xi, 
                                                MatLDimNN & dphiDxi ) const
{
    ShapeFun shapeFun;
    shapeFun.evaluateGradient( xi, dphiDxi );
}

//------------------------------------------------------------------------------
/** Compute the 2nd derivative w.r.t xi of the shape functions at the given xi:
 *
 * \param[in]  xi         Local coordinate at which the gradient is evaluated
 *                        \f$\xi^\alpha\f$
 * \param[out] ddPhiDDXi  2nd derivative of shape functions
 *                        \f$\varphi^K{}_{,\alpha\beta}\f$
 *                        \f$=\partial^2\varphi^K / \partial\xi^\alpha\partial\xi^\beta\f$.
 */
template<typename NODE, typename SFUN>
void corlib::ElementBasic<NODE,SFUN>::sfunHess( const VecLDim & xi,
                                                MatVecNNLDimLDim & ddPhiDDXi ) const
{
    ShapeFun shapeFun;
    shapeFun.evaluateHessian( xi, ddPhiDDXi );
}

//------------------------------------------------------------------------------
/** Computes the determinant of the element's Jacobi matrix at given xi.
 *  At first the covariant basis vectors 
 *  \f[
 *        g_\alpha = \frac{ \partial x }{\partial \xi_\alpha }
 *  \f]
 *  are computed. These form the columns of some matrix \f$ G \f$ whose
 *  (generalised) determinant is computed and returned.
 *  \param[in]  xi   Local coordinate at which det G is computed
 *  \return          The determinant \f$ |G(\xi)| \f$
 */
template<typename NODE, typename SFUN>
double corlib::ElementBasic<NODE,SFUN>::jacobian( const VecLDim & xi ) const
{
    // nodal coordinates
    MatDimNN X; 
    this -> nodalCoordinates_( X );

    // shape function derivatives
    MatLDimNN dphiDxi; 
    this -> sfunGrad( xi, dphiDxi );

    // covariant basis
    const eigenX::MatrixSd<dim,localDim> coG 
        = X * dphiDxi.transpose( );

    // return its determinant (metric of the mapping)
    const double detJ = corlib::metric( coG );

    FTL_VERIFY_DESCRIPTIVE( detJ > 0., "Negative Jacobian." );

    return detJ;
}


//------------------------------------------------------------------------------
/** The derivatives with respect to the reference coordinates evaluated at the
 *  given local coordinate \f$ xi \f$ are computed. For this purpose, first
 *  the co-variant basis vectors
 *  \f[
 *      g_\alpha[i] = \frac{\partial x_i}{\partial \xi^\alpha},
 *  \f]
 *  and then the contra-variant basis vectors are computed
 *  \f[
 *      g^\alpha[i] = \frac{\partial \xi^\alpha}{\partial x_i}
 *  \f]
 *
 *  The derivatives of the shape function with respect to the
 *  global coordinates \f$ x \f$ are now
 *  \f[
 *         \frac{\partial \phi^K}{\partial x} = 
 *            g^\alpha \frac{\partial \phi^K}{\partial \xi^\alpha}
 *  \f]
 *  \param[in]  xi       Local coordinate at which the derivates are computed
 *  \param[out] dphiDx   dphiDx[i,K] \f$= d \varphi^K / d x_i \f$
 *  \return              The Jacobi determinant \f$ \sqrt{|g(\xi)|} \f$
 */
template<typename NODE, typename SFUN>
double corlib::ElementBasic<NODE,SFUN>::globalDerivatives( const VecLDim & xi, 
                                                           MatDimNN & dphiDx ) const
{
    // nodal coordinates
    MatDimNN X; 
    this -> nodalCoordinates_( X );

    // shape function derivatives
    MatLDimNN dphiDxi; 
    this -> sfunGrad( xi, dphiDxi );

    // covariant basis
    const eigenX::MatrixSd<dim,localDim> coBase 
        = X * dphiDxi.transpose( ); 

    // contravariant basis (inverse of the covariant basis)
    eigenX::MatrixSd<dim,localDim> contraBase;
    // specialised call for determinant and contra-variant basis
    const double detJ = 
        corlib::ContraVariantBasis<dim,localDim>::compute( coBase, contraBase );

    FTL_VERIFY_DESCRIPTIVE( detJ > 0., "Negative Jacobian." );

    // global derivatives
    dphiDx = contraBase * dphiDxi;
    // return determinant
    return detJ;
}

//------------------------------------------------------------------------------
/** Computation of the global derivatives of second-order
 *  This function implements the identity
 *  \f[
 *    \frac{\partial^2 \phi^K}{\partial x_i \partial x_j} =
 *    \left(
 *      \frac{\partial^2 \phi^K}{\partial \xi^\alpha \partial \xi^{\beta}} -
 *      \frac{\partial^2 x_m}{\partial \xi^\alpha \partial \xi^\beta}
 *      \frac{\partial \phi^K}{\partial x_m}
 *    \right)
 *      \frac{\partial \xi^\alpha}{\partial x_i} 
 *      \frac{\partial \xi^\beta}{ \partial x_j} 
 *  \f]
 *  Summation over \f$ \alpha \f$, \f$ \beta \f$ and \f$ m \f$ is implied.
 *
 *  \param[in]  xi       Local coordinate at which the derivates are computed
 *  \param[out] ddphiDdx ddphiDdx[i,j](K) \f$= d^2 \varphi^K / d x_i d x_j \f$
 *  \return              The Jacobi determinant \f$ \sqrt{|g(\xi)|} \f$
 */
template<typename NODE, typename SFUN>
double corlib::ElementBasic<NODE,SFUN>::
globalSecondDerivatives( const VecLDim & xi, 
                         MatVecNNDimDim & ddphiDdx ) const
{
    // nodal coordinates
    MatDimNN X; this -> nodalCoordinates_( X );

    // shape function derivatives
    MatLDimNN dphiDxi; this -> sfunGrad( xi, dphiDxi );

    // shape function Hessian
    MatVecNNLDimLDim ddphiDDxi; this -> sfunHess( xi, ddphiDDxi );

    // covariant basis
    const eigenX::MatrixSd<dim,localDim> coBase 
        = X *  dphiDxi.transpose( ); 

    // contravariant basis (inverse of the covariant basis)
    eigenX::MatrixSd<dim,localDim> contraBase;
    // specialised call for determinant and contra-variant basis
    const double detJ = 
        corlib::ContraVariantBasis<dim,localDim>::compute( coBase, contraBase );

    FTL_VERIFY_DESCRIPTIVE( detJ > 0., "Negative Jacobian." );

    //--------------------------------------------------------------------------
    // compute the Hessian of the coordinate representation
    Eigen::Array<VecDim, localDim, localDim> xHess;
    {
        for ( unsigned alpha = 0; alpha < localDim; alpha ++) {
            for ( unsigned beta = 0; beta < localDim; beta ++ ) {
                xHess(alpha,beta).setZero();
            
                for ( unsigned L = 0; L < numNodes; L ++ ) {
                    xHess(alpha, beta) += 
                        X.col( L ) * ddphiDDxi( alpha, beta )( L );
                }
            }
        }

    }
    
    //--------------------------------------------------------------------------
    // compute the global first-order derivative 
    const MatDimNN dphiDx = contraBase * dphiDxi;

    //--------------------------------------------------------------------------
    // go through the cartesian directions of the partial global derivatives
    for ( unsigned i = 0; i < dim; i ++ ) { 
        for ( unsigned j = 0; j < dim; j ++ ) {

            ddphiDdx(i,j).setZero();

            // do the sums
            for ( unsigned alpha = 0; alpha < localDim; alpha ++ ) {
                for ( unsigned beta = 0; beta < localDim; beta ++ ) {

                    // for all shape functions
                    for ( unsigned K = 0; K < numNodes; K ++ ) {

                        // second order  local derivative
                        double secDeriv = ddphiDDxi( alpha, beta )( K );
                        // subtract mixed terms
                        for ( unsigned m = 0; m < dim; m ++ ) 
                            secDeriv -= xHess( alpha, beta )( m ) * dphiDx( m, K );

                        // added transformed terms to result
                        ddphiDdx( i, j )( K ) 
                            += secDeriv * contraBase( i, alpha ) * contraBase( j, beta );
                    }
                    // finished all shape functions

                } // beta
            }// alpha

        } // j
    } // i

    // return determinant
    return detJ;
}

//------------------------------------------------------------------------------
/** Evaluate the geometry description at local coordinate xi
 *  \f[
 *           x_h = \sum_k \phi^k(\xi) x^k
 *  \f]
 *  \param[in]  xi Local coordinate at which the geometry is evaluated
 *  \param[out] x  Global coordinate x
 */
template<typename NODE, typename SFUN>
void corlib::ElementBasic<NODE,SFUN>::geometry( const VecLDim & xi, 
                                                VecDim  & x ) const
{
    // nodal coordinates
    MatDimNN X; this -> nodalCoordinates_( X ); 

    // shape function evaluation
    VecNN phi;  
    this -> sfun( xi, phi );

    // compute product
    x = X * phi;
    return;
}

//------------------------------------------------------------------------------
/** Avoid dangling pointers                                                   */
template<typename NODE, typename SFUN>
void corlib::ElementBasic<NODE,SFUN>::clearNodes_( ) 
{
    for ( unsigned i = 0; i < numNodes; i ++ ) nodes_[ i ] = NULL;
    return;
}

