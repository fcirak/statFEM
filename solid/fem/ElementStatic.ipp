// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//
//                 Computational Structural Mechanics Lab
//                         University of Cambridge
//
// This software is copyrighted and all rights are retained by the University of
// Cambridge. All use, disclosure, and/or reproduction of any part not expressly
// authorised by the University of Cambridge is prohibited. (C) 2010
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

//! @file   ElementStatic.ipp

#include <corlib/linalg.hpp>
#include <functional>

//------------------------------------------------------------------------------
/** Accessor to the indices of the nodal degrees of freedom
 *  \param[in] dofIndices Vector storing the dof indices
 */
template<typename BELEMENT, typename MAT>
void solid::fem::ElementStatic<BELEMENT,MAT>::
getDofIndices( std::vector<unsigned> & dofIndices ) const
{
    // resize array (important for assembler!)
    dofIndices.resize( numNodes * dof );
    std::vector<unsigned>::iterator iter = dofIndices.begin();

    typename BasisElement::NodeConstIterator first = BasisElement::nodesBegin();
    typename BasisElement::NodeConstIterator last  = BasisElement::nodesEnd();
    for ( ; first != last; ++first ) {
        std::vector<unsigned> aux;
        (*first) -> copyDofArray( aux );
        iter = std::copy( aux.begin(), aux.end(), iter );
    }
    return;
}

//------------------------------------------------------------------------------
/** Accessor to a general nodal quantitiy accessed by the pass operator
 *  \tparam     OP  Node functor which access the desired nodal quantity
 *  \param[out] Q   Matrix (dim x number of nodes) containing the result
 *  \param[in]  op  Node functor for access of nodal quantity
 */
template<typename BELEMENT, typename MAT>
template<typename OP>
void solid::fem::ElementStatic<BELEMENT,MAT>::nodalQuantity_( MatDofNN & Q, 
                                                              OP  op ) const 
{
    typename BasisElement::NodeConstIterator first = BasisElement::nodesBegin();
    typename BasisElement::NodeConstIterator last  = BasisElement::nodesEnd();
    for ( unsigned c = 0; first != last; ++first, c++ ) {
        ublas::column( Q, c ) = op( *first );
    }
    return;
}

//------------------------------------------------------------------------------
/** \param[out] U Matrix with the nodal displacements as columns              */
template<typename BELEMENT, typename MAT>
void solid::fem::ElementStatic<BELEMENT,MAT>::
nodalDisplacements( MatDofNN & U ) const
{
    this -> nodalQuantity_( U, boost::bind( &Node::giveDisplacements, _1 ) );
    return;
}

//------------------------------------------------------------------------------
/** \param[out] deltaU Matrix with the nodal increments as columns            */
template<typename BELEMENT, typename MAT>
void solid::fem::ElementStatic<BELEMENT,MAT>::
nodalIncrements( MatDofNN & deltaU ) const
{
    this -> nodalQuantity_( deltaU, boost::bind( &Node::giveIncrement, _1 ) );
    return;
}

//------------------------------------------------------------------------------
/** \param[out] deltaU Matrix with the nodal increments as columns            */
template<typename BELEMENT, typename MAT>
void solid::fem::ElementStatic<BELEMENT,MAT>::
nodalCoordinates( MatDimNN & X ) const
{
    unsigned c = 0;
    for ( typename BasisElement::NodeConstIterator first = BasisElement::nodesBegin();
          first != BasisElement::nodesEnd(); ++ first ) {
        ublas::column( X, c++ ) = (*first) -> giveCoordinates();
    }
    return;
}

//------------------------------------------------------------------------------
/** All kinematic quantities necessary for the stiffness matrix and internal 
 *  forces are computed and the Jacobi determinant is returned.
 *  \param[in]  xi       Local coordinate at which kinematics are computed
 *  \param[out] dphiDX   dphiDX[i,j] \f$= d \varphi^j / d X_i \f$
 *  \param[out] F        Deformation gradient: \f$ F[i,J] = d x_i / d X_J \f$
 *  \return              Determinant of Jacobi matrix
 */
template<typename BELEMENT, typename MAT>
double solid::fem::ElementStatic<BELEMENT,MAT>::kinematics_( const VecLDim & xi,
                                                             MatDimNN & dphiDX, 
                                                             Mat3x3 & F ) const
{
    // compute jacobian and global derivatives
    const double detG = BasisElement::globalDerivatives( xi, dphiDX );
    // nodal coordinates
    MatDimNN X; this -> nodalCoordinates( X );
    // nodal displacements
    MatDofNN U; this -> nodalDisplacements( U );
    // compute deformation gradient F = Grad (X+U)
    F = ublas::identity_matrix<double>(3); // needed if for instance dim=2,dof=2
    // Grad X : the following way is needed if dim=3,dof=3,localDim=1 (string)
    //          then Grad X is not the identity matrix
    ublas::subrange( F, 0, dim, 0, dim )  = ublas::prod( X, ublas::trans( dphiDX ) );
    // Grad U
    ublas::subrange( F, 0, dof, 0, dim ) += ublas::prod( U, ublas::trans( dphiDX ) );
    // return determinant
    return detG;
}

//------------------------------------------------------------------------------
/** Compute global derivative of deformation gradient \f$F_{iJ,K}(\xi^\alpha)\f$
 *  at given parametric co-ordinate \f$\xi^\alpha\f$
 *
 *  \param[in]  xi       Local coordinate at which the derivates are computed
 *  \param[out] fGrad    Gradient of deformation gradient
 *  \return              Jacobi determinant of deformation gradient
 */
template<typename BELEMENT, typename MAT>
double solid::fem::ElementStatic<BELEMENT,MAT>::
deformationGradientGlobalDerivatives_( const VecLDim & xi,
                                       MatDimDimDim_ & fGrad ) const
{
    // nodal coordinates
    MatDimNN X; this -> nodalCoordinates( X );
    // nodal displacements
    MatDofNN U; this -> nodalDisplacements( U );

    // shape function derivatives
    typename BasisElement::MatLDimNN dphiDxi;
    this -> BasisElement::sfunGrad( xi, dphiDxi );
    typename BasisElement::MatVecNNLDimLDim ddphiDDxi;
    this -> BasisElement::sfunHess( xi, ddphiDDxi );

    // covariant basis
    const ublas::bounded_matrix<double,dim,localDim> coG 
        = ublas::prod( X, ublas::trans( dphiDxi ) );

    // covariant metric coefficients
    ublas::bounded_matrix<double,localDim,localDim> gramG 
        = ublas::prod( ublas::trans( coG ), coG );
 
    // inverse of the metric coefficients (in place)
    const double detG = std::sqrt( corlib::inverse( gramG ) );

    // contravariant base vectors
    const ublas::bounded_matrix<double,dim,localDim> contraG
        = ublas::prod( coG, gramG );

    // first parametric derivative of displacements
    const ublas::bounded_matrix<double,dim,localDim> uGrad 
        = ublas::prod( U, ublas::trans( dphiDxi ) );

    // second parametric derivatives of coordinates and displacements
    ublas::bounded_matrix<VecDim,localDim,localDim> xHess, uHess;
    for ( unsigned a = 0; a < localDim; ++a ) {
        for ( unsigned b = 0; b < localDim; ++b ) {
            xHess( a, b ) = ublas::prod( X, ublas::trans( ddphiDDxi( a, b ) ) );
            uHess( a, b ) = ublas::prod( U, ublas::trans( ddphiDDxi( a, b ) ) );
        }
    }

    // gradient of deformation gradient F_{iJ,K} = x_{i,JK}
    for ( unsigned i = 0; i < dim; ++i ) {
        for ( unsigned j = 0; j < dim; ++j ) {
            for ( unsigned k = 0; k < dim; ++k ) {
                double fGrad_ijk = 0.;
                for ( unsigned a = 0; a < localDim; ++a ) {
                    for ( unsigned b = 0; b < localDim; ++b ) {
                        fGrad_ijk += uHess( a, b )( i ) *
                            contraG( k, b ) * contraG( j, a );
                    }
                }
                for ( unsigned a = 0; a < localDim; ++a ) {
                    for ( unsigned b = 0; b < localDim; ++b ) {
                        for ( unsigned c = 0; c < localDim; ++c ) {
                            for ( unsigned l = 0; l < dim; ++l ) {
                                fGrad_ijk -= uGrad( i, a ) * contraG( l, a ) * 
                                    contraG( k, c ) * contraG( j, b ) * xHess( b, c )( l );
                            }
                        }
                    }
                }
                fGrad( i )( j, k ) = fGrad_ijk;
            }
        }
    }

    // return determinant
    return detG;
}

//------------------------------------------------------------------------------
/** The first Piola-Kirchhoff stress tensor \b P at local coordinate xi
 *  \param[in] xi Local Coordinate
 *  \return       First Piola-Kirchhoff stress tensor
 */
template<typename BELEMENT, typename MAT>
boost::numeric::ublas::bounded_matrix<double,3,3>
solid::fem::ElementStatic<BELEMENT,MAT>::
firstPiolaKirchhoff( const VecLDim & xi ) const
{
    MatDimNN dphiDX;
    Mat3x3   F;
    this -> kinematics_( xi, dphiDX, F );

    Mat3x3   P;
    material_ -> FPKstress( F, P );
    return P;
}

//------------------------------------------------------------------------------
/** The second Piola-Kirchhoff stress tensor \b S at local coordinate xi
 *
 *  Formula: \f$S_{IJ} = (F^{-1})_{Ik} P_{kJ}\f$
 *
 *  \param[in] xi Local Coordinate
 *  \return       Second Piola-Kirchhoff stress tensor
 */
template<typename BELEMENT, typename MAT>
boost::numeric::ublas::bounded_matrix<double,3,3>
solid::fem::ElementStatic<BELEMENT,MAT>::
secondPiolaKirchhoff( const VecLDim & xi ) const
{
    MatDimNN dphiDX;
    Mat3x3   F;
    this -> kinematics_( xi, dphiDX, F );

    Mat3x3   P;
    material_ -> FPKstress( F, P );

    corlib::inverse( F );
    Mat3x3   S = ublas::prod( F, P );
    return S;
}

//------------------------------------------------------------------------------
/** The Cauchy stress tensor \f$ \sigma \f$ at local coordinate xi
 *
 *  Formula: \f$\sigma_{ij} = 1/det(F) P_{ik} (F^T)_{kJ} \f$
 *
 *  \param[in] xi   Local Coordinate
 *  \return         Cauchy stress tensor
 */
template<typename BELEMENT, typename MAT>
boost::numeric::ublas::bounded_matrix<double,3,3>
solid::fem::ElementStatic<BELEMENT,MAT>::
cauchyStress( const VecLDim & xi ) const
{
    MatDimNN dphiDX;
    Mat3x3   F;
    this -> kinematics_( xi, dphiDX, F );

    Mat3x3   P;
    material_ -> FPKstress( F, P );

    double detF = corlib::determinant( F );
    Mat3x3  sig = ublas::prod( P, ublas::trans( F ) ) / detF;
    return sig;
}

//------------------------------------------------------------------------------
/** The integral kernel for the computation of the element stiffness matrix. 
 *  Evaluates the kernel function at the local coordinate xi, multiplies the 
 *  result with the weight and stores the result in a matrix.
 *  \f[
 *    K[m d + i, n d + k] = 
 *          \int (\varphi_{,J}^m  C_{iJkL}(u)  \varphi_{,L}^n)  \det G dX
 *  \f]
 *  where \f$ d \f$ refers to the number of dofs per node.
 *  \param[in] xi     Local coordinate at which this function is evaluated
 *  \param[in] weight Corresponding quadrature weight
 */
template<typename BELEMENT, typename MAT>
void solid::fem::ElementStatic<BELEMENT,MAT>::
stiffnessIntegrand( const VecLDim & xi, const double & weight,
                    ublas::matrix<double> & result ) const
{
    // debug check
    assert( (result.size1() == numNodes * dof) and
            (result.size2() == numNodes * dof) );

    // get kinematic quantities: jacobian, Dphi/DX, F
    MatDimNN dphiDX; 
    Mat3x3   F;
    const double detG = this -> kinematics_( xi, dphiDX, F );

    // compute second derivative of strain energy ( d^2 W(F) ) / ( dF dF )
    ublas::bounded_matrix<Mat3x3, 3, 3> Cmat;
    material_ -> elasticityTensor( F, Cmat );

    // product Dphi : Cmat : Dphi
    for (unsigned m = 0; m < numNodes; m++ ) {
        for (unsigned i = 0; i < dof; i++ ) {
            for (unsigned n = 0; n < numNodes; n++ ) {
                for (unsigned k = 0; k < dof; k++ ) {
                    double kentry = 0;
                    for (unsigned J = 0; J < dof; J++ ) {
                        for (unsigned L = 0; L < dof; L++ ) {
                            kentry += Cmat( i, J )( k, L ) * 
                                dphiDX( J, m ) * dphiDX( L, n );
                        }
                    }
                    result( m * dof + i, n * dof + k ) += kentry * detG * weight;
                }
            }
        }
    }
    
    return;
}

//------------------------------------------------------------------------------
/** Computation of internal force due to stress residuals. Evaluates the kernel
 *  function at given quadrature point xi, multiplies the result with the weight
 *  and passes it to the nodes.
 *  \f[
 *       F^{int}[ m d + i ] = - \int P(i,J) \varphi_{,J}^m \det G dX
 *  \f]
 *  where \f$ d \f$ refers to the number of dofs per node.
 *  \param[in] xi      Local coordinate at which this function is evaluated
 *  \param[in] weight  Corresponding quadrature weight
 */
template<typename BELEMENT, typename MAT>
void solid::fem::ElementStatic<BELEMENT,MAT>::
internalForceIntegrand( const VecLDim & xi, const double & weight ) 
{
    // get kinematic quantities: jacobian, Dphi/DX, F
    MatDimNN dphiDX; 
    Mat3x3   F;
    const double detG = this -> kinematics_( xi, dphiDX, F );

    // compute first derivative of strain energy ( d W(F) ) / ( dF )
    Mat3x3  FPK;
    material_ -> FPKstress( F, FPK );

    // compute matrix of nodal forces 
    const MatDofNN nodalF = - ublas::prod( ublas::subrange( FPK, 0, dof, 0, dim ), 
                                           dphiDX ) * detG * weight;

    // pass each nodal force to the nodes
    typename BasisElement::NodeIterator  first = BasisElement::nodesBegin();
    typename BasisElement::NodeIterator  last  = BasisElement::nodesEnd();
    for ( unsigned n = 0; first != last; ++first, n++ ) {
        (*first) -> addToForce( ublas::column( nodalF, n ) );
    }

    return;
}

//------------------------------------------------------------------------------
/** \brief Evaluate the traction term used for the Nitsche technique at xi
 *  \details For the weak treatment of Dirichlet boundary conditions (with
 *  datum g) the required boundary integrals are
 *  \f[
 *      \int_\Gamma g \cdot t(v) ds \text{ and } \int_\Gamma u \cdot t(v) ds
 *  \f]
 *  Therefore, this routine returns the matrix with entries
 *  \f[
 *        T[i, N d + k ] = C_{ijkl} \phi^N_{,l} n_j
 *  \f]
 *  \param[in]  xi       Integration point
 *  \param[in]  normal   Normal to interface (reference configuration)
 *  \param[out] dphiDn   result matrix T
 */
template<typename BELEMENT, typename MAT>
void solid::fem::ElementStatic<BELEMENT,MAT>::
nitscheTraction( const VecLDim & xi, const VecDim & normal,
                 ublas::bounded_matrix<double,dof,numNodes*dof> & dphiDn ) const
{
    // get kinematic quantities: jacobian, Dphi/DX, F
    MatDimNN dphiDX;
    Mat3x3   F;
    this -> kinematics_( xi, dphiDX, F );

    ublas::bounded_matrix<Mat3x3, 3, 3> Cmat;
    material_ -> elasticityTensor( F, Cmat );

    for ( unsigned i = 0; i < dof; i ++ ) {
    	for ( unsigned N = 0; N < numNodes; N ++ ) {
            for ( unsigned k = 0; k < dof;  k ++ ) {
                double entry = 0.;
                for ( unsigned j = 0; j < dof; j ++ ) {
                    for ( unsigned el = 0; el < dof; el ++ ) {
                       entry += Cmat( i, j )( k, el ) * dphiDX(el, N) * normal[j];
                    }
                }
                dphiDn( i, N * dof + k ) = entry;
            }
        }
    }
    return;
}

//------------------------------------------------------------------------------
/** Computation of the internal strain energy. Evaluation of kernel function at
 *  given local coordinate xi, multiplication with a weight and addition on the
 *  energy scalar: energy =  \f$ \int W(F(\xi)) \det G dX \f$.
 *  \param[in] xi     Local coordinate at which this function is evaluated
 *  \param[in] weight  Corresponding quadrature weight
 */
template<typename BELEMENT, typename MAT>
void solid::fem::ElementStatic<BELEMENT,MAT>::
strainEnergyIntegrand( const VecLDim & xi, const double & weight ) 
{
    // get kinematic quantities: jacobian, Dphi/DX, F
    MatDimNN dphiDX; 
    Mat3x3   F;
    const double detG = this -> kinematics_( xi, dphiDX, F );

    // pass back the strain energy density due to F 
    energy_ += ( material_ -> strainEnergy( F ) ) * detG * weight;

    return;
}

//------------------------------------------------------------------------------
/** Compute nodal forces due to an applied body force. This function is given
 *  the coordinates and weight of a quadrature point. Moreover, a function object
 *  is passed with the operator of type VecDof = func(VecDim), where the 
 *  argument refers to the reference coordinate of the integration point. 
 *  \tparam FUNC       Type of function object
 *  \param[in] xi      Local coordinate of quadrature point
 *  \param[in] weight  Weight of quadrature point
 *  \param[in] func    The force function object
 *  \param[in] factor  Scalar multiplier
 */
template<typename BELEMENT, typename MAT>
template<typename FUNC>
void solid::fem::ElementStatic<BELEMENT,MAT>::
bodyForce( const VecLDim & xi, const double & weight, FUNC func, 
           const double & factor )
{
    // get global coordinate of xi
    MatDimNN nods; this -> nodalCoordinates( nods );
    VecNN phi;     BasisElement::sfun( xi, phi );
    const VecDim X = ublas::prod( nods, phi );

    // compute value of body force function
    const VecDof fofX = func( X );

    // get jacobian of this element
    const double detG = BasisElement::jacobian( xi );

    // nodal forces
    const MatDofNN F = detG * weight * factor * ublas::outer_prod( fofX, phi );

    // Compute the nodal forces 
    typename BasisElement::NodeIterator  first = BasisElement::nodesBegin();
    typename BasisElement::NodeIterator  last  = BasisElement::nodesEnd();
    for ( unsigned n = 0; first != last; ++first, n++ ) {
        (*first) -> addToForce( ublas::column( F, n ) );
    }
    return;
    
}

//------------------------------------------------------------------------------
/** Evaluate stress residual at given parametric point
 *
 *  Stress residual is
 *  \f[
 *     r_i(\xi^\alpha) = P_{iJ,J}(\xi^\alpha) + f_i(\xi^\alpha)
 *  \f]
 *  for \f$\xi^\alpha \in (0,1)^{dim} \f$.
 *
 * \tparam       FUNC       Function type: X_J -> f_i
 *
 * \param[in]    xi         Parametric evaluation point
 * \param[in]    func       Body load density function
 * \param[out]   stressRes  Stress residual
 */
template<typename BELEMENT, typename MAT>
template<typename FUNC>
void solid::fem::ElementStatic<BELEMENT,MAT>::stressResidual( const VecLDim & xi,
                                                              FUNC func,
                                                              VecDim & stressRes )
{
    // get global coordinate of xi
    MatDimNN nods; this -> nodalCoordinates( nods );
    VecNN phi;     BasisElement::sfun( xi, phi );
    const VecDim X = ublas::prod( nods, phi );

    // compute value of body force function
    const VecDof fofX = func( X );

    // compute gradient of def.grad.
    MatDimDimDim_ fGrad;
    this -> deformationGradientGlobalDerivatives_( xi, fGrad );

    // identity matrix
    const Mat3x3 id3 = ublas::identity_matrix< double >( 3 );

    // compute stress resultant vector
    for ( unsigned d = 0; d < dim; ++d ) {
        Mat3x3 fGrad3;  fGrad3.clear();
        ublas::subrange( fGrad3, 0,dim, 0,dim ) = fGrad( d );
        Mat3x3 pGrad3;  
        material_ -> FPKstress( fGrad3, pGrad3 );  // this is incomplete for non-linear material
        stressRes[ d ] = corlib::inner_prod( pGrad3, id3 ) + fofX[ d ];
    }

    return;
}
