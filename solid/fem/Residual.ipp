// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//
//                 Computational Structural Mechanics Lab
//                         University of Cambridge
//
// This software is copyrighted and all rights are retained by the University of
// Cambridge. All use, disclosure, and/or reproduction of any part not expressly
// authorised by the University of Cambridge is prohibited. (C) 2010
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

//! @file Residual.ipp

#include <boost/functional.hpp>

//------------------------------------------------------------------------------
/** Constructor with a mesh pointer.
 *  \param[in] mesh   Pointer to the mesh
 */
template<typename MESH, typename QUAD>
solid::fem::Residual<MESH,QUAD>::Residual( MESH * mesh )
    : mesh_( mesh ), increment_( NULL ), residual_( NULL )
{
    
}

//-----------------------------------------------------------------------------
/** Query the constraints from the mesh object and fill the given container
 *  \param[in,out] constraints Vector of <dofId,value>-pairs
 */
template<typename MESH, typename QUAD>
void solid::fem::Residual<MESH,QUAD>::
fillConstraintVec( std::vector< std::pair<unsigned,double> > & constraints )
{
    typedef typename MESH::Node  Node;
    mesh_ -> iterateOverNodes( boost::bind( &Node::giveConstraints, _1, 
                                            boost::ref(constraints) ) );
    return;
}

//-----------------------------------------------------------------------------
/** Distribute the current solution to the nodes. This is usually the value of
 *  the increment corresponding to a global load step.
 *  \param[in] X  Current solution
 */
template<typename MESH, typename QUAD>
void solid::fem::Residual<MESH,QUAD>::distributeSolution( DVec_  & X ) 
{
    increment_ = &X;

    //! distribute nodal solutions to nodes
    mesh_ -> iterateOverNodes( corlib::distributorFun( this, 
                                                       &MESH::Node::copyDofArray, 
                                                       &MESH::Node::setIncrement ) );

    return;
}

//------------------------------------------------------------------------------
template<typename MESH, typename QUAD>
void solid::fem::Residual<MESH,QUAD>::giveSolution( const std::vector<unsigned> & dofIndices, 
                                                    DVec_ & result ) const
{
    for ( unsigned i = 0; i < dofIndices.size( ); i ++ ) {
        result[ i ] = (*increment_)( dofIndices[i] );
    }
    return;
}


//-----------------------------------------------------------------------------
/** Compute the residual forces at each node */
template<typename MESH, typename QUAD>
void solid::fem::Residual<MESH,QUAD>::computeNodalResiduals( ) 
{
    // clear the nodal residuals
    mesh_ -> iterateOverNodes( boost::bind( &MESH::Node::clearForce, _1 ) );

    //! compute new forces
    typedef boost::function<void( typename MESH::Element*, 
                                  const typename MESH::Element::VecLDim &,
                                  const double & ) > Integrand;
    Integrand internalForce( &MESH::Element::internalForceIntegrand );
    corlib::Integrator<QUAD,Integrand> integrator( quadrature_, internalForce );
    mesh_ -> iterateOverElements( integrator );

    return;
}

//-----------------------------------------------------------------------------
/** Collect the residual nodal forces. 
 *  \param[out] F           Vector of residual forces
 */
template<typename MESH, typename QUAD>
void solid::fem::Residual<MESH,QUAD>::collectResidual( DVec_  & F )
{
    residual_ = &F;
    residual_ -> clear();

    //! assemble the residual forces to the system matrices rhs
    mesh_ -> iterateOverNodes( corlib::vectorAssemblerFun( &MESH::Node::giveForce, 
                                                           &MESH::Node::copyDofArray, 
                                                           this ) );
 
    return;
}

//------------------------------------------------------------------------------
/** Insert values into 'RHS' which is the vector of residual forces
 *  \param[in]  force       Nodal residual force
 *  \param[in]  dofIndices  Numbers of corresponding degrees of freedom
 */
template<typename MESH, typename QUAD>
void solid::fem::Residual<MESH,QUAD>::insertToRhs( const DVec_ &   force, 
                                                   const DofVec_ & dofIndices )
{
    for ( unsigned i = 0; i < dofIndices.size( ); i ++ ) {
        (*residual_)( dofIndices[ i ] ) += force( i );
    }
    return;
}


