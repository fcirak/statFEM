// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//
//                 Computational Structural Mechanics Lab
//                         University of Cambridge
//
// This software is copyrighted and all rights are retained by the University of
// Cambridge. All use, disclosure, and/or reproduction of any part not expressly
// authorised by the University of Cambridge is prohibited. (C) 2010
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

//! @file fem.hpp

#ifndef del2_apps_reference_laplaceInverse_fem_h
#define del2_apps_reference_laplaceInverse_fem_h

//！boost lib headers
#include <boost/function.hpp>

namespace eigenX = corlib::eigenX;

//------------------------------------------------------------------------------
template<unsigned DIM>
eigenX::VectorSd<1>
force( const eigenX::VectorSd<DIM> x )
{
    //! analytical solution: quadratic function
    eigenX::VectorSd<1> f;
    f[0] = 1. ;
    return f;
}

//------------------------------------------------------------------------------
// homogeneous Dirichlet boundary conditions
template<unsigned DIM>
corlib::NodalConstraint
dirichletFun( const eigenX::VectorSd<DIM> x )
{
    corlib::NodalConstraint constraint;
    if ( corlib::fuzzyEqual( x[0], 0. ) or
         corlib::fuzzyEqual( x[0], 1. ) ) {
        constraint.setComponent( 0, 0. );
    }
    return constraint;
}


//------------------------------------------------------------------------------
template< typename MESH, typename QUAD, typename SYSSOLVE >
void fem( MESH & mesh, QUAD quadrature, SYSSOLVE * sysmat )
{

    typedef MESH                      Mesh;
    typedef QUAD                      Quad;
    typedef SYSSOLVE                  SysSolve;
    typedef typename Mesh::Node       Node;
    typedef typename Mesh::Element    Element;
    typedef typename Node::VecDim     VecDim;

    constexpr auto dim = Node::dim;
    constexpr auto dof = Node::dof;

    //! Apply dirichlet boundary conditions
    typedef boost::function<corlib::NodalConstraint( const VecDim )>     VecDim2Constraint;
    typedef corlib::ConstraintsFromFunction<Node,VecDim2Constraint>      Dirichlet;
    Dirichlet dirichlet( boost::bind( dirichletFun<dim>, _1 ) );
    mesh.iterateOverNodes( dirichlet );

    //! compute element stiffness matrices
    corlib::MatrixComputeAndAssembleFun<const Element, Quad, SysSolve>
        maf( &Element::stiffnessIntegrand,
             &Element::getDofIndices,
             &Element::getDofIndices,
             quadrature, sysmat );
    mesh.iterateOverElements( maf );

    //! Clear force arrays
    mesh.iterateOverNodes( boost::bind( &Node::clearForce, _1 ) );

    //! body force
    typedef eigenX::VectorSd<1>                                         Vec1;
    typedef boost::function<Vec1( const VecDim )>                       VecDim2Vec1;
    typedef corlib::BodyForce<Element,VecDim2Vec1>                      BodyForce;
    BodyForce bodyForce( boost::bind( force<dim>, _1 ) );
    corlib::Integrator<Quad,BodyForce> bodyForceIntegrator( quadrature, bodyForce );
    mesh.iterateOverElements( bodyForceIntegrator );

    //! assemble the  forces to the system matrices rhs
    mesh.iterateOverNodes( corlib::vectorAssemblerFun( &Node::getForce,
                                                       &Node::copyDofArray,
                                                       sysmat ) );

    //! collect and apply constraints
    sysmat -> finishAssembly( );
    corlib::ConstraintFunctor< Node, SysSolve > constraint( sysmat );
    constraint = mesh.iterateOverNodes( constraint );
    constraint.applyConstraints( );

}

#endif
