// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//
//                 Computational Structural Mechanics Lab
//                         University of Cambridge
//
// This software is copyrighted and all rights are retained by the University of
// Cambridge. All use, disclosure, and/or reproduction of any part not expressly
// authorised by the University of Cambridge is prohibited. (C) 2010
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

//! @file main.cpp

//------------------------------------------------------------------------------
#include <iostream>
#include <fstream>
#include <string>
#include <cstdlib>

#include <Eigen/Core>
#include <corlib/eigenX.hpp>

#include <corlib/Mesh.hpp>
#include <corlib/Shapefun.hpp>
#include <corlib/Quadrature.hpp>
#include <corlib/SystemSolveEigenSparse.hpp>
#include <corlib/fuzzyEqual.hpp>
#include <corlib/verify.hpp>
#include <corlib/NodeBasic.hpp>
#include <corlib/ElementBasic.hpp>

#include <del2/fem/NodePotential.hpp>
#include <del2/fem/ElementPotential.hpp>
#include <del2/fem/DriverPotential.hpp>

namespace eigenX = corlib::eigenX;

//------------------------------------------------------------------------------
eigenX::VectorSd<1>
analytical( const eigenX::VectorSd<2> & X )
{
    eigenX::VectorSd<1> u;
    u[0] = std::sin( 6. * X[0] ) * std::sin( 8. * X[1] );
    return u;
}

//------------------------------------------------------------------------------
// Body force function (has to be given the coordinate as argument)
eigenX::VectorSd<1>
force( const eigenX::VectorSd<2> & X )
{
    eigenX::VectorSd<1> U = analytical(X);
    return 100. * U;
}

// homogeneous Dirichle boundary conditions
corlib::NodalConstraint dirichletFun( const eigenX::VectorSd<2> & X )
{
    corlib::NodalConstraint constraint;
    if ( corlib::fuzzyEqual( X[0], 0. ) or
         corlib::fuzzyEqual( X[1], 0. ) or
         corlib::fuzzyEqual( X[0], 1. ) or
         corlib::fuzzyEqual( X[1], 1. ) ) {
        
        eigenX::VectorSd<1> U = analytical(X);
        constraint.setComponent( 0, U[0] );
    }
    return constraint;
}

//==============================================================================
// everything starts here
int main( int argc, char* argv[] )
{ 

    //! set some typedefs
    typedef corlib::Shapefun< corlib::QUADRILATERAL, 9 >              SfunType;

    typedef corlib::NodeBasic<2>                                      BasisNode;
    typedef del2::fem::NodePotential<BasisNode>                       NodeType;

    typedef corlib::ElementBasic<NodeType, SfunType>                  BasisElement;
    typedef del2::fem::ElementPotential< BasisElement >               ElementType;
    typedef corlib::Mesh< ElementType >                               MeshType;

    // the global linear system and its linear solver
    typedef corlib::SystemSolveEigenSparse                            SysSolve;

    // the driver
    typedef del2::fem::DriverPotential< MeshType, SysSolve >       DriverPotential;

    // create driver
    DriverPotential driver( "input.dat" );

    // Dirichlet constraints
    driver.createDirichletConstraints( std::bind( dirichletFun, std::placeholders::_1 ) );

    // body forces
    driver.createBodyForce( std::bind( force, std::placeholders::_1 ) );

    // number DOFs
    const unsigned numdofs = driver.numberDofs( 0 );
    std::cout << numdofs << " degrees of freedom" << std :: endl;

    // create system solver and initialise it
    driver.createSystem( );
    driver.clearSystem( );  // only needed in case of several steps

    // assemble residual and Jacobian
    driver.assembleResidual( );
    driver.assembleJacobian( );
    driver.assembleFinish( );

    // apply constraints due Dirichlet BCs
    driver.dirichletConstrainSystem( );

    // solve global system
    driver.solveSystem( );

    // update iteration and step
    driver.updateIncrement( );
    driver.updateStep( );

    // write result to file
    driver.writeVtu( "laplace.vtu" );

    // destroy system
    driver.destroySystem( );

    // finish
    return 0;
}
