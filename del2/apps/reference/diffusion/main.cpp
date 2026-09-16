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

#include <iostream>
#include <fstream>
#include <string>
#include <cstdlib>
#include <algorithm>

#include <boost/function.hpp>

#include <Eigen/Core>
#include <corlib/eigenX.hpp>

#include <corlib/Mesh.hpp>
#include <corlib/Shapefun.hpp>
#include <corlib/DofNumberer.hpp>
#include <corlib/Quadrature.hpp>
#include <corlib/SystemSolveEigenSparse.hpp>
#include <corlib/Distributor.hpp>
#include <corlib/Accessor.hpp>
#include <corlib/MatrixAssembler.hpp>
#include <corlib/VectorAssembler.hpp>
#include <corlib/Constraints.hpp>
#include <corlib/PropertiesParser.hpp>
#include <corlib/VTUwriter.hpp>
#include <corlib/linalg.hpp>
#include <corlib/fuzzyEqual.hpp>
#include <corlib/verify.hpp>
#include <corlib/NodeBasic.hpp>
#include <corlib/ElementBasic.hpp>
#include <corlib/Integrator.hpp>

#include <del2/fem/NodePotential.hpp>
#include <del2/fem/ElementDiffusion.hpp>

namespace eigenX = corlib::eigenX;

//------------------------------------------------------------------------------
// inhomogeneous Dirichle boundary conditions
template<unsigned DIM>
corlib::NodalConstraint 
dirichletFun( const eigenX::VectorSd<DIM> & X )
{
    // note the reduced precision in the fuzzyEqual call
    corlib::NodalConstraint constraint;
    if ( corlib::fuzzyEqual( X[0]*X[0]+X[1]*X[1], 0.25, 1e-6 ) )
        constraint.setComponent( 0, 1. );
    
    return constraint;
}

//------------------------------------------------------------------------------
template<typename MESH>
void writeMeshData( const std::string & vtuFile, MESH * mesh );


//------------------------------------------------------------------------------
int main( int argc, char* argv[] )
{ 
    // variables to be read from input file
    std::string meshFile;
    double stepSize; 
    unsigned numSteps;
    double rho, kappa;

    // feed properties parser with the variables to read
    corlib::PropertiesParser *prop = new corlib::PropertiesParser;
    prop -> registerPropertiesVar( "stepSize",         stepSize   );
    prop -> registerPropertiesVar( "numSteps",         numSteps   );
    prop -> registerPropertiesVar( "meshFile",         meshFile   );
    prop -> registerPropertiesVar( "rho",              rho        );
    prop -> registerPropertiesVar( "kappa",            kappa      );

    // read variables from the input.dat file
    std::ifstream inputFile( "./input.dat" );
    FTL_VERIFY( inputFile.is_open() );
    prop -> readValues( inputFile );
    delete prop;
    inputFile.close( );

    //! problem dimension
    const unsigned dim = 2;

    //! set some typedefs
    typedef corlib::Quadrature<corlib::TRIANGLE, 3>             Quad;
    typedef corlib::Shapefun<  corlib::TRIANGLE, 6>             Sfun;
    typedef corlib::NodeBasic<dim>                              BasisNode;
    typedef del2::fem::NodePotential<BasisNode>                 Node;
    typedef corlib::ElementBasic<Node,Sfun>                     BasisElement;
    typedef del2::fem::ElementDiffusion<BasisElement>           Element;
    typedef corlib::Mesh<Element>                               Mesh;
    typedef boost::function<void( Element*,
                                  const Node::VecDim &, 
                                  const double & ) >            Integrand;
    typedef corlib::SystemSolveEigenSparse                      SysSolve;

    // read smf file and construct mesh
    std::ifstream smf( meshFile.c_str( ) );
    FTL_VERIFY( smf.is_open() );
    Mesh mesh( smf );
    smf.close( );

    //! cache material to elements
    mesh.iterateOverElements( std::bind( &Element::setConductivity,
                                         std::placeholders::_1, kappa ));
    mesh.iterateOverElements( std::bind( &Element::setDensity,
                                         std::placeholders::_1, rho   ));

    //! Apply dirichlet boundary conditions 
    typedef eigenX::VectorSd<dim>                                        VecDim;
    typedef boost::function< corlib::NodalConstraint( const VecDim)>     VecDimToConstraint;
    typedef corlib::ConstraintsFromFunction<Node,VecDimToConstraint>     Dirichlet;
    Dirichlet dirichlet( std::bind( dirichletFun<dim>, std::placeholders::_1 ) );
    mesh.iterateOverNodes( dirichlet );

    //! number dofs
    unsigned numdofs = mesh.iterateOverNodes( corlib::dofNumberFun( &Node::numberDOFs, 0 ) );
    std :: cout << numdofs << " degrees of freedom" << std :: endl;

    //! system matrix object 
    SysSolve * sysmat  = new SysSolve( numdofs );

    //--------------------------------------------------------------------------
    // set up an animation writer
    int check = system( "mkdir -p animation" ); FTL_VERIFY( !check );
    std::string basename = "./animation/diffuse";
    corlib::VTUanim animation( basename, "vtu" );

    // write data to VTU file
    std :: string vtuFile = animation.snapshotName( 0., 0 );
    writeMeshData( vtuFile, &mesh );

    //! construct a quadrature
    Quad quadrature;

    //--------------------------------------------------------------------------
    // Time loop
    for ( unsigned step = 0; step < numSteps; step ++ ) {

        //! Compute inertial forces
        mesh.iterateOverNodes(    std::bind( &Node::clearForce, std::placeholders::_1 ) );
        Integrand inertialForceIntegrand = 
            std::bind( &Element::computeInertialForces, std::placeholders::_1,
                       std::placeholders::_2, std::placeholders::_3, stepSize );
        corlib::Integrator<Quad,Integrand> inertialForceIntegrator( quadrature, 
                                                                    inertialForceIntegrand );
        mesh.iterateOverElements( inertialForceIntegrator );

        //! compute element stiffness matrices   
        corlib::MatrixComputeAndAssembleFun<const Element,Quad,SysSolve>
            mafStiffness( &Element::stiffnessIntegrand, 
                          &Element::getDofIndices,
                          &Element::getDofIndices,
                          quadrature, sysmat );
        mesh.iterateOverElements( mafStiffness );

        //! compute element mass matrices   
        corlib::MatrixComputeAndAssembleFun<const Element,Quad,SysSolve>
            mafMass( &Element::massIntegrand, 
                     &Element::getDofIndices,
                     &Element::getDofIndices,
                     quadrature, sysmat, 1./stepSize );
        mesh.iterateOverElements( mafMass );


        //! assemble the  forces to the system matrices rhs
        mesh.iterateOverNodes( corlib::vectorAssemblerFun( &Node::getForce, 
                                                           &Node::copyDofArray,
                                                           sysmat ) );
    
        //! collect and apply constraints
        sysmat -> finishAssembly();
        corlib::ConstraintFunctor< Node, SysSolve > constraint( sysmat );
        constraint = mesh.iterateOverNodes( constraint );
        constraint.applyConstraints();

        sysmat -> solveSystem( );

        //! add nodal solutions to increment 
        mesh.iterateOverNodes( corlib::distributorFun( sysmat, 
                                                       &Node::copyDofArray, 
                                                       &Node::setIncrement ) );

        //! clear the stiffness matrix and RHS
        sysmat -> clearMatrix( );
        sysmat -> clearRhs( );

        mesh.iterateOverNodes( std::bind( &Node::updatePotential, std::placeholders::_1 ) );

        // write data to VTU file
        std :: string vtuFile = animation.snapshotName( (step+1)*stepSize, step+1 );
        writeMeshData( vtuFile, &mesh );
        
    }

    // write the animation file
    std::string animationFile( "animation.pvd" );
    std::ofstream anim( animationFile.c_str( ) );
    animation.writeAnimationFile( anim );
    anim.close();

    delete sysmat;

    return 0;
}

//------------------------------------------------------------------------------
template< typename MESH>
void writeMeshData( const std::string & vtuFile, MESH * mesh )
{
    // set up the VTU writer and open a file
    std::ofstream vtu( vtuFile.c_str( ) );
    corlib::VTUwriter<MESH> vtuwriter( mesh, vtu );

    //! write mesh data
    vtuwriter.writeMesh();
    
    //-- write point data
    vtuwriter.openPointData();
    vtuwriter.writeNodalQuantity( corlib::accessorFun( &MESH::Node::getPotential, "Potential" ) );
    vtuwriter.closePointData( );

    // finish writing
    vtuwriter.finalize();
    vtu.close( );
}
