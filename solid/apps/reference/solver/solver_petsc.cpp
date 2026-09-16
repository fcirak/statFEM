// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//
//                 Computational Structural Mechanics Lab
//                         University of Cambridge
//
// This software is copyrighted and all rights are retained by the University of
// Cambridge. All use, disclosure, and/or reproduction of any part not expressly
// authorised by the University of Cambridge is prohibited. (C) 2010
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

//! @file solver_petsc.cpp

#include <corlib/PropertiesParser.hpp>
#include <corlib/Constraints.hpp>
#include <corlib/Mesh.hpp>
#include <corlib/Shapefun.hpp>
#include <corlib/Quadrature.hpp>
#include <corlib/Integrator.hpp>
#include <corlib/MatrixAssembler.hpp>
#include <corlib/VectorAssembler.hpp>
#include <corlib/Distributor.hpp>
#include <corlib/DofNumberer.hpp>
#include <corlib/VTUwriter.hpp>
#include <corlib/SystemSolvePetsc.hpp>
#include <corlib/verify.hpp>

#include <solid/fem/NodeStatic.hpp>
#include <solid/fem/ElementStatic.hpp>

#include <solid/material/MaterialContainer.hpp>


#include <iostream>
#include <fstream>
#include <string>
#include <cstdlib>

#include "./helpers.hpp"


//------------------------------------------------------------------------------
int main( int argc, char* argv[] )
{ 
    // initialize PETSC
    PetscInitialize(&argc,&argv,0,0);

    // variables to be read from input file
    std::string meshFile, materialFile, constraintsFile;
    unsigned numberOfSubsteps = 0;
    readFromInput( meshFile, materialFile, constraintsFile, numberOfSubsteps );

    // read material
    solid::material::MaterialContainer * mc = solid::material::MaterialContainer::instance( );
    std::ifstream mf( materialFile.c_str( ) );
    if( !mf.is_open( ) ) FTL_VERIFY(!"Cannot open material file" ); 
    mc -> readMaterialStream( mf );
    mf.close( );

    //--------------------------------------------------------------------------
    // define attributes
    const corlib::shape elemShape  = corlib::QUADRILATERAL;
    const unsigned      numNodesPE = 4;
    const unsigned      numGaussP  = 4;
    typedef Attributes<elemShape,numNodesPE,numGaussP>  Attr;

    // read smf file and construct mesh
    std::ifstream smf( meshFile.c_str( ) );
    FTL_VERIFY( smf.is_open() );
    Attr::Mesh mesh( smf );
    smf.close( );

    // obtain the 1st material
    Attr::Material * material = mc -> getMaterial( 0 );

    // cache material to elements
    mesh.iterateOverElements( boost::bind( &Attr::Element::cacheMaterial, _1, material ) );

    // read the constraints
    std::ifstream cstrFile( constraintsFile.c_str( ) );
    FTL_VERIFY( cstrFile.is_open( ) );
    corlib::ConstraintsFromFile<Attr::Node>  constraintsHandler( cstrFile );
    cstrFile.close();
    mesh.iterateOverNodes( constraintsHandler );
    
    // number dofs
    const unsigned numdofs = mesh.iterateOverNodes( corlib::dofNumberFun( &Attr::Node::numberDOFs, 0 ) );
    std::cout << numdofs << " degrees of freedom" << std::endl;

    // set up an animation writer
    FTL_VERIFY( not system( "mkdir -p animation" ) );
    std::string basename = "./animation/" + meshFile;
    corlib::VTUanim animation( basename, "vtu" );
    
    // write data to VTU file
    std::string vtuFile = animation.snapshotName( 0., 0 );
    writeMeshData( vtuFile, &mesh );

    // solver
    corlib::SystemSolvePetsc * sysSolve = new corlib::SystemSolvePetsc( numdofs );

    // use Petsc defaults, make sure you run it with mpirun -np >=2
    sysSolve -> insertOptions( " " );

    // quadrature
    Attr::Quad quadrature;

    // some tolerances
    const double tolerance = 1e-12;
    const unsigned maxiter = 10;

    for ( unsigned step = 0; step < numberOfSubsteps; step ++ ) {
        std::cout << "SUBSTEP: " << step;

        for ( unsigned i = 0; i < maxiter; i ++ ) {
            std :: cout << "- Iteration number " << i;

            mesh.iterateOverNodes( boost::bind( &Attr::Node::clearForce, _1 ) );

            // internal force functor
            Attr::Integrand internalForce( &Attr::Element::internalForceIntegrand );
            Attr::Integrator integrator( quadrature, internalForce );
            mesh.iterateOverElements( integrator );

            // assemble the forces to the system matrices rhs
            mesh.iterateOverNodes( corlib::vectorAssemblerFun( &Attr::Node::giveForce, 
                                                               &Attr::Node::copyDofArray, 
                                                               sysSolve ) );


            // Compute and assemble element stiffness matrices
            corlib::MatrixComputeAndAssembleFun<const Attr::Element,Attr::Quad,corlib::SystemSolvePetsc>
                mafStiffness( &Attr::Element::stiffnessIntegrand, 
                              &Attr::Element::getDofIndices,
                              &Attr::Element::getDofIndices,
                              quadrature, sysSolve );
            mesh.iterateOverElements( mafStiffness );

            // let petsc collect itself
            sysSolve -> finishAssembly();

            // collect and apply constraints
            const double factor = ( i==0 ? 1./static_cast<double>(numberOfSubsteps) : 0. );
            corlib::ConstraintFunctor<Attr::Node,corlib::SystemSolvePetsc> constraint( sysSolve, factor );
            constraint = mesh.iterateOverNodes( constraint );
            constraint.applyConstraints();

            // norm of the residual
            const double normResidual = sysSolve -> normRhs( );
            std :: cout << "  ||residual|| = " << normResidual;
            if ( normResidual < tolerance ) break;

            // solve the system
            sysSolve -> solveSystem( );

            const double normDeltaU = sysSolve -> normRhs( );
            std :: cout << "  ||delta U|| = " << normDeltaU << std :: endl;
            if ( normDeltaU < tolerance ) break;

            // collect nodal solutions
            mesh.iterateOverNodes( corlib::distributorFun( sysSolve, 
                                                           &Attr::Node::copyDofArray, 
                                                           &Attr::Node::addToIncrement ) );

            // clear the stiffness matrix and RHS
            sysSolve -> clearMatrix( );
            sysSolve -> clearRhs( );

        }

        // add nodal increments onto nodal displacements
        mesh.iterateOverNodes( boost::bind( &Attr::Node::updateDisplacements, _1 ) );

        // write data to VTU file
        std::string vtuFile = animation.snapshotName( step+1, step+1 );
        writeMeshData( vtuFile, &mesh );

        std::cout << std::endl;
    }

    // write the animation file
    std::string animationFile( "animation.pvd" );
    std::ofstream anim( animationFile.c_str( ) );
    animation.writeAnimationFile( anim );
    anim.close();

    delete sysSolve;

    PetscFinalize();

    return 0;
}

//------------------------------------------------------------------------------
