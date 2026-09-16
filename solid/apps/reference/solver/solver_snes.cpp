// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//
//                 Computational Structural Mechanics Lab
//                         University of Cambridge
//
// This software is copyrighted and all rights are retained by the University of
// Cambridge. All use, disclosure, and/or reproduction of any part not expressly
// authorised by the University of Cambridge is prohibited. (C) 2010
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

//! @file solver_snes.cpp

#include <corlib/PropertiesParser.hpp>
#include <corlib/Mesh.hpp>
#include <corlib/Shapefun.hpp>
#include <corlib/Quadrature.hpp>
#include <corlib/Constraints.hpp>
#include <corlib/DofNumberer.hpp>
#include <corlib/VTUwriter.hpp>
#include <corlib/PetscSNESwrapper.hpp>
#include <corlib/verify.hpp>

#include <solid/fem/NodeStatic.hpp>
#include <solid/fem/ElementStatic.hpp>
#include <solid/fem/Residual.hpp>

#include <solid/material/MaterialContainer.hpp>

#include <iostream>
#include <fstream>
#include <string>
#include <cstdlib>

#include <petscksp.h>
#include <petscsnes.h>

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

    // residual object for callback to SNES solver
    typedef solid::fem::Residual<Attr::Mesh, Attr::Quad> Residual;
    Residual residual(&mesh);

    // Petsc's SNES without matrices: call with -snes_mf
    std::cout << "Using the SNES module of PETSC without matrix (-snes_mf option)"
              << std::endl;
    corlib::PetscSNESwrapper<Residual> solver( &residual, numdofs, 1./numberOfSubsteps );

    // initialize, set initial guess etc.
    solver.initialize();

    for ( unsigned step = 0; step < numberOfSubsteps; step ++ ) {
        std::cout << "SUBSTEP: " << step << std::endl;
        
        solver.solve();

        // copy increment into displacements
        mesh.iterateOverNodes( boost::bind( &Attr::Node::updateDisplacements, _1 ) );
      
        // write data to VTU file
        std :: string vtuFile = animation.snapshotName( step+1, step+1 );
        writeMeshData( vtuFile, &mesh );
    }

    solver.clear();

    // write the animation file
    std::string animationFile( "animation.pvd" );
    std::ofstream anim( animationFile.c_str( ) );
    animation.writeAnimationFile( anim );
    anim.close();

    PetscFinalize();

    return 0;
}

//------------------------------------------------------------------------------
