// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//
//                 Computational Structural Mechanics Lab
//                         University of Cambridge
//
// This software is copyrighted and all rights are retained by the University of
// Cambridge. All use, disclosure, and/or reproduction of any part not expressly
// authorised by the University of Cambridge is prohibited. (C) 2010
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

//! @file staticDrivenFine.cpp

#include <iostream>
#include <fstream>
#include <string>
#include <cstdlib>

#include <boost/function.hpp>

// corlib::
#include <corlib/Shapefun.hpp>
#include <corlib/ElementBasic.hpp>
#include <corlib/NodeBasic.hpp>
#include <corlib/PropertiesParser.hpp>
#include <corlib/verify.hpp>

// decide upon solver
#if   defined PETSC
    #include <corlib/SystemSolvePetsc.hpp>
#elif defined SUPERLU
    #include <corlib/SystemSolveLU.hpp>
#elif defined PARDISO_MKL
    #include <corlib/SystemSolvePardiso.hpp>
#else
    #include <corlib/SystemSolve.hpp>
#endif
// solid::material
#include <solid/material/MaterialContainer.hpp>
#include <solid/material/PlaneStressDecorator.hpp>
// solid::fem
#include <solid/fem/NodeStatic.hpp>
#include <solid/fem/ElementStatic.hpp>

#include <solid/driver/Static.hpp>
#include <solid/fem/NodalForces.hpp>


//==============================================================================
namespace apps {

    namespace ublas = boost::numeric::ublas;

    //--------------------------------------------------------------------------
    // Dirichlet boundary condition
    template<unsigned DIM>
    corlib::NodalConstraint dirichletFun( const ublas::bounded_vector<double,DIM> X )
    {
        corlib::NodalConstraint constraint;
        // wall boundaries
        if ( corlib::fuzzyEqual( X[1], 0. ) ) {  // fix bottom boundary
            constraint.setComponent( 0, 0. );
            constraint.setComponent( 1, 0. );
            if (DIM == 3 ) constraint.setComponent( 2, 0. );
        }
        else if ( corlib::fuzzyEqual( X[1], 1. ) ) { // pull top boundary
            constraint.setComponent( 1, 1.0 );
        }
        return constraint;
    }

    //--------------------------------------------------------------------------
    // Body force function (has to be given the coordinate as argument)
    template<unsigned DIM>
    ublas::bounded_vector<double,DIM>
    gravity( const ublas::bounded_vector<double,DIM> X )
    {
        boost::numeric::ublas::bounded_vector<double,DIM> G;
        G.clear();
        G[DIM-1] = -100.;
        return G;
    }

    //--------------------------------------------------------------------------
    template<unsigned DIM>
    bool nodalForceFun( const ublas::bounded_vector<double,DIM> & X,
                        ublas::bounded_vector<double,DIM> & F )
    {
        bool applyAforce = false;
        if ( corlib::fuzzyEqual( X[0], 1. ) and corlib::fuzzyEqual( X[1], 1. ) ) {
            F.clear();
            applyAforce = true;
            F[0] = 30.;
        }
        return applyAforce;
    }

}

//==============================================================================
// everything starts here
int main( int argc, char* argv[] )
{
#ifdef PETSC
    // initialize PETSC
    PetscInitialize( &argc, &argv, PETSC_NULL, PETSC_NULL );
#endif

    // Basic (variable) attributes
    const corlib::shape elemShape  = corlib::QUADRILATERAL;
    const unsigned      numNodesPE = 9;

    // Basic (derived) attributes
    const unsigned dim = corlib::ShapeTraits<elemShape>::dim;
    typedef corlib::Shapefun<elemShape,numNodesPE>            Sfun;
    typedef corlib::NodeBasic<dim>                            BasisNode;
    typedef solid::fem::NodeStatic<BasisNode>                 Node;
    typedef solid::material::MaterialBase                     Material;
    typedef corlib::ElementBasic<Node,Sfun>                   BasisElement;
    typedef solid::fem::ElementStatic<BasisElement,Material>  Element;
    typedef corlib::Mesh<Element>                             Mesh;

    // variables to be read from input file
    std::string meshFile, materialFile, constraintsFile, nodalForcesFile;
    unsigned numberOfSubSteps = 0, numberOfSubIterations = 0;

    // feed properties parser with the variables to read
    corlib::PropertiesParser * prop = new corlib::PropertiesParser;
    prop -> registerPropertiesVar( "meshFile",              meshFile   );
    prop -> registerPropertiesVar( "materialFile",          materialFile );
    prop -> registerPropertiesVar( "constraintsFile",       constraintsFile );
    prop -> registerPropertiesVar( "numberOfSubSteps",      numberOfSubSteps );
    prop -> registerPropertiesVar( "nodalForcesFile",       nodalForcesFile );
    prop -> registerPropertiesVar( "numberOfSubIterations", numberOfSubIterations );

    // read variables from the input.dat file
    const std::string inputFile = "./input.dat";
    std::ifstream inp( inputFile.c_str()  );
    FTL_VERIFY_DESCRIPTIVE( inp.is_open(),
                            "Cannot open input file %s\n", inputFile.c_str() );
    prop -> readValues( inp );
    delete prop;
    inp.close( );

    // read material
    solid::material::MaterialContainer * mc = solid::material::MaterialContainer::instance( );
    std::ifstream mf( materialFile.c_str( ) );
    FTL_VERIFY( mf.is_open() );
    mc -> readMaterialStream( mf );
    mf.close( );

    // the solver
#ifdef PETSC
    typedef corlib::SystemSolvePetsc                             SysSolve;
#elif defined SUPERLU
    typedef corlib::SystemSolveLU                                SysSolve;
#elif defined PARDISO_MKL
    typedef corlib::SystemSolvePardiso                           SysSolve;
#else
    typedef corlib::SystemSolve                                  SysSolve;
#endif

    // the driver --- travel safely
    typedef solid::driver::Static<Element,SysSolve>              Driver;

    typedef solid::fem::NodalForces<Node>                        NodalForces;

    // read smf file and construct mesh
    std::ifstream smf( meshFile.c_str( ) );
    FTL_VERIFY( smf.is_open() );
    Driver driver( smf );
    smf.close( );

    // cache 1st material to elements
    Material * material = mc -> getMaterial( 0 );
    driver.setMaterial( material );

    // set Dirichlet BCs
    driver.setNodeConstraints( boost::bind( &apps::dirichletFun<dim>, _1 ) );

    // set up an animation writer
    FTL_VERIFY( not system( "mkdir -p animation" ) );
    const std::string basename = "./animation/" + meshFile;
    corlib::VTUanim animation( basename, "vtu" );

    // write data to VTU file
    const std::string vtuFile = animation.snapshotName( 0., 0 );
    driver.writeMeshData( vtuFile );

    // nodal forces from file
    std::ifstream nf( nodalForcesFile.c_str() );
    FTL_VERIFY( nf.is_open() );
    driver.addNodalLoads( nf );
    nf.close();

    // nodal forces from function
    driver.addNodalLoads( boost::bind( &apps::nodalForceFun<dim>, _1, _2 ) );

    // number DOFs
    const unsigned numdofs = driver.numberDofs();
    std::cout << numdofs << " degrees of freedom" << std :: endl;

    // set up solver
    driver.allocateSolver( numdofs );
#ifdef PETSC
    // use Petsc defaults, make sure you it with mpirun -np >=2
    driver.getSolver() -> insertOptions( " " );
#endif
 
    //==========================================================================
    // displacement control loop
    for ( unsigned step = 0; step < numberOfSubSteps; step ++ ) {
        std::cout << "-----------------------------------------" << std::endl
                  << "Substep " << step +1 << " of " << numberOfSubSteps << std::endl;

        //----------------------------------------------------------------------
        // Newton-Raphson iteration
        const double tolerance = 1.e-5;

        for ( unsigned i = 0; i < numberOfSubIterations; i ++ ) {
            std :: cout << "  Iteration number " << i << std :: endl;

            // blank right hand side of linearised system
            driver.clearForces();

            // apply nodal forces
            const double loadFactor = double(step+1)/double(numberOfSubSteps);
            driver.computeNodalLoads( loadFactor );

            // apply body forces
            driver.computeBodyLoad( boost::bind( &apps::gravity<dim>, _1 ), loadFactor );

            // compute nodal  forces
            driver.computeInternalForces();

            // assemble the  forces to the system matrices rhs
            driver.assembleForces();

            // assemble the element stiffness matrices to the system matrix
            driver.computeAndAssembleStiffness();

            // process the matrix and right hand side
            driver.getSolver() -> finishAssembly();

            // collect and apply constraints (depending on load step)
            const double factor = ( i==0 ? 1./double(numberOfSubSteps) : 0. );
            std::cout << "    Constraint factor " << factor << std::endl;
            driver.applyConstraints( factor );

            // norm of the residual
            const double normResidual = driver . normRhs();
            std::cout << "    Norm of residual  " << normResidual << std :: endl;
            if ( normResidual < tolerance ) {
                break;
            }

            // solve the system
            driver.solveSystem();

            // add nodal solutions (ie residual displacements) to incremental displ
            driver.distributeSolution();

            // check norm of residual displacements
            const double normDeltaU = driver . normSolution();
            std::cout << "    Norm of delta U  " << normDeltaU << std :: endl;
            if ( normDeltaU < tolerance ) {
                break;
            }

            // clear the stiffness matrix and RHS
            driver.clearSystem();

        }

        // update displacements
        driver.updateSolution();

        // write data to VTU file
        const std::string vtuFile = animation.snapshotName( step+1, step+1 );
        driver.writeMeshData( vtuFile );

    }

    //------------------------------------------------------------
    // write animation file
    const std::string animationFile( "animation.pvd" );
    std::ofstream anim( animationFile.c_str( ) );
    animation.writeAnimationFile( anim );
    anim.close();

    // clean-up
    mc -> destroy();
#ifdef PETSC
    driver.deAllocateSolver();
    PetscFinalize();
#endif
    return 0;
}

