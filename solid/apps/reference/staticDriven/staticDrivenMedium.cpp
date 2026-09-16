// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//
//                 Computational Structural Mechanics Lab
//                         University of Cambridge
//
// This software is copyrighted and all rights are retained by the University of
// Cambridge. All use, disclosure, and/or reproduction of any part not expressly
// authorised by the University of Cambridge is prohibited. (C) 2010
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

//! @file staticDrivenMedium.cpp

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
    PetscInitialize( &argc, &argv, 0, 0 );
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
    typedef corlib::Mesh<Element>                              Mesh;

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
    // Dirichlet BCs factor
    Driver::DoubleToDouble constrFactorFun = 
        solid::driver::detail_::ReturnConstant<double,double>( 1. / numberOfSubSteps );
    driver.setNodeConstraintsFactor( constrFactorFun );

    // nodal forces from file
    std::ifstream nf( nodalForcesFile.c_str() );
    FTL_VERIFY( nf.is_open() );
    driver.addNodalLoads( nf );
    nf.close();

    // nodal forces from function
    driver.addNodalLoads( boost::bind( &apps::nodalForceFun<dim>, _1, _2 ) );

    // set body load density function
    driver.setBodyLoad( boost::bind( &apps::gravity<dim>, _1 ) );

    // load factor function
    Driver::DoubleToDouble loadFactorFun = 
        boost::bind( std::divides<double>(), _1, numberOfSubSteps );
    driver.setNodalLoadsFactor( loadFactorFun );
    driver.setBodyLoadFactor( loadFactorFun );

    // set up an animation writer
    FTL_VERIFY( not system( "mkdir -p animation" ) );
    const std::string basename = "./animation/" + meshFile;
    driver.allocateAnimation( basename );

    // write data to VTU file
    driver.writeAnimationStep( 0., 0 );

#ifdef PETSC
    // create solver
    driver.allocateSolver();
    // use Petsc defaults, make sure you run it with mpirun -np >=2
    driver.getSolver() -> insertOptions( " " );
#endif

    //==========================================================================
    // displacement and load control loop
    for ( unsigned step = 0; step < numberOfSubSteps; step ++ ) {
        std::cout << "-----------------------------------------" << std::endl
                  << "Substep " << step +1 << " of " << numberOfSubSteps << std::endl;

        //----------------------------------------------------------------------
        // Newton-Raphson iteration
        driver.iterate( numberOfSubIterations,
                        boost::bind( std::less<double>(), _1, 1.e-5 ),
                        boost::bind( std::less<double>(), _1, 1.e-5 ) );

        // update displacements
        driver.updateSolution();

        // write data to VTU file
        driver.writeAnimationStep( driver.getNextTime(), step+1 );

        // update control 
        driver.updateTime( );
    }

    //------------------------------------------------------------
    // write animation file
    driver.writeAnimationFile( "animation.pvd" );

    // clean-up
    mc -> destroy();
#ifdef PETSC
    driver.deAllocateSolver();
    PetscFinalize();
#endif
    return 0;
}

