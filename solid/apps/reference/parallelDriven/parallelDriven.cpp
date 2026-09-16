// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//
//                 Computational Structural Mechanics Lab
//                         University of Cambridge
//
// This software is copyrighted and all rights are retained by the University of
// Cambridge. All use, disclosure, and/or reproduction of any part not expressly
// authorised by the University of Cambridge is prohibited. (C) 2010
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

//! @file parallelDriven.cpp

//------------------------------------------------------------------------------
/** \brief Parallel application of elasticity. 
 *  \detail Supports PETSc and BDDCML solvers.
 */

#include <iostream>
#include <fstream>

#include <boost/function.hpp>
#include <boost/bind.hpp>
#include <boost/lexical_cast.hpp>

#include <corlib/PropertiesParser.hpp>
#include <corlib/verify.hpp>
#include <corlib/UniqueFilename.hpp>
#include <corlib/ElementBasic.hpp>
#include <corlib/NodeBasic.hpp>

#if   defined(PETSC)
    #include "petsc.h"
    #include <corlib/SystemSolvePetsc.hpp>
#elif defined(BDDCML)
    #include <mpi.h>
    #include <corlib/SystemSolveBddc.hpp>
#endif
#include <corlib/Shapefun.hpp>
#include <corlib/DistributedNode.hpp>
#include <corlib/DistributedElement.hpp>

// solid::material::
#include <solid/material/MaterialContainer.hpp>
// solid::fem
#include <solid/fem/NodeStatic.hpp>
#include <solid/fem/ElementStatic.hpp>

#include <solid/driver/Static.hpp>

#if   defined(PETSC)
    #include <solid/driver/ParallelPetsc.hpp>
#elif defined(BDDCML)
    #include <solid/driver/ParallelBddc.hpp>
#endif

//==============================================================================
namespace apps {

    namespace ublas = boost::numeric::ublas;

    //--------------------------------------------------------------------------
    // Dirichlet boundary conditions
    corlib::NodalConstraint dirichletFun( const ublas::bounded_vector<double,2> & X )
    {
        corlib::NodalConstraint constraint;
        // wall boundaries
        if ( corlib::fuzzyEqual( X[0], 0. ) ) {  // fix bottom boundary
            constraint.setComponent( 0, 0. );
            constraint.setComponent( 1, 0. );
        }
        else if ( corlib::fuzzyEqual( X[0], 1. ) ) { // pull top boundary
            constraint.setComponent( 0, 0.5 );
            constraint.setComponent( 1, 0.0 );
        }
        return constraint;
    }
    
}

//==============================================================================
// everything starts here
int main( int argc, char *argv[] )
{

    // initialize MPI
    MPI_Init( &argc, &argv );
    MPI_Comm commAll = MPI_COMM_WORLD;

#ifdef PETSC
    // initialize PETSc
    PetscInitialize( &argc, &argv, PETSC_NULL, PETSC_NULL );
    commAll = PETSC_COMM_WORLD;
#endif

    // start timer
    const double clockStart = MPI_Wtime();

    // variables to be read from input file
    std::string basename, materialFile;
    unsigned numberOfSubSteps = 0, numberOfProcesses = 0;
    
    // feed properties parser with the variables to read
    corlib::PropertiesParser * prop = new corlib::PropertiesParser;
    prop -> registerPropertiesVar( "basename",           basename  );
    prop -> registerPropertiesVar( "materialFile",       materialFile );
    prop -> registerPropertiesVar( "numberOfSubSteps",   numberOfSubSteps );
    prop -> registerPropertiesVar( "numberOfProcesses",  numberOfProcesses );

    // read variables from the input.dat file
    std::ifstream inputFile( "./input.dat" );
    FTL_VERIFY( inputFile.is_open() );
    prop -> readValues( inputFile );
    delete prop;
    inputFile.close();
    
    // check the number of processes
    int nProc;
    MPI_Comm_size( commAll, &nProc );
    if ( numberOfProcesses > 0 )
        FTL_VERIFY_DESCRIPTIVE( static_cast<unsigned>(nProc) == numberOfProcesses,
                                "Number of MPI processors : %d, "
                                "Number of defined processors in input file : numberOfProcesses=%d\n",
                                static_cast<unsigned>(nProc), numberOfProcesses );

    // read material
    solid::material::MaterialContainer * mc = 
        solid::material::MaterialContainer::instance();
    std::ifstream mf( materialFile.c_str() );
    FTL_VERIFY( mf.is_open() );
    mc -> readMaterialStream( mf );
    mf.close();

    // some typedefs
    typedef corlib::NodeBasic<2>                                        BasisNode;
    typedef solid::fem::NodeStatic<BasisNode>                           SolidNode;
    typedef corlib::DistributedNode<SolidNode>                          Node;

    typedef corlib::Shapefun<corlib::TRIANGLE,      3>                  Sfun;
    //typedef corlib::Shapefun<corlib::QUADRILATERAL, 9>                  Sfun;
    //typedef corlib::Shapefun<corlib::HEXAHEDRON,    8>                  Sfun;
    //typedef corlib::Shapefun<corlib::TETRAHEDRON,  10>                  Sfun;

    typedef solid::material::MaterialBase                               Material;

    typedef corlib::ElementBasic<Node,Sfun>                             BasisElement;
    typedef solid::fem::ElementStatic<BasisElement,Material>            StaticElement;
    typedef corlib::DistributedElement<StaticElement>                   Element;

#if   defined(PETSC)
    typedef corlib::SystemSolvePetsc                                    SysSolve;
#elif defined(BDDCML)
    typedef corlib::SystemSolveBddc                                     SysSolve;
#endif

    // the driver --- look out for icebergs
    typedef solid::driver::Static<Element,SysSolve>                     DriverStatic;
    typedef solid::driver::Parallel<DriverStatic,SysSolve>              Driver;


    // unique filename generator
    corlib::UniqueFilename filenameGen(3);

    // read smf file and construct mesh
    int rank;
    MPI_Comm_rank( commAll, &rank );
    Driver * driverPtr = NULL;
    if ( numberOfProcesses > 0 ) {  // pre-partitioned mesh
        std::string meshFile = filenameGen( basename, rank, "psmf" );
        std::ifstream smf( meshFile.c_str() );
        FTL_VERIFY( smf.is_open() );
        driverPtr = new Driver( smf );
        smf.close();
    }
    else {                          // partition mesh on-the-fly with RCB
        std::string meshFile = basename + ".smf";
        std::ifstream smf( meshFile.c_str() );
        FTL_VERIFY( smf.is_open() );
        driverPtr = new Driver( smf, 0 );
        smf.close();
    }
    MPI_Barrier( commAll );
    FTL_VERIFY( driverPtr );
    Driver & driver = *driverPtr;

    // cache material to elements
    Material * material = mc -> getMaterial( 0 );
    driver.setMaterial( material );

    // read the constraints
    //std::string constraintsFile = filenameGen( basename, rank, "constraints" );
    //std::ifstream cstrFile( constraintsFile.c_str() );
    //FTL_VERIFY( cstrFile.is_open() );
    //driver.setNodeConstraints( cstrFile );
    //cstrFile.close();
    
    // Dirichlet boundary conditions
    driver.setNodeConstraints( boost::bind( &apps::dirichletFun, _1 ) );
    // Dirichlet BCs factor
    Driver::DoubleToDouble constrFactorFun = 
        solid::driver::detail_::ReturnConstant<double,double>( 1. / numberOfSubSteps );
    driver.setNodeConstraintsFactor( constrFactorFun );


    // container for the filenames
    driver.allocateAnimation( basename );

    // write VTU file for this process
    driver.writeAnimationStep( 0., 0 );

    // store snap shot on proc 0
    if ( rank == 0 )
        driver.writeAnimationStep( 0., 0 );

    // create linear solver
    driver.allocateSolver();
#ifdef PETSC
    // use Petsc defaults, make sure you run it with mpirun -np >=2
    driver.getSolver() -> insertOptions( " " );
#endif

    //==========================================================================
    // displacement control loop
    for ( unsigned step = 0; step < numberOfSubSteps; step ++ ) {
        if ( rank == 0 ) std::cout << "SUBSTEP: " << step << std::endl;

        //----------------------------------------------------------------------
        // Newton-Raphson iteration
        const unsigned maxiter = 10;
        const double tolerance = 1.e-5;
        driver.iterate( maxiter, corlib::Negative(),
                        boost::bind( std::less<double>(), _1, tolerance ) );

        // set increment as current displacement
        driver.updateSolution();
    
        // write VTU file for this process
        driver.writeAnimationStep( driver.getNextTime(), step+1 );

        // update displacement control
        driver.updateTime();
    }

    //--------------------------------------------------------------------------
    // write pvd file
    if ( rank == 0 ) {
        const std::string animFile = basename + ".pvd";
        driver.writeAnimationFile( animFile );
    }

    // check overall used time
    const double elapsed = MPI_Wtime() - clockStart;
    if ( rank == 0 )
        std::cout << "Elapsed time: " << corlib::secondsToHMMSS( elapsed ) << std::endl;

    // need to deallocate solver here _before_ MPI is finalised
    driver.deAllocateSolver();
    delete driverPtr;

#ifdef PETSC
    PetscFinalize();
#endif
    MPI_Finalize();

    return 0;
}

