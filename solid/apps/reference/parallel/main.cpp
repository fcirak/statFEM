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
/** \brief Parallel application of elasticity. 
 *  \detail Supports PETSc and BDDCML solvers.
 */

#include <iostream>
#include <fstream>

#include <boost/function.hpp>
#include <boost/bind.hpp>

#include <corlib/PropertiesParser.hpp>
#include <corlib/verify.hpp>
#include <corlib/VTUwriter.hpp>
#include <corlib/UniqueFilename.hpp>
#include <corlib/Mesh.hpp>
#include <corlib/Shapefun.hpp>
#include <corlib/Quadrature.hpp>
#include <corlib/Integrator.hpp>
#include <corlib/GenerateDofStarts.hpp>
#include <corlib/DofNumbererFromStarts.hpp>
#include <corlib/ElementBasic.hpp>
#include <corlib/NodeBasic.hpp>
#if   defined(PETSC)
    #include "petsc.h"
    #include <corlib/SystemSolvePetsc.hpp>
#elif defined(BDDCML)
    #include <mpi.h>
    #include <corlib/SystemSolveBddc.hpp>
#endif
#include <corlib/MatrixAssembler.hpp>
#include <corlib/VectorAssembler.hpp>
#include <corlib/Distributor.hpp>
#include <corlib/Accessor.hpp>
#include <corlib/DistributedNode.hpp>
#include <corlib/DistributedElement.hpp>
#include <corlib/Constraints.hpp>
#include <corlib/SmfHead.hpp>
#include <corlib/AccumulateQuantity.hpp>
#include <corlib/misc.hpp>

// solid::material::
#include <solid/material/MaterialContainer.hpp>
// solid::fem
#include <solid/fem/NodeStatic.hpp>
#include <solid/fem/ElementStatic.hpp>


//------------------------------------------------------------------------------
template<typename MESH>
void writeMeshData( const std::string & vtuFile, MESH * mesh, const unsigned rank );

//------------------------------------------------------------------------------
// Dirichlet boundary conditions
corlib::NodalConstraint
dirichletFun( const boost::numeric::ublas::bounded_vector< double, 2 > & X )
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

//------------------------------------------------------------------------------
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
    const double seconds = MPI_Wtime();

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
    inputFile.close( );
    
    // check the number of processes
    int nProc;
    MPI_Comm_size( commAll, &nProc );
    FTL_VERIFY_DESCRIPTIVE( static_cast<unsigned>(nProc) == numberOfProcesses,
                            "Number of MPI processors : %d, "
                            "Number of defined processors in input file : numberOfProcesses=%d\n",
                            static_cast<unsigned>(nProc), numberOfProcesses );

    // read material
    solid::material::MaterialContainer * mc = solid::material::MaterialContainer::instance();
    std::ifstream mf( materialFile.c_str() );
    FTL_VERIFY( mf.is_open() );
    mc -> readMaterialStream( mf );
    mf.close();

    // set up a material
    typedef solid::material::MaterialBase                               MaterialType;
    MaterialType * material = mc -> getMaterial( 0 );


    // set some typedefs
    typedef corlib::Quadrature<corlib::TRIANGLE, 1>                     QuadType;
    typedef corlib::Shapefun<  corlib::TRIANGLE, 3>                     SfunType;
    //typedef corlib::Quadrature<corlib::QUADRILATERAL,4>                 QuadType;
    //typedef corlib::Shapefun<  corlib::QUADRILATERAL,4>                 SfunType;

    typedef corlib::NodeBasic<2>                                        BasisNode;
    typedef solid::fem::NodeStatic<BasisNode>                           SolidNode;
    typedef corlib::DistributedNode<SolidNode>                          NodeType;


    typedef corlib::ElementBasic<NodeType,SfunType>                     BasisElement;
    typedef solid::fem::ElementStatic<BasisElement, MaterialType>       StaticElement;
    typedef corlib::DistributedElement<StaticElement>                   ElementType;
    typedef corlib::Mesh<ElementType>                                   Mesh;
    typedef boost::function<void( ElementType *, const NodeType::VecDim &, 
                                  const double & ) >                    Integrand;

    // unique filename generator
    corlib::UniqueFilename filenameGen(3);

    // read smf file and construct mesh
    int rank;
    MPI_Comm_rank( commAll, &rank );
    std::string meshFile = filenameGen( basename, rank, "psmf" );
    std::ifstream smf( meshFile.c_str( ) );
    FTL_VERIFY( smf.is_open() );
    // skim validator
    corlib::SmfHead smfHead;
    smfHead.readValidated( smf, SfunType::myShape, SfunType::numFunctions );
    // read header for number of local nodes
    unsigned numSubNodes, numSubElements;
    smf >> numSubNodes >> numSubElements;
    smf.close();
    smf.open( meshFile.c_str( ) );
    Mesh mesh( smf );
    smf.close( );

    // cache material to elements
    mesh.iterateOverElements( boost::bind( &ElementType::cacheMaterial, _1, material ) );

    // read the constraints
    //std::string constraintsFile = filenameGen( basename, rank, "constraints" );
    //std::ifstream cstrFile( constraintsFile.c_str( ) );
    //FTL_VERIFY( cstrFile.is_open( ) );
    //corlib::ConstraintsFromFile<NodeType>  constraintsHandler( cstrFile );
    //cstrFile.close();
    //mesh.iterateOverNodes( constraintsHandler );
    
    // Dirichlet boundary conditions
    typedef corlib::NodalConstraint                                     SConstraint;
    typedef boost::function< SConstraint( const SolidNode::VecDim & ) > SFunType;
    typedef corlib::ConstraintsFromFunction< NodeType, SFunType >       Dirichlet;
    Dirichlet dirichlet( boost::bind( &dirichletFun, _1 ) );
    mesh.iterateOverNodes( dirichlet );

    // create global vector of dof starts
    std::map<int,int> numNodeDofsGlobal;
    std::pair<int,int> nums = corlib::generateDofStarts<Mesh> ( mesh, commAll, numNodeDofsGlobal );

    // find number of global dofs
    const unsigned numTotalNodes = nums.first;
    const unsigned numTotalDofs  = nums.second;

    // number dofs with the given array
    mesh.iterateOverNodes( corlib::dofNumberFromStartsFun( &NodeType::numberDOFs, &numNodeDofsGlobal ) );
    // free memory needed
    numNodeDofsGlobal.clear();

    // inform user about basic parameters
    if ( rank == 0 ) {
        std::cout << "number of processes " << numberOfProcesses << std::endl;
        std::cout << "total number of DOFs:          " << numTotalDofs << std::endl;
        std::cout.flush();
    }
    
#if defined(PETSC)
    // count number of local nodes 
    corlib::AccumulateQuantity<NodeType, unsigned> 
        numNodesLocCounter( boost::bind( corlib::getOne ), 0 );
    const int numLocalNodes = 
        mesh.iterateOverNodesWithPredicate( numNodesLocCounter,
                                            boost::bind( &NodeType::isDDActive, _1 ) );

    // number of local dofs
    const unsigned numLocalDofs = numLocalNodes * NodeType::dof;

    // Set up a system matrix + solver
    typedef corlib::SystemSolvePetsc SysSolve;
    SysSolve * sysSolve = new SysSolve( numTotalDofs, numLocalDofs );
    
    // default Petsc options
    sysSolve -> insertOptions( " " );

#elif defined(BDDCML)
    // count number of nodes on subdomain
    corlib::AccumulateQuantity<NodeType, unsigned> 
        numNodesSubCounter( boost::bind( corlib::getOne ), 0 );
    int numNodesSub = mesh.iterateOverNodes( numNodesSubCounter );

    // number of subdomain dofs
    unsigned numDofsSub = numNodesSub * NodeType::dof;

    // Set up a system matrix + solver
    typedef corlib::SystemSolveBddc SysSolve;

    //specify type of matrix ( GENERAL, SPD, or SYMMETRICGENERAL )
    SysSolve::MatrixType matrixType = SysSolve::GENERAL;

    SysSolve * sysSolve = new SysSolve( numTotalDofs, numDofsSub, matrixType ); 

    // load subdomain mesh into BDDCML solver
    sysSolve -> loadMesh( mesh );
#endif

    // quadrature object
    QuadType quadrature;

    // container for the filenames
    corlib::VTUparaAnim animWriter( basename, "vtu", nProc );

    // write VTU file for this process
    const std::string vtuFile = filenameGen( basename, rank, 0, "vtu" );
    writeMeshData( vtuFile, &mesh, rank );

    // store snap shot on proc 0
    if ( rank == 0 ) animWriter.storeSnapshots( 0, 0 );

    for ( unsigned step = 0; step < numberOfSubSteps; step ++ ) {
        if ( rank == 0 ) std::cout << "SUBSTEP: " << step << std::endl;

        const double tolerance = 1.e-5;
        const unsigned maxiter = 10;

        for ( unsigned i = 0; i < maxiter; i ++ ) {
            if ( rank == 0 ) std :: cout << "Iteration number " << i << std :: endl;

            mesh.iterateOverNodes( boost::bind( &NodeType::clearForce, _1 ) );

            // compute nodal  forces
            Integrand internalForce( &ElementType::internalForceIntegrand );
            corlib::Integrator<QuadType,Integrand> integrator2( quadrature, internalForce );
            mesh.iterateOverElements( integrator2 );

            // assemble the  forces to the system matrices rhs
            mesh.iterateOverNodes( corlib::vectorAssemblerFun( &NodeType::giveForce, 
                                                               &NodeType::copyDofArray, 
                                                               sysSolve ) );

            // Compute and assemble element stiffness matrices
            corlib::MatrixComputeAndAssembleFun<const ElementType,QuadType,SysSolve>
                mafStiffness( &ElementType::stiffnessIntegrand, 
                              &ElementType::getDofIndices,
                              &ElementType::getDofIndices,
                              quadrature, sysSolve );
            mesh.iterateOverElements( mafStiffness );


#ifdef PETSC
            // process the matrix and right hand side
            sysSolve -> finishAssembly( );
#endif
            
            // collect and apply constraints
            double factor = ( i==0 ? 1./double(numberOfSubSteps) : 0. );
            corlib::ConstraintFunctor< NodeType, SysSolve > constraint( sysSolve, factor );
            constraint = mesh.iterateOverNodes( constraint );
            constraint.applyConstraints();
            
            // norm of the residual
            const double normResidual = sysSolve -> normRhs( );
            if ( rank == 0 ) std::cout << "Norm of residual  " << normResidual << std::endl;
            //if ( normResidual < tolerance ) break;

            sysSolve -> solveSystem( );
            // get info on solution
            const unsigned numberOfIterations = sysSolve -> giveNumIterations();
            const int convergedReason         = sysSolve -> giveConvergedReason();
            if ( rank == 0 ) {
                std::cout << "solver converged in   " << numberOfIterations << " iterations"
                          << " for reason   " << convergedReason << std::endl;
            }

            // add nodal solutions to increment 
            mesh.iterateOverNodes( corlib::distributorFun( sysSolve, 
                                                           &NodeType::copyDofArray, 
                                                           &NodeType::addToIncrement ) );

            // clear the stiffness matrix and RHS
            sysSolve -> clearMatrix( );
            sysSolve -> clearRhs( );


            const double normDeltaU = sysSolve -> normSol( );
            if ( rank == 0 ) std::cout << "Norm of delta U  " << normDeltaU << std::endl;
            if ( normDeltaU < tolerance ) break;

            sysSolve -> clearSol( );
        }

        // set increment as current displacement
        mesh.iterateOverNodes( boost::bind( &NodeType::updateDisplacements, _1 ) );
    
        // write VTU file for this process
        std::string vtuFile = filenameGen( basename, rank, step+1, "vtu" );
        writeMeshData( vtuFile, &mesh, rank );

        // store snap shot on proc 0
        if ( rank == 0 ) animWriter.storeSnapshots( step+1, step+1 );

    }

    // write pvd file
    if ( rank == 0 ) {
        const std::string animFile = basename + ".pvd";
        std::ofstream anim( animFile.c_str() );
        animWriter.writeAnimationFile( anim );
        anim.close();
    }

    delete sysSolve;

    const double elapsed = MPI_Wtime() - seconds;
    if ( rank == 0 ) std::cout << "Elapsed time: " << corlib::secondsToHMMSS( elapsed ) << std::endl;

#ifdef PETSC
    PetscFinalize();
#endif

    MPI_Finalize();

    return 0;
}

//------------------------------------------------------------------------------
template< typename MESH>
void writeMeshData( const std::string & vtuFile, MESH * mesh, const unsigned rank )
{
    // set up the VTU writer and open a file
    std::ofstream vtu( vtuFile.c_str( ) );
    FTL_VERIFY( vtu.is_open( ) );
    corlib::VTUwriter<MESH> vtuwriter( mesh, vtu );


    // write mesh data
    vtuwriter.writeMesh();
    
    // write point data
    vtuwriter.openPointData();
    vtuwriter.writeNodalQuantity( corlib::accessorFun( &MESH::Node::giveDisplacements, "Displacements" ) );
    vtuwriter.closePointData( );

    // write element data
    vtuwriter.openCellData( );
    corlib::AssignConstantNumber<const typename MESH::Element> procId( rank, "procID" );
    vtuwriter.writeElementQuantity( procId );
    typedef boost::function< typename MESH::Element::Mat3x3( const typename MESH::Element * ) > Mat3x3Fun;
    Mat3x3Fun giveStress = boost::bind( &MESH::Element::firstPiolaKirchhoff, _1,
                                        corlib::ShapeTraits< MESH::Element::myShape >::centroid( ) );
    corlib::Accessor< Mat3x3Fun > accessStress( giveStress, "Stress" );
    vtuwriter.writeElementQuantity( accessStress );
    vtuwriter.closeCellData( );

    // finish writing
    vtuwriter.finalize();
    vtu.close( );
}
