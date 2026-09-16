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

// corlib
#include <corlib/Mesh.hpp>
#include <corlib/ElementBasic.hpp>
#include <corlib/Shapefun.hpp>
#include <corlib/Quadrature.hpp>
#include <corlib/Constraints.hpp>
#include <corlib/Distributor.hpp>
#include <corlib/DofNumberer.hpp>
#include <corlib/Accessor.hpp>
#include <corlib/MatrixAssembler.hpp>
#include <corlib/VectorAssembler.hpp>
#include <corlib/BodyForce.hpp>
#include <corlib/PropertiesParser.hpp>
#include <corlib/VTUwriter.hpp>
#include <corlib/linalg.hpp>
#include <corlib/fuzzyEqual.hpp>
#include <corlib/verify.hpp>
#include <corlib/ElementBasic.hpp>
#include <corlib/NodeBasic.hpp>
#include <corlib/Integrator.hpp>

// decide upon solver
#ifdef PETSC
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
#include <solid/fem/NodalForces.hpp>

#include <iostream>
#include <fstream>
#include <string>
#include <cstdlib>

#include <boost/timer.hpp>
#include <boost/function.hpp>

// helper function to write the VTU output
template<typename MESH>
void writeMeshData( const std::string & vtuFile, MESH * mesh );

//------------------------------------------------------------------------------
// Dirichlet boundary conditions
template<unsigned DIM>
corlib::NodalConstraint
dirichletFun( const boost::numeric::ublas::bounded_vector<double,DIM> X )
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

//------------------------------------------------------------------------------
// Body force function (has to be given the coordinate as argument)
template<unsigned DIM>
boost::numeric::ublas::bounded_vector<double,DIM> 
gravity( const boost::numeric::ublas::bounded_vector<double,DIM> X )
{
    boost::numeric::ublas::bounded_vector<double,DIM> G;
    G.clear();
    G[DIM-1] = -100.;
    return G;
}

//------------------------------------------------------------------------------
template<unsigned DIM>
bool nodalForceFun( const boost::numeric::ublas::bounded_vector<double,DIM> & X,
                    boost::numeric::ublas::bounded_vector<double,DIM> & F )
{
    bool applyAforce = false;
    if ( corlib::fuzzyEqual( X[0], 1. ) and corlib::fuzzyEqual( X[1], 1. ) ) {
        F.clear();
        applyAforce = true;
        F[0] = 30.;
    }
    return applyAforce;
}

//------------------------------------------------------------------------------
int main( int argc, char* argv[] )
{ 
    //! Variable attributes
    const corlib::shape elemShape  = corlib::QUADRILATERAL;
    const unsigned      numNodesPE = 9;
    const unsigned      numGaussP  = 9;

    //! Fixed attributes (some derived from the above)
    const unsigned dim = corlib::ShapeTraits<elemShape>::dim;
    typedef corlib::Shapefun<elemShape,numNodesPE>               Sfun;
    typedef corlib::Quadrature<elemShape,numGaussP>              Quad;
    typedef corlib::NodeBasic<dim>                               BasisNode;
    typedef solid::fem::NodeStatic<BasisNode>                    Node;
    typedef solid::material::MaterialBase                        Material;
    typedef corlib::ElementBasic<Node,Sfun>                      BasisElement;
    typedef solid::fem::ElementStatic<BasisElement,Material>     Element;
    typedef corlib::Mesh<Element>                                Mesh;
    typedef boost::function<void( Element *, 
                                  const Node::VecDim &, 
                                  const double & ) >             Integrand;

    //! Timer for the total time of execution
    boost::timer totalTimer;
    
#ifdef PETSC
    //! Initialize PETSC
    PetscInitialize( &argc, &argv, PETSC_NULL, PETSC_NULL );
#endif

    //! Variables to be read from input file
    std::string meshFile, materialFile, constraintsFile, nodalForcesFile;
    unsigned numberOfSubSteps = 0, numberOfSubIterations = 0;

    //! Feed properties parser with the variables to read
    corlib::PropertiesParser *prop = new corlib::PropertiesParser;
    prop -> registerPropertiesVar( "meshFile",              meshFile   );
    prop -> registerPropertiesVar( "materialFile",          materialFile );
    prop -> registerPropertiesVar( "constraintsFile",       constraintsFile );
    prop -> registerPropertiesVar( "numberOfSubSteps",      numberOfSubSteps );
    prop -> registerPropertiesVar( "nodalForcesFile",       nodalForcesFile );
    prop -> registerPropertiesVar( "numberOfSubIterations", numberOfSubIterations );

    //! Read variables from the input.dat file
    const std::string inputFile = "./input.dat";
    std::ifstream inp( inputFile.c_str()  );
    FTL_VERIFY( inp.is_open() );
    prop -> readValues( inp );
    delete prop;
    inp.close( );

    //! Read material
    solid::material::MaterialContainer * mc = solid::material::MaterialContainer::instance( );
    std::ifstream mf( materialFile.c_str( ) );
    FTL_VERIFY( mf.is_open() );
    mc -> readMaterialStream( mf );
    mf.close( );

    //! Read smf file and construct mesh
    std::ifstream smf( meshFile.c_str( ) );
    FTL_VERIFY( smf.is_open() );
    Mesh mesh( smf );
    smf.close( );

    //! Obtain the 1st material
    Material * material = mc -> getMaterial( 0 );

    //! Cache material to elements
    mesh.iterateOverElements( boost::bind( &Element::cacheMaterial, _1, material ));

    //! Dirichlet boundary conditions
    typedef boost::numeric::ublas::bounded_vector<double,dim>            VecDim;
    typedef boost::function<corlib::NodalConstraint( VecDim ) >          VecDim2Constraint;
    typedef corlib::ConstraintsFromFunction<Node,VecDim2Constraint>      Dirichlet;
    Dirichlet dirichlet( boost::bind( &dirichletFun<dim>, _1 ) );
    mesh.iterateOverNodes( dirichlet );

    //! Number dofs
    const unsigned numdofs = mesh.iterateOverNodes( corlib::dofNumberFun( &Node::numberDOFs, 0 ) );
    std::cout << numdofs << " degrees of freedom" << std :: endl;

    //! Create a solver
#ifdef PETSC
    typedef corlib::SystemSolvePetsc                                   SysSolve;
#elif defined SUPERLU
    typedef corlib::SystemSolveLU                                      SysSolve;
#elif defined PARDISO_MKL
    typedef corlib::SystemSolvePardiso                                 SysSolve;
#else
    typedef corlib::SystemSolve                                        SysSolve;
#endif
    SysSolve * sysmat  = new SysSolve( numdofs );
#ifdef PETSC
    sysmat -> insertOptions( "-ksp_type preonly -pc_type lu"
                             " -pc_factor_mat_solver_package mumps" );
#endif

    //! The quadrature
    Quad quadrature;

    //! Set up an animation writer
    FTL_VERIFY( not system( "mkdir -p animation" ) );
    std::string basename = "./animation/" + meshFile;
    corlib::VTUanim animation( basename, "vtu" );
    
    //! Write data to VTU file
    std::string vtuFile = animation.snapshotName( 0., 0 );
    writeMeshData( vtuFile, &mesh );

    //! Nodal forces from file
    std::ifstream nf( nodalForcesFile.c_str() );
    FTL_VERIFY( nf.is_open() );
    solid::fem::NodalForces<Node> nodalForces( nf );
    nf.close();

    //! Nodal forces from function
    boost::function< bool( const boost::numeric::ublas::bounded_vector<double,dim> &,
                           boost::numeric::ublas::bounded_vector<double,dim> &  ) > 
        nff( boost::bind( nodalForceFun<dim>, _1, _2  ) );
    solid::fem::NodalForces<Node> nodalForces2( mesh, nff );

    //! Show user what nodal force there are
    std::cout << "These are the nodal forces: " << std::endl;
    std::cout << "  "; nodalForces.write(  std::cout );
    std::cout << "  "; nodalForces2.write( std::cout );

    //! Define a body force
    typedef boost::function< VecDim( VecDim ) >              VecDim2VecDim;
    typedef corlib::BodyForce<Element,VecDim2VecDim>         BodyForce;
    BodyForce bodyForce( boost::bind( &gravity<dim>, _1 ) );
    corlib::Integrator<Quad,BodyForce>  bodyForceIntegrator( quadrature, bodyForce );

    //--------------------------------------------------------------------------
    //! Loop over load increments
    for ( unsigned step = 0; step < numberOfSubSteps; step ++ ) {
        std::cout << "-----------------------------------------" << std::endl
                  << "Substep " << step + 1 << " of " << numberOfSubSteps
                  << std::endl;

        const double tolerance = 1e-5;
        
        //----------------------------------------------------------------------
        //! Loop until convergence
        for ( unsigned i = 0; i < numberOfSubIterations; i ++ ) {
            std::cout << "  Iteration number " << i << std :: endl;

            boost::timer timerSetup;
            
            //! Purge the nodal force storage
            mesh.iterateOverNodes( boost::bind( &Node::clearForce, _1 ) );

            //! Apply nodal forces
            const double loadFactor = 
                static_cast<double>(step+1) / static_cast<double>(numberOfSubSteps);
            mesh.iterateOverNodes( boost::bind(  nodalForces, _1, loadFactor ) );
            mesh.iterateOverNodes( boost::bind( nodalForces2, _1, loadFactor ) );

            //! Apply body forces
            bodyForce.setFactor( loadFactor );
            mesh.iterateOverElements( bodyForceIntegrator );

            //! Compute nodal  forces (equation residual)
            Integrand internalForce( &Element::internalForceIntegrand );
            corlib::Integrator<Quad,Integrand> integrator2( quadrature, internalForce );
            mesh.iterateOverElements( integrator2 );

            //! Compute and assemble element stiffness matrices
            corlib::MatrixComputeAndAssembleFun<const Element,Quad,SysSolve>
                mafStiffness( &Element::stiffnessIntegrand, 
                              &Element::getDofIndices,
                              &Element::getDofIndices,
                              quadrature, sysmat );
            mesh.iterateOverElements( mafStiffness );

            //! Assemble the  forces to the system matrices rhs
            mesh.iterateOverNodes( corlib::vectorAssemblerFun( &Node::giveForce, 
                                                               &Node::copyDofArray, 
                                                               sysmat ) );

            //! For petsc only (else a dummy call)
            sysmat -> finishAssembly();

            //! Collect and apply constraints
            const double factor = ( i==0 ? 1./ static_cast<double>(numberOfSubSteps) : 0. );
            std::cout << "    Constraint factor " << factor << std::endl;
            corlib::ConstraintFunctor<Node,SysSolve> constraint( sysmat, factor );
            constraint = mesh.iterateOverNodes( constraint );
            constraint.applyConstraints();

            std::cout << "    " << timerSetup.elapsed() << " s for setup" << std::endl;

            //! Check size of residual
            const double normResidual = sysmat -> normRhs( );
            std::cout << "    Norm of residual  " << normResidual << std :: endl;
            if ( normResidual < tolerance ) break;

            //! Solve the system
            boost::timer timerSolve;
            sysmat -> solveSystem( );
            std::cout << "    " << timerSolve.elapsed() << " s for solve" << std::endl;

            //! Add nodal solutions to increment 
            mesh.iterateOverNodes( corlib::distributorFun( sysmat, 
                                                           &Node::copyDofArray, 
                                                           &Node::addToIncrement ) );

            //! Check size of increment
            const double normDeltaU = sysmat -> normSol( );
            std::cout << "    Norm of delta U  " << normDeltaU << std :: endl;
            if ( normDeltaU < tolerance ) break;

            //! Clear the stiffness matrix and RHS
            sysmat -> clearMatrix( );
            sysmat -> clearRhs( );

        }

        //! New displacement field
        mesh.iterateOverNodes( boost::bind( &Node::updateDisplacements, _1 ) );

        //! Write data to VTU file
        std::string vtuFile = animation.snapshotName( step+1, step+1 );
        writeMeshData( vtuFile, &mesh );

    }

    //------------------------------------------------------------
    //! write animation file
    std::string animationFile( "animation.pvd" );
    std::ofstream anim( animationFile.c_str( ) );
    animation.writeAnimationFile( anim );
    anim.close();

    std::cout << "-----------------------------------------" << std::endl
              << "Total elapsed time in seconds: " << totalTimer.elapsed()
              << std::endl;

    delete sysmat;

#ifdef PETSC
    PetscFinalize();
#endif

    return 0;
}

//------------------------------------------------------------------------------
template<typename MESH>
void writeMeshData( const std::string & vtuFile, MESH * mesh )
{
    // set up the VTU writer and open a file
    std::ofstream vtu( vtuFile.c_str( ) );
    corlib::VTUwriter<MESH> vtuwriter( mesh, vtu );

    //! write mesh data
    vtuwriter.writeMesh();
    
    //-- write point data
    vtuwriter.openPointData();
    vtuwriter.writeNodalQuantity( corlib::accessorFun( &MESH::Node::giveDisplacements, "Displacements" ) );
    vtuwriter.closePointData( );

    // -- write element data
    vtuwriter.openCellData( );
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
