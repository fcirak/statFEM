// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//
//                 Computational Structural Mechanics Lab
//                         University of Cambridge
//
// This software is copyrighted and all rights are retained by the University of
// Cambridge. All use, disclosure, and/or reproduction of any part not expressly
// authorised by the University of Cambridge is prohibited. (C) 2010
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

//! @file DriverPotential.ipp

//------------------------------------------------------------------------------
template< typename MESH, typename SOLVER >
del2::fem::DriverPotential< MESH, SOLVER >::
DriverPotential( const std::string & inputFileName ) 
    : fileBaseName_( "" ),
      conductivity_( 0.0 ),
      mesh_( NULL ),
      dirichlet_( NULL ),
      bodyForce_( NULL ),
      numDofs_( 0 ),
      sysMat_( NULL )
{
    // variables to be read from input file
    std::string meshFile;

    // feed properties parser with the variables to read
    corlib::PropertiesParser * prop = new corlib::PropertiesParser;
    prop -> registerPropertiesVar( "meshFile",         meshFile   );
    prop -> registerPropertiesVar( "conductivity",     conductivity_ );

    // read variables from the input.dat file
    std::ifstream inputFile( inputFileName.c_str( ) );
    FTL_VERIFY_DESCRIPTIVE( inputFile.is_open(), "Cannot open input file" );
    prop->readValues( inputFile );

    // finish reading
    delete prop;
    inputFile.close( );

    // set file base name
    fileBaseName_ = meshFile.substr( 0, meshFile.rfind( "." ) );

    // create mesh
    // read smf file and construct mesh
    std::ifstream smf( meshFile.c_str( ) );
    FTL_VERIFY_DESCRIPTIVE( smf.is_open(), "Canot open SMF file" );
    mesh_ = new Mesh( smf );
    FTL_VERIFY_DESCRIPTIVE( mesh_, "Failed to create mesh" );
    smf.close( );

    //! cache material in elements
    mesh_->iterateOverElements( boost::bind2nd( boost::mem_fun( &Element::setConductivity ), conductivity_ ));

    return;
}

//------------------------------------------------------------------------------
template< typename MESH, typename SOLVER >
del2::fem::DriverPotential< MESH, SOLVER >::~DriverPotential( )
{
    if ( mesh_      ) { delete mesh_;      mesh_      = NULL; }
    if ( dirichlet_ ) { delete dirichlet_; dirichlet_ = NULL; }
    if ( bodyForce_ ) { delete bodyForce_; bodyForce_ = NULL; }
    this->destroySystem( );
    return;
}

//------------------------------------------------------------------------------
template< typename MESH, typename SOLVER >
void del2::fem::DriverPotential< MESH, SOLVER >::
createDirichletConstraints( DirichFun dirichFun )
{
    FTL_VERIFY_DESCRIPTIVE( !dirichlet_, "Cannot allocate dirichlet constraints : already allocated" );
    dirichlet_ = new Dirichlet( dirichFun );
    return;
}

//------------------------------------------------------------------------------
template< typename MESH, typename SOLVER >
void del2::fem::DriverPotential< MESH, SOLVER >::
createBodyForce( BodyForceFun bodyForceFun )
{
    FTL_VERIFY_DESCRIPTIVE( !bodyForce_, "Cannot allocate bodyForce_ : already allocated" );
    bodyForce_ = new BodyForce( bodyForceFun );
    return;
}

//------------------------------------------------------------------------------
template< typename MESH, typename SOLVER >
unsigned del2::fem::DriverPotential< MESH, SOLVER >::
numberDofs( const unsigned seedDof )
{
    numDofs_ = 
        mesh_ -> iterateOverNodes( corlib::dofNumberFun( &Node::numberDOFs ) ) - seedDof;
    return numDofs_;
}

//------------------------------------------------------------------------------
template< typename MESH, typename SOLVER >
void del2::fem::DriverPotential< MESH, SOLVER >::createSystem( )
{
    FTL_VERIFY_DESCRIPTIVE( !sysMat_, "Cannot allocate solver : already allocated" );
    FTL_VERIFY_DESCRIPTIVE( numDofs_, "No degree-of-freedom" );
    sysMat_ = new SysSolve( numDofs_ );
    return;
}

//------------------------------------------------------------------------------
template< typename MESH, typename SOLVER >
void del2::fem::DriverPotential< MESH, SOLVER >::destroySystem( )
{
    if ( sysMat_ ) {
        delete sysMat_;
        sysMat_ = NULL;
    }
    return;
}

//------------------------------------------------------------------------------
template< typename MESH, typename SOLVER >
void del2::fem::DriverPotential< MESH, SOLVER >::clearSystem( )
{
    // clear the stiffness matrix and RHS
    sysMat_->clearMatrix( );
    sysMat_->clearRhs( );
    return;
}

//------------------------------------------------------------------------------
template< typename MESH, typename SOLVER >
void del2::fem::DriverPotential< MESH, SOLVER >::assembleResidual( )
{
    // clear nodally distributed forces
    mesh_->iterateOverNodes( std::mem_fun( &Node::clearForce ) );

    // external forces due to body load
    if ( bodyForce_ ) {
        corlib::Integrator<Quad,BodyForce> bodyForceIntegrator( quadrature_, 
                                                                *bodyForce_ );
        mesh_->iterateOverElements( bodyForceIntegrator );
    }

    // internal forces -- do nothing, since linear problem,
    // and we start from flux-free configuration
    ;

    // assemble the  forces to the system matrices rhs
    mesh_ -> iterateOverNodes( corlib::vectorAssemblerFun( &Node::getForce, 
                                                           &Node::copyDofArray, 
                                                           sysMat_ ) );
    // done
    return;
}

//------------------------------------------------------------------------------
template< typename MESH, typename SOLVER >
void del2::fem::DriverPotential< MESH, SOLVER >::assembleJacobian( )
{
    // tangent of body force --  do nothing, since assumed constant
    ; 

    // tangent of internal forces
    corlib::MatrixComputeAndAssembleFun<const Element,Quad,SysSolve>
        maf( &Element::stiffnessIntegrand, 
             &Element::getDofIndices,
             &Element::getDofIndices,
             quadrature_, sysMat_ );
    mesh_ -> iterateOverElements( maf );

    // done
    return;
}

//------------------------------------------------------------------------------
template< typename MESH, typename SOLVER >
void del2::fem::DriverPotential< MESH, SOLVER >::assembleFinish( )
{
    sysMat_->finishAssembly( );
    return;
}

//------------------------------------------------------------------------------
template< typename MESH, typename SOLVER >
void del2::fem::DriverPotential< MESH, SOLVER >::dirichletConstrainSystem( )
{
    if ( dirichlet_ ) {
        // distribute constraints due to Dirichlet function onto nodes
        mesh_->iterateOverNodes( *dirichlet_ );
        // apply Dirichlet constraints to global linear system
        corlib::ConstraintFunctor< Node, SysSolve > constraint( sysMat_ );
        constraint = mesh_->iterateOverNodes( constraint );
        constraint.applyConstraints( );
    }
    return;
}

//------------------------------------------------------------------------------
template< typename MESH, typename SOLVER >
void del2::fem::DriverPotential< MESH, SOLVER >::solveSystem( )
{
    sysMat_->solveSystem( ); // has to conform
    return;
}

//------------------------------------------------------------------------------
template< typename MESH, typename SOLVER >
void del2::fem::DriverPotential< MESH, SOLVER >::updateIncrement( )
{
    //! add nodal solutions to increment 
    mesh_ -> iterateOverNodes( corlib::distributorFun( sysMat_, 
                                                       &Node::copyDofArray, 
                                                       &Node::setIncrement ) );
    return;
}

//------------------------------------------------------------------------------
template< typename MESH, typename SOLVER >
void del2::fem::DriverPotential< MESH, SOLVER >::updateStep( )
{
    mesh_->iterateOverNodes( std::mem_fun( &Node::updatePotential ) );
    return;
}

//------------------------------------------------------------------------------
template< typename MESH, typename SOLVER >
void del2::fem::DriverPotential<MESH, SOLVER>::
writeVtu( const std::string & vtuFile )
{
    // set up the VTU writer and open a file
    std::ofstream vtu( vtuFile.c_str() );
    FTL_VERIFY_DESCRIPTIVE( vtu.is_open( ), "Cannot open VTU file for writing" );
    corlib::VTUwriter< Mesh > vtuwriter( mesh_, vtu );

    // write mesh data
    vtuwriter.writeMesh( );
    
    // write point data
    vtuwriter.openPointData( );
    vtuwriter.writeNodalQuantity( corlib::accessorFun( &Node::getPotential, "Potential" ) );
    vtuwriter.closePointData( );
    
    // finish writing
    vtuwriter.finalize( );
    vtu.close( );

    // done
    return;
}

