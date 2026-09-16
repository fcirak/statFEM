// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//
//                 Computational Structural Mechanics Lab
//                         University of Cambridge
//
// This software is copyrighted and all rights are retained by the University of
// Cambridge. All use, disclosure, and/or reproduction of any part not expressly
// authorised by the University of Cambridge is prohibited. (C) 2010
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

//! @file Static.ipp

//------------------------------------------------------------------------------
template< typename MESH, typename SOLVER, unsigned NGP >
void solid::driver::Static<MESH,SOLVER,NGP>::computeForces()
{
    const double targetTime = this -> getNextTime();

    double loadFactor = nodalForcesFac_( targetTime );
    this -> computeNodalLoads( loadFactor );
    
    if ( bodyForceFun_ != boost::bind( &detail_::zeroVecFun<VecDof,VecDim>, _1 ) ) {
        loadFactor = bodyForceFac_( targetTime );
        this -> computeBodyLoad( bodyForceFun_, loadFactor, bodyForcePred_ );
    }
    
    this -> computeInternalForces();
}

//------------------------------------------------------------------------------
template< typename MESH, typename SOLVER, unsigned NGP >
bool solid::driver::Static<MESH,SOLVER,NGP>::iterate( const unsigned maxIterations,
                                                      DoubleToBool convCheckRhs,
                                                      DoubleToBool convCheckSolution )
{
    if ( not solver_ )
        this -> allocateSolver();

    const double targetTime = this -> getNextTime();

    // iteration loop
    bool converged = false;
    for ( unsigned i = 0; i < maxIterations; ++i ) {
            
        this -> accessOutStream_()
            << "  Iteration number " << i << std::endl;

        this -> clearForces();

        this -> computeAndAssembleForces();

        this -> computeAndAssembleStiffness();

        this -> finishAssembly();

        const double factor = ( i==0 ) ? constraintsFac_( targetTime ) : 0.;
        this -> applyConstraints( factor );

        // norm of the residual
        const double normResidual = this -> normRhs();
        this -> accessOutStream_()
            << "    Norm of residual forces (RHS) " << normResidual << std::endl;
        if ( convCheckRhs( normResidual ) ) {
            converged = true;
            break;
        }

        this -> solveSystem();

        this -> informAboutSolve();

        this -> distributeSolution();

        const double normDeltaU = this -> normSolution();
        this -> accessOutStream_()
            << "    Norm of residual displacements  " << normDeltaU << std::endl;
        if ( convCheckSolution( normDeltaU ) ) {
            converged = true;
            break;
        }

        this -> clearSystem();

    }

    // inform user
    if ( not converged ) {
        this -> accessOutStream_()
            << "WARNING: Not converged in " << maxIterations << " steps."
            << " Continuing." << std::endl;
    }
    return converged;
}

//------------------------------------------------------------------------------
template< typename MESH, typename SOLVER, unsigned NGP >
void solid::driver::Static<MESH,SOLVER,NGP>::advance( const double stepSize,
                                                      const unsigned maxSteps,
                                                      const unsigned maxIterations,
                                                      DoubleToBool convCheckRhs,
                                                      DoubleToBool convCheckSolution )
{
    if ( not solver_ )
        this -> allocateSolver();

    timeStepSize_ = stepSize;

    for ( unsigned step = 0; step < maxSteps; ++step ) {

        this -> accessOutStream_()
            << "-----------------------------------------" << std::endl
            << "Step " << step+1 << " of " << maxSteps << std::endl;
        this -> accessOutStream_()
            << "    Targeting control factor "
            << this -> getNextTime() << std::endl;

        this -> iterate( maxIterations, convCheckRhs, convCheckSolution );

        this -> updateSolution();

        this -> writeAnimationStep( this -> getNextTime(), step+1 );

        this -> updateTime( stepSize );

    }
}

//------------------------------------------------------------------------------
template< typename MESH, typename SOLVER, unsigned NGP >
void solid::driver::Static<MESH,SOLVER,NGP>::writeMeshData( const std::string & vtuFile )
{
    // set up the VTU writer and open a file
    std::ofstream vtu( vtuFile.c_str( ) );
    FTL_VERIFY( vtu.is_open() );
    VTUwriter vtuwriter( &mesh_, vtu );

    // write mesh data
    vtuwriter.writeMesh();
    
    // write point data
    vtuwriter.openPointData();
    this -> writeNodeData( vtuwriter );
    vtuwriter.closePointData( );

    // write element data
    vtuwriter.openCellData( );
    this -> writeElementData( vtuwriter );
    vtuwriter.closeCellData( );

    // finish writing
    vtuwriter.finalize();
    vtu.close();
}

//------------------------------------------------------------------------------
template< typename MESH, typename SOLVER, unsigned NGP >
void solid::driver::Static<MESH,SOLVER,NGP>::writeNodeData( VTUwriter & vtuwriter ) const
{
    vtuwriter.writeNodalQuantity( corlib::accessorFun( &Node::giveDisplacements,
                                                       "Displacements" ) );
}

//------------------------------------------------------------------------------
template< typename MESH, typename SOLVER, unsigned NGP >
void solid::driver::Static<MESH,SOLVER,NGP>::writeElementData( VTUwriter & vtuwriter ) const
{
    typedef boost::function< typename Element::Mat3x3( const Element * ) > Mat3x3Fun;
    Mat3x3Fun giveStress = boost::bind( &Element::firstPiolaKirchhoff, _1,
                                        corlib::ShapeTraits< Element::myShape >::centroid( ) );
    corlib::Accessor< Mat3x3Fun > accessStress( giveStress, "Stress" );
    vtuwriter.writeElementQuantity( accessStress );
}
