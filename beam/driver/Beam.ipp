// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//
//                 Computational Structural Mechanics Lab
//                         University of Cambridge
//
// This software is copyrighted and all rights are retained by the University of
// Cambridge. All use, disclosure, and/or reproduction of any part not expressly
// authorised by the University of Cambridge is prohibited. (C) 2010
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

//! @file Beam.ipp

#include <fstream>

//------------------------------------------------------------------------------
//! Compute elements' mass matrices
template<typename ELEMENT, typename SOLVER, unsigned NGP>
void beam::driver::Beam<ELEMENT,SOLVER,NGP>::computeElementMass()
{
    mesh_.iterateOverElements( std::bind( &Element::clearMass, std::placeholders::_1 ) );
    Integrand mass( &Element::massIntegrand );
    Integrator massIntegrator( quadrature_, mass );
    mesh_.iterateOverElements( massIntegrator );
}

//------------------------------------------------------------------------------
//! Compute elements' stiffness matrices
template<typename ELEMENT, typename SOLVER, unsigned NGP>
void beam::driver::Beam<ELEMENT,SOLVER,NGP>::computeElementStiffness()
{
    mesh_.iterateOverElements( std::bind( &Element::clearStiffness, std::placeholders::_1 ) );
    Integrand stiffness( &Element::stiffnessIntegrand );
    Integrator stiffnessIntegrator( quadrature_, stiffness );
    mesh_.iterateOverElements( stiffnessIntegrator );
}

//------------------------------------------------------------------------------
//! Assemble stiffness matrix
template<typename ELEMENT, typename SOLVER, unsigned NGP>
void beam::driver::Beam<ELEMENT,SOLVER,NGP>::assembleStiffness( )
{
    solver_ -> clearMatrix();

    corlib::MatrixGiveAndAssembleFun<const ELEMENT,Solver>
        giveStiffness( &Element::addStiffnessMatrix, 
                       &Element::getDofIndices, 
                       &Element::getDofIndices,
                       solver_ );
    mesh_.iterateOverElements( giveStiffness );
}

//------------------------------------------------------------------------------
//! Assemble matrix for Newmark method
template<typename ELEMENT, typename SOLVER, unsigned NGP>
void beam::driver::Beam<ELEMENT,SOLVER,NGP>::assembleMatricesNewmark( const double stepSize )
{
    solver_ -> clearMatrix();
    MatrixDonator massMat(      &Element::addMassMatrix );
    MatrixDonator stiffnessMat( &Element::addStiffnessMatrix );
    corlib::NewmarkMatrixAssembler<Solver,Element,MatrixDonator,MatrixDonator> 
        nma( newmarkBeta_, newmarkGamma_, stepSize, massMat, stiffnessMat, 
             rayleighM_, rayleighK_ );

    mesh_.iterateOverElements( std::bind( nma, std::placeholders::_1, solver_ )  );
}

//------------------------------------------------------------------------------
//! Assemble matrix for Euler Backward method
template<typename ELEMENT, typename SOLVER, unsigned NGP>
void beam::driver::Beam<ELEMENT,SOLVER,NGP>::assembleMatricesBE( const double stepSize )
{
    solver_ -> clearMatrix();
    MatrixDonator massMat( &Element::addMassMatrix );
    MatrixDonator stiffnessMat( &Element::addStiffnessMatrix );
    corlib::BackwardEulerMatrixAssembler<Solver,Element,MatrixDonator,MatrixDonator> 
        bea( stepSize, massMat, stiffnessMat, rayleighM_, rayleighK_ );

    mesh_.iterateOverElements( std::bind( bea, std::placeholders::_1, solver_ )  );
}

//------------------------------------------------------------------------------
template<typename ELEMENT, typename SOLVER, unsigned NGP>
void beam::driver::Beam<ELEMENT,SOLVER,NGP>::
computeAndAssembleLHS( const double stepSize )
{
    // compute LHS
    this->computeElementMass();
    this->computeElementStiffness();

    // assemble LHS (these methods clear sys matrix first)
    switch ( timeIntegrator_ ) {
    case TI_STATIC:
        this -> assembleStiffness();
        break;
    case TI_BACKWARDEULER:
        FTL_VERIFY( stepSize > 0. );
        this -> assembleMatricesBE(      stepSize );
        break;
    case TI_NEWMARK:
        FTL_VERIFY( stepSize > 0. );
        this -> assembleMatricesNewmark( stepSize );
        break;
    default:
        FTL_VERIFY( false );
        break;
    }
}

//==============================================================================

//------------------------------------------------------------------------------
//! Compute forces (internal and body force)
template<typename ELEMENT, typename SOLVER, unsigned NGP>
void beam::driver::Beam<ELEMENT,SOLVER,NGP>::clearForces( )
{
    mesh_.iterateOverNodes( std::bind( &Node::clearForce, std::placeholders::_1 ) );
}

//------------------------------------------------------------------------------
//! Compute nodal forces due to external applied forces
template<typename ELEMENT, typename SOLVER, unsigned NGP>
void beam::driver::Beam<ELEMENT,SOLVER,NGP>::
computeExternalForces( )
{
    if ( bodyForceFun_ ) {
        BodyForce bodyForce( bodyForceFun_ );
        corlib::Integrator<Quadrature,BodyForce> 
            bodyForceIntegrator( quadrature_, bodyForce );
        mesh_.iterateOverElements( bodyForceIntegrator );
    }

    if ( externalTractionFun_ ) {
        corlib::Integrator<Quadrature,Integrand>
            etIntegrator( quadrature_, externalTractionFun_ );
        mesh_.iterateOverElements( etIntegrator );
    }
}

//------------------------------------------------------------------------------
//! Compute residual forces
template<typename ELEMENT, typename SOLVER, unsigned NGP>
void beam::driver::Beam<ELEMENT,SOLVER,NGP>::
computeInternalForces( )
{
    Integrand internalForce( &Element::internalForceIntegrand );
    Integrator integrator( quadrature_, internalForce );
    mesh_.iterateOverElements( integrator );
}

//------------------------------------------------------------------------------
//! Compute dynamic forces a la Newmark
template<typename ELEMENT, typename SOLVER, unsigned NGP>
void beam::driver::Beam<ELEMENT,SOLVER,NGP>::
computeForcesNewmark( const double stepSize )
{
    MatrixDonator massMat(  &Element::addMassMatrix );
    MatrixDonator stiffMat( &Element::addStiffnessMatrix );
    corlib::NewmarkInertiaForceAssembler<Solver,Element,MatrixDonator,MatrixDonator> 
        nfa( newmarkBeta_, newmarkGamma_, stepSize, massMat, stiffMat,
             rayleighM_, rayleighK_ );
    mesh_.iterateOverElements( std::bind( nfa, std::placeholders::_1, solver_ ) );
}

//------------------------------------------------------------------------------
//! Compute dynamic forces a la Euler backward
template<typename ELEMENT, typename SOLVER, unsigned NGP>
void beam::driver::Beam<ELEMENT,SOLVER,NGP>::computeForcesBE( const double stepSize )
{
    MatrixDonator massMat(  &Element::addMassMatrix );
    MatrixDonator stiffMat( &Element::addStiffnessMatrix );
    corlib::BackwardEulerInertiaForceAssembler<Solver,Element,MatrixDonator,MatrixDonator> 
        befa( stepSize, massMat, stiffMat, rayleighM_, rayleighK_ );
    mesh_.iterateOverElements( std::bind( befa, std::placeholders::_1, solver_ ) );
}

//------------------------------------------------------------------------------
//! Assemble forces from nodes to the system solver
template<typename ELEMENT, typename SOLVER, unsigned NGP>
void beam::driver::Beam<ELEMENT,SOLVER,NGP>::assembleForces()
{
    solver_ -> clearRhs( );

    mesh_.iterateOverNodes( corlib::vectorAssemblerFun( &Node::giveForce,
                                                        &Node::copyDofArray,
                                                        solver_ ) );
}

//------------------------------------------------------------------------------
template<typename ELEMENT, typename SOLVER, unsigned NGP>
void beam::driver::Beam<ELEMENT,SOLVER,NGP>::
computeAndAssembleRHS( const double stepSize )
{
    // compute and assemble RHS
    this -> clearForces( );
    this -> computeInternalForces();
    this -> computeExternalForces();
    this -> assembleForces();  // clears global RHS first

    switch ( timeIntegrator_ ) {
    case TI_BACKWARDEULER:
        FTL_VERIFY( stepSize > 0. );
        this -> computeForcesBE(      stepSize );
        break;
    case TI_NEWMARK:
        FTL_VERIFY( stepSize > 0. );
        this -> computeForcesNewmark( stepSize );
        break;
    default:
        break;  // do nothing
    }
}

//------------------------------------------------------------------------------
//! Assemble link related matrices and vectors
template<typename ELEMENT, typename SOLVER, unsigned NGP>
void beam::driver::Beam<ELEMENT,SOLVER,NGP>::assembleLinks()
{
    //! assemble link matrices to the system matrix
    corlib::MatrixGiveAndAssembleFun<const Link,Solver>
        couplingMat( &Link::giveCouplingMatrix, 
                     &Link::getDofIndices, 
                     &Link::getDofIndicesP,
                     solver_ );
    mesh_.iterateOverLinks( couplingMat );

    corlib::MatrixGiveAndAssembleFun<const Link,Solver>
        couplingMatT( &Link::giveCouplingMatrixT, 
                      &Link::getDofIndicesP, 
                      &Link::getDofIndices,
                      solver_ );
    mesh_.iterateOverLinks( couplingMatT );

    mesh_.iterateOverLinks( corlib::vectorAssemblerFun( &Link::giveMultiplierRhs,
                                                        &Link::getDofIndices,
                                                        solver_ ) );
    mesh_.iterateOverLinks( corlib::vectorAssemblerFun( &Link::giveLinkResiduum,
                                                        &Link::getDofIndicesP,
                                                        solver_ ) );
}

//==============================================================================

//------------------------------------------------------------------------------
//! Pass back the solution to the nodes and links
template<typename ELEMENT, typename SOLVER, unsigned NGP>
void beam::driver::Beam<ELEMENT,SOLVER,NGP>::distributeSolution()
{
    mesh_.iterateOverNodes( corlib::distributorFun( solver_,
                                                    &Node::copyDofArray, 
                                                    &Node::addToIncrement ) );
    mesh_.iterateOverLinks( corlib::distributorFun( solver_, 
                                                    &Link::getDofIndicesP, 
                                                    &Link::addToIncrement ) );
}

//------------------------------------------------------------------------------
template<typename ELEMENT, typename SOLVER, unsigned NGP>
double beam::driver::Beam<ELEMENT,SOLVER,NGP>::maxDisplacementIncrement() const
{
    corlib::AccumulateQuantity<const Node,double,corlib::Max> 
        findMax( std::bind( &beam::fem::normMax<VecDof>,
                              std::bind( &Node::giveIncrement, std::placeholders::_1 ) ), 0. );
    const double maxIncr = std::for_each( mesh_.nodesBegin(), mesh_.nodesEnd(),
                                          findMax );
    return maxIncr;
}

//------------------------------------------------------------------------------
//! Update a la Newmark
template<typename ELEMENT, typename SOLVER, unsigned NGP>
void beam::driver::Beam<ELEMENT,SOLVER,NGP>::updateNewmark( const double stepSize )
{
    // collect nodal solutions and update
    corlib::NewmarkSolutionDistributor<Node> collector( newmarkBeta_, 
                                                        newmarkGamma_, 
                                                        stepSize  );
    mesh_.iterateOverNodes( collector );
}

//------------------------------------------------------------------------------
//! Update a la Euler backward
template<typename ELEMENT, typename SOLVER, unsigned NGP>
void beam::driver::Beam<ELEMENT,SOLVER,NGP>::updateBE( const double stepSize )
{
    corlib::BackwardEulerSolutionDistributor<Node> collector( stepSize  );
    mesh_.iterateOverNodes( collector );
}

//------------------------------------------------------------------------------
//! Update displacements, velocities, accelerations and other primary unknowns
template<typename ELEMENT, typename SOLVER, unsigned NGP>
void beam::driver::Beam<ELEMENT,SOLVER,NGP>::dynamicUpdate( const double stepSize )
{
    if      ( timeIntegrator_ == TI_BACKWARDEULER ) this -> updateBE(      stepSize );
    else if ( timeIntegrator_ == TI_NEWMARK )       this -> updateNewmark( stepSize );
}

//------------------------------------------------------------------------------
//! Update nonlinear solutions
template<typename ELEMENT, typename SOLVER, unsigned NGP>
void beam::driver::Beam<ELEMENT,SOLVER,NGP>::updateSolution( )
{
    //! update displacements
    mesh_.iterateOverNodes( std::bind( &Node::updateDisplacements, std::placeholders::_1 ) );
    // update Lagrange multipliers
    mesh_.iterateOverLinks( std::bind( &Link::updateMultiplier, std::placeholders::_1 ) );

}

//==============================================================================

//------------------------------------------------------------------------------
//! Write a VTU file
template<typename ELEMENT, typename SOLVER, unsigned NGP>
void beam::driver::Beam<ELEMENT,SOLVER,NGP>::writeMeshData( const std::string & vtuFile )
{
    // set up the VTU writer and open a file
    std::ofstream vtu( vtuFile.c_str( ) );
    FTL_VERIFY( vtu.is_open( ) );

    corlib::VTUwriter<Mesh> vtuwriter( &mesh_, vtu );

    //! write mesh data
    vtuwriter.writeMesh( );
    
    //-- write point data
    vtuwriter.openPointData( );
    vtuwriter.writeNodalQuantity( corlib::accessorFun( &Node::giveDisplacements, "Displacements" ) );
    detail_::WriteDynamicSolution<corlib::VTUwriter<Mesh>,Node,
                                  detail_::IsDynamic<Node>::result>()( vtuwriter );
    vtuwriter.writeNodalQuantity( corlib::accessorFun( &Node::giveForce,    "Forces" ) );
    vtuwriter.closePointData( );

    // finish writing
    vtuwriter.finalize( );
    vtu.close( );

    return;
}

//------------------------------------------------------------------------------
template<typename ELEMENT, typename SOLVER, unsigned NGP>
bool beam::driver::Beam<ELEMENT,SOLVER,NGP>::iterate( const double timeStepSize,
                                                      const unsigned maxIterations,
                                                      DoubleToBool convCheckRhs,
                                                      DoubleToBool convCheckSolution )
{
    if ( not this->getSolver() ) {
        std::pair<unsigned,unsigned> aux = this->numberDofs();
        this->allocateSolver( aux.first + aux.second );
    }

    std::cout << "    Beam solve:" << std::endl;

    // iteration loop
    bool converged = false;
    for ( unsigned i = 0; i < maxIterations and not converged; ++i ) {
            
         std::cout << "  Iteration " << i << std::endl;

         //! Compute (dynamic) LHS 
         this -> computeAndAssembleLHS( timeStepSize );
         //! Compute (dynamic) RHS
         this -> computeAndAssembleRHS( timeStepSize );
         //! Assemble Lagrange multiplier terms
         this -> assembleLinks();

         // norm of the residual
         const double normResidual = this->getSolver()->normRhs();
         std::cout << "    Norm of residual forces (RHS) " << normResidual << std::endl;
         converged = convCheckRhs( normResidual );

         // solve system
         this -> solveSystem();

         this -> distributeSolution();

         // check residual displacements
         const double normDeltaU = this->normSolution();
         std::cout << "    Norm of residual displacements " << normDeltaU << std::endl;
         converged = converged or convCheckSolution( normDeltaU );
    }

    // inform user
    if ( not converged ) {
        std::cout
            << "WARNING: Not converged in " << maxIterations << " iterations."
            << " Continuing." << std::endl;
    }    

    return converged;
}
