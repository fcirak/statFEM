// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//
//                 Computational Structural Mechanics Lab
//                         University of Cambridge
//
// This software is copyrighted and all rights are retained by the University of
// Cambridge. All use, disclosure, and/or reproduction of any part not expressly
// authorised by the University of Cambridge is prohibited. (C) 2010
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

//! @file Static.hpp

#ifndef solid_driver_static_h
#define solid_driver_static_h

#include <iostream>
#include <fstream>

// includes from corlib
#include <corlib/verify.hpp>
#include <corlib/Accessor.hpp>
#include <corlib/VTUwriter.hpp>
#include <corlib/fuzzyEqual.hpp>
#include <corlib/Quadrature.hpp>
#include <corlib/Integrator.hpp>
#include <corlib/DofNumberer.hpp>
#include <corlib/Distributor.hpp>
#include <corlib/MatrixAssembler.hpp>
#include <corlib/VectorAssembler.hpp>
#include <corlib/BodyForce.hpp>
#include <corlib/Constraints.hpp>
#include <corlib/UniqueFilename.hpp>
#include <corlib/Mesh.hpp>

#include <solid/fem/QuadratureTraitsStatic.hpp>
#include <solid/fem/NodalForces.hpp>

//==============================================================================
namespace solid {
    namespace driver {

        namespace ublas = boost::numeric::ublas;

        template< typename ELEMENT, typename SOLVER, unsigned NGP >
        class Static;

        /// namespace for helper functions
        namespace detail_ {
            
            /// Helper function returning a N-dim zero vector given any M-dim vector
            template< typename VECOUT, typename VECIN >
            VECOUT zeroVecFun( VECIN x )
            {
                VECOUT zero; zero.clear(); return zero;
            }

            /// Helper function delivering a constant value stored at constructing
            template< typename X, typename Y >
            struct ReturnConstant : std::unary_function<X,Y>
            {
                ReturnConstant( const Y y ) : y_( y ) { }
                Y operator()( const X x ) const { return y_; }
                const Y y_;  ///< the constant to return
            };

            /// Helper method to allocate global system / solver
            ///
            /// This method is needed because the PETSc and BDDCML solver does not have
            /// a constructor with the used signature. Thus specialisation
            /// of this helper allows to hide the call to a driver
            /// using the corlib::SystemSolveBddc solver.
            template< typename SOLVER >
            void allocateSolver( SOLVER * & solver, const unsigned numDof )
            {
                solver = new SOLVER( numDof );
            }

        }

    }
}

//==============================================================================
/// Basic driver for static solid computations
///
/// \tparam ELEMENT  Solid element
/// \tparam SOLVER   Linear algebraic solver
template< typename ELEMENT, typename SOLVER, unsigned NGP = 0 >
class solid::driver::Static
{
public:
    typedef ELEMENT                                           Element;
    static const enum corlib::shape myShape  = Element::myShape;
    static const unsigned           numNodes = Element::numNodes;

    typedef typename Element::Node                            Node;
    typedef typename Node::VecDim                             VecDim;
    typedef typename Node::VecDof                             VecDof;
    static const unsigned dof = Node::dof;
    static const unsigned dim = Node::dim;

    typedef corlib::Mesh<Element>                             Mesh;

    // predicates (conditionals)
    typedef boost::function<bool( Element * )>                ElementPred;
    typedef boost::function<VecDof( const VecDim )>           VecDimToVecDof;
    typedef boost::function<double( const double )>           DoubleToDouble;
    typedef boost::function<bool( const double )>             DoubleToBool;

    // constraints
    typedef corlib::NodalConstraint                           NodalConstraint;
    typedef boost::function<NodalConstraint( const VecDim )>  DiriFun;

    // load
    typedef solid::fem::NodalForces<Node>                     NodalForces;
    typedef boost::function<bool( const VecDim &, VecDim & )> NodalForceFun;
    typedef corlib::BodyForce<Element,VecDimToVecDof>         BodyForce;

    // quadrature
    typedef solid::fem::QuadratureTraitsStatic<myShape,numNodes> DriverTraits;
    static const unsigned numQuadPoints = (NGP != 0 ? NGP : DriverTraits::numPoints);
    typedef corlib::Quadrature<myShape,numQuadPoints>         Quadrature;
    typedef boost::function<void( Element *,
                                  const VecDim &,
                                  const double & )>           Integrand;
    typedef corlib::Integrator<Quadrature,Integrand,Element>  Integrator;

    typedef SOLVER                                            Solver;

    typedef typename Element::Material                        Material;

    typedef typename corlib::VTUwriter<Mesh>                  VTUwriter;
    typedef corlib::VTUanim                                   VTUanim;

private:
    typedef solid::driver::Static<ELEMENT,SOLVER>             DriverStatic_;

public:

    //--------------------------------------------------------------------------
    /// @name Coarse level methods
    //@{

    /// Default constructor
    Static( std::istream & smf )
        : mesh_( smf )
    { 
        this -> initialise_();
    }

    /// Constructor without mesh
    Static()
        : mesh_()
    { 
        this -> initialise_();
    }


protected:
    /// Common function to initialise member variables
    void initialise_()
    {
        constraintsFac_ = detail_::ReturnConstant<double,double>( 1. );
        nodalForcesFac_ = detail_::ReturnConstant<double,double>( 1. );
        bodyForceFun_   = boost::bind( &detail_::zeroVecFun<VecDof,VecDim>, _1 );
        bodyForceFac_   = detail_::ReturnConstant<double,double>( 1. );
        bodyForcePred_  = corlib::Positive();
        numDofs_        = 0;
        solver_         = NULL;
        time_           = 0.;
        timeStepSize_   = 1.;
        animation_      = NULL;
    }

public:
    /// Destructor
    ~Static()
    { 
        if ( solver_    ) this -> deAllocateSolver();
        if ( animation_ ) this -> deAllocateAnimation();
    }

    /// Set material in elements
    void setMaterial( Material * material )
    {
        mesh_.iterateOverElements( boost::bind( &Element::cacheMaterial, _1, material ) );
    }

    /// Apply Dirichlet BCs by 'constraining' DOFs on nodes.
    /// The constrained DOFs are read from file.
    void setNodeConstraints( std::istream & cstrFile )
    {
        typedef corlib::ConstraintsFromFile<Node>  ConstraintsFF;
        ConstraintsFF constraintsHandler( cstrFile );
        mesh_.iterateOverNodes( constraintsHandler );
    }

    /// Apply Dirichlet BCs by 'constraining' DOFs on nodes
    /// The constrained are chosen by their physical co-ordinate using a function.
    void setNodeConstraints( DiriFun diriFun )
    {
        typedef corlib::ConstraintsFromFunction<Node,DiriFun> Dirichlet;
        Dirichlet dirichlet( diriFun );
        mesh_.iterateOverNodes( dirichlet );
    }

    /// Set Dirichlet BCs factor (displacement control)
    void setNodeConstraintsFactor( DoubleToDouble factorFun )
    {
        constraintsFac_ = factorFun;
    }

    /// Add nodally concentrated loads from file stream
    void addNodalLoads( std::istream & inp )
    {
        nodalForces_.add( inp );
    }

    /// Add nodally concentrated loads with function
    void addNodalLoads( NodalForceFun nodalForceFun )
    {
        nodalForces_.add( mesh_, nodalForceFun );
    }

    /// Set function to scale concentrated nodal loads
    void setNodalLoadsFactor( DoubleToDouble factorFun )
    {
        nodalForcesFac_ = factorFun;
    }

    /// Set body load
    void setBodyLoad( VecDimToVecDof forceFun,
                      ElementPred condition = corlib::Positive() )
    {
        bodyForceFun_ = forceFun;
        bodyForcePred_ = condition;
    }

    /// Set body load factor
    void setBodyLoadFactor( DoubleToDouble factorFun )
    {
        bodyForceFac_ = factorFun;
    }

    /// Create animation writer
    virtual void allocateAnimation( const std::string basename )
    {
        animation_ = new VTUanim( basename, "vtu" );
    }

    /// Destroy animation writer
    virtual void deAllocateAnimation() 
    {
        delete animation_; animation_ = NULL;
    }

    /// Do Newton-Raphson iteration (NRI)
    ///
    /// \param[in]  maxIterations      Maximally permitted iterations
    /// \param[in]  convCheckRhs       Convergence check applied to residual of RHS
    /// \param[in]  convCheckSolution  Convergence check applied to residual displacements
    /// \return                        True if NRI converged
    bool iterate( const unsigned maxIterations,
                  DoubleToBool convCheckRhs,
                  DoubleToBool convCheckSolution );

    /// Do displacement/load control step loop
    ///
    /// \param[in]  stepSize           (Time) step size to advance (eg load increment)
    /// \param[in]  maxSteps           Number of steps to go
    /// \param[in]  maxIterations      Maximally permitted NR iterations
    /// \param[in]  convCheckRhs       Convergence check applied to residual of RHS
    /// \param[in]  convCheckSolution  Convergence check applied to residual displacements
    void advance( const double stepSize,
                  const unsigned maxSteps,
                  const unsigned maxIterations,
                  DoubleToBool convCheckRhs,
                  DoubleToBool convCheckSolution );

    //@}

    //--------------------------------------------------------------------------
    /// @name Medium level methods
    //@{ 

    /// Set step size
    void setTimeStepSize( const double stepSize ) { timeStepSize_ = stepSize; } 

    /// Return time (last converged state)
    double getTime() const { return time_; }

    /// Return target time \f$t_{n+1} = t_n + \Delta t_n\f$
    double getNextTime() const { return time_ + timeStepSize_; }

    /// Update control factor (ie time)
    void updateTime() { this -> updateTime( timeStepSize_ ); }

    /// Compute and store (negative) force residual (ie RHS) on nodes
    virtual void computeForces();

    /// Compute and assemble force residual
    virtual void computeAndAssembleForces()
    {
        this -> computeForces();

        this -> assembleForces();
    }

    /// Update displacements by adding displacement increment of step
    virtual void updateSolution()
    {
        this -> updateDisplacements();
    }

    /// Write step (ie VTU file)
    ///
    /// \param[in]  time  Time or characteristic load
    /// \param[in]  step  Time step or load step index
    virtual void writeAnimationStep( const double factor, const unsigned step )
    {
        if ( animation_ ) {
            const std::string vtuFile = animation_ -> snapshotName( factor, step );
            this -> writeMeshData( vtuFile );
        }
    }

    /// Write animation file (ie PVD file)
    void writeAnimationFile( const std::string animationFileName ) const
    {
        if ( animation_ ) {
            const std::string animationFile( animationFileName.c_str() );
            std::ofstream anim( animationFile.c_str( ) );
            animation_ -> writeAnimationFile( anim );
            anim.close();
        }
    }

    /// Define what node data is written to VTU
    virtual void writeNodeData( VTUwriter & vtuwriter ) const;

    /// Define what element data is written to VTU
    virtual void writeElementData( VTUwriter & vtuwriter ) const;

    //@}

    //--------------------------------------------------------------------------
    /// @name Fine level methods
    //@{

    /// number the DOFs
    virtual unsigned numberDofs()
    {
        numDofs_ = mesh_.iterateOverNodes( corlib::dofNumberFun( &Node::numberDOFs, 0 ) );
        return numDofs_;
    }

    /// Clear constrained DOFs on nodes
    void clearNodeConstraints()
    {
        mesh_.iterateOverNodes( boost::bind( &Node::clearConstraints, _1 ) );
    }

    /// Blank forces (distributed RHS) on nodes
    void clearForces()
    {
        mesh_.iterateOverNodes( boost::bind( &Node::clearForce, _1 ) );
    }

    /// Compute (ie multiply with factor) external forces due to concentrated loads
    void computeNodalLoads( const double loadFactor )
    {
        mesh_.iterateOverNodes( boost::bind( nodalForces_, _1, loadFactor ) );
    }

    /// Compute external forces due to body load
    void computeBodyLoad( VecDimToVecDof forceFun, const double factor,
                          ElementPred condition = corlib::Positive() )
    {
        BodyForce bf( forceFun );
        bf.setFactor( factor );
        corlib::Integrator<Quadrature,BodyForce> bfIntegrator( quadrature_, bf );
        mesh_.iterateOverElementsWithPredicate( bfIntegrator, condition );
    }

    /// Compute internal forces
    void computeInternalForces()
    {
        Integrand internalForceIntegrand( &Element::internalForceIntegrand );
        Integrator internalForceIntegrator( quadrature_, internalForceIntegrand );
        mesh_.iterateOverElements( internalForceIntegrator );
    }

    /// Create linear algebraic solver
    virtual void allocateSolver()
    {
        if ( numDofs_ == 0 ) {
            numDofs_ = this -> numberDofs();
            this -> accessOutStream_()
                << numDofs_ << " degrees of freedom" << std::endl;
        }
        this -> allocateSolver( numDofs_ );
    }

    virtual void allocateSolver( const unsigned size )
    {
        detail_::allocateSolver( solver_, size );
    }

    /// Destroy linear algebraic solver
    void deAllocateSolver()
    {
        delete solver_;  solver_ = NULL;
    }

    /// Clear global system of equations: matrix and RHS
    virtual void clearSystem()
    {
        solver_->clearRhs();
        solver_->clearMatrix();
    }

    /// Assemble element stiffness matrices in global matrix
    virtual void computeAndAssembleStiffness()
    {
        //! Compute and assemble element stiffness matrices
        corlib::MatrixComputeAndAssembleFun<const Element,Quadrature,Solver>
            mafStiffness( &Element::stiffnessIntegrand, 
                          &Element::getDofIndices,
                          &Element::getDofIndices,
                          quadrature_, solver_ );
        mesh_.iterateOverElements( mafStiffness );
    }

    /// Assemble nodally distributed residual forces into global RHS vector
    void assembleForces( )
    {
        mesh_.iterateOverNodes( corlib::vectorAssemblerFun( &Node::giveForce,
                                                            &Node::copyDofArray,
                                                            solver_ ) );
    }

    /// Finalise assembly of global system
    virtual void finishAssembly()
    {
        solver_ -> finishAssembly();
    }

    /// Introduce constrained (dropped) DOFs in linearised system
    void applyConstraints( const double factor )
    {
        corlib::ConstraintFunctor<Node,Solver> constraint( solver_, factor );
        constraint = mesh_.iterateOverNodes( constraint );
        constraint.applyConstraints();
    }

    /// Solve linear system of equations arising from discretisation and linearisation
    void solveSystem()
    {
        solver_ -> solveSystem();
    }

    /// Inform user why/how linear system solve did (may be iterative solver)
    virtual void informAboutSolve() const { } // do nothing

    /// Return norm on right-hand-side (out-of-balance-force)
    double normRhs() const
    {
        return solver_ -> normRhs();
    }

    /// Return norm on residual displacements
    double normSolution() const
    {
        return solver_ -> normSol();
    }

    /// Iterative update of solution increment stored on nodes
    void distributeSolution()
    {
        mesh_.iterateOverNodes( corlib::distributorFun( solver_,
                                                        &Node::copyDofArray,
                                                        &Node::addToIncrement ) );
    }

    /// Update displacements by adding displacement increment of step
    void updateDisplacements()
    {
        mesh_.iterateOverNodes( boost::bind( &Node::updateDisplacements, _1 ) );
    }

    /// Update control factor (ie time)
    void updateTime( const double stepSize )
    {
        time_ += stepSize;
    }

    /// Write VTU file
    virtual void writeMeshData( const std::string & vtuFile );

    //@}

    //--------------------------------------------------------------------------
    /// @name Accessors
    //@{

    /// Grant access to mechanical body
    Mesh & accessMesh() { return mesh_; }

    /// Provide pointer to solver
    Solver * getSolver() { return solver_; }

    //@}

protected:
    /// Access to output stream
    ///
    /// This is replaced with non std::cout in parallel environment
    virtual std::ostream & accessOutStream_() const { return std::cout; }

protected:
    /// Solid mesh
    Mesh               mesh_;

    /// Constraints scaling factor function
    DoubleToDouble     constraintsFac_;

    /// Concentrated nodal loads
    NodalForces        nodalForces_;
    /// Concentracted nodal loads factor function
    DoubleToDouble     nodalForcesFac_;
    /// Body force function
    VecDimToVecDof     bodyForceFun_;
    /// Body force load factor function
    DoubleToDouble     bodyForceFac_;
    /// Body force conditional
    ElementPred        bodyForcePred_;

    //! Quadrature object
    Quadrature         quadrature_;

    /// Number of degrees of freedeom
    unsigned           numDofs_;
    /// Linear algebraic solver
    Solver *           solver_;

    /// Load/displacement control factor
    double             time_;
    /// Increments of load/displacement control factor
    double             timeStepSize_;
    
    /// Output
    VTUanim *          animation_;

};

//------------------------------------------------------------------------------
#include "Static.ipp"

#endif
