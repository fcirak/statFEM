// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//
//                 Computational Structural Mechanics Lab
//                         University of Cambridge
//
// This software is copyrighted and all rights are retained by the University of
// Cambridge. All use, disclosure, and/or reproduction of any part not expressly
// authorised by the University of Cambridge is prohibited. (C) 2010
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

//! @file Dynamic.hpp

#ifndef solid_driver_dynamic_h
#define solid_driver_dynamic_h

// includes from corlib
#include <corlib/BackwardEulerAssembler.hpp>
#include <corlib/NewmarkAssembler.hpp>

// includes from solid::driver
#include <solid/driver/Static.hpp>

//==============================================================================
namespace solid {
    namespace driver {

        namespace ublas = boost::numeric::ublas;

        template< typename ELEMENT, typename SOLVER,
                  unsigned NGPSTIFF, unsigned NGPMASS > class Dynamic;

        template< typename ELEMENT, typename SOLVER,
                  unsigned NGPSTIFF, unsigned NGPMASS > class Newmark;

        template< typename ELEMENT, typename SOLVER,
                  unsigned NGPSTIFF, unsigned NGPMASS > class BackwardEuler;

    }
}

//------------------------------------------------------------------------------
/// Base class driver for dynamic solid computations
///
/// \tparam ELEMENT  Solid element
/// \tparam SOLVER   Linear algebraic solver
/// \tparam NGPSTIFF Number of Gauss quadrature points per element
///                  (defaults to 0: the number is deduced from traits)
/// \tparam NGPMASS  Number of Gauss quadrature points per element for mass matrix
///                  (defaults to 0: the number is deduced from traits)
template< typename ELEMENT, typename SOLVER,
          unsigned NGPSTIFF = 0, unsigned NGPMASS = 0 >
class solid::driver::Dynamic
    : public solid::driver::Static< ELEMENT, SOLVER, NGPSTIFF >
{
protected:
    typedef solid::driver::Static< ELEMENT, SOLVER,
                                   NGPSTIFF >                 DriverStatic_;

public:
    typedef typename DriverStatic_::Element                   Element;
    typedef typename DriverStatic_::Node                      Node;

    typedef typename DriverStatic_::VecDim                    VecDim;
    typedef typename Element::VecLDim                         VecLDim;

    typedef typename DriverStatic_::Mesh                      Mesh;

    static const enum corlib::shape myShape = DriverStatic_::myShape;
    static const unsigned numQuadPointsMass = NGPMASS;
    BOOST_STATIC_ASSERT( numQuadPointsMass > 0 );
    typedef corlib::Quadrature< myShape, numQuadPointsMass >  QuadratureMass;
    typedef boost::function<void( const Element *, 
                                  const VecLDim &, 
                                  const double &, 
                                  ublas::matrix<double> & ) > Integrand2;
    typedef corlib::ComputeElementMatrix< QuadratureMass,
                                          Integrand2 >        MatrixComputerMass;

    typedef typename DriverStatic_::Quadrature                QuadratureStiff;
    typedef corlib::ComputeElementMatrix< QuadratureStiff,
                                          Integrand2 >        MatrixComputerStiff;

    typedef typename DriverStatic_::Solver                    Solver;

    typedef typename DriverStatic_::VTUwriter                 VTUwriter;

public:
    //--------------------------------------------------------------------------
    /// @name Coarse level methods
    //@{

    /// Default constructor
    Dynamic( std::istream & smf ) : DriverStatic_( smf ) { }

    //@}

    //--------------------------------------------------------------------------
    /// @name Fine level methods
    //@{

    /// Write VTU file containing also dynamical data
    virtual void writeNodeData( VTUwriter & vtuwriter ) const
    {
        vtuwriter.writeNodalQuantity( corlib::accessorFun( &Node::giveDisplacements,
                                                           "Displacements" ) );
        vtuwriter.writeNodalQuantity( corlib::accessorFun( &Node::giveVelocities,
                                                           "Velocities" ) );
        //vtuwriter.writeNodalQuantity( corlib::accessorFun( &Node::giveAccelerations,
        //                                                   "Accelerations" ) );
    }

    //@}

protected:
    /// Shell mesh
    using DriverStatic_::mesh_;
    /// Linear algebraic solver
    using DriverStatic_::solver_;
    ///
    using DriverStatic_::quadrature_;
    ///
    QuadratureMass quadratureMass_;

    /// control parameters
    using DriverStatic_::timeStepSize_;
};


//==============================================================================
/// Driver for dynamic solid computations using Newmark's time stepping
/// \tparam ELEMENT  Solid element
/// \tparam SOLVER   Linear algebraic solver
/// \tparam NGPSTIFF Number of Gauss quadrature points per element
///                  (defaults to 0: the number is deduced from traits)
/// \tparam NGPMASS  Number of Gauss quadrature points per element for mass matrix
///                  (defaults to 0: the number is deduced from traits)
template< typename ELEMENT, typename SOLVER,
          unsigned NGPSTIFF = 0, unsigned NGPMASS = 0 >
class solid::driver::Newmark
    : public solid::driver::Dynamic< ELEMENT, SOLVER, NGPSTIFF, NGPMASS >
{
protected:
    typedef solid::driver::Dynamic< ELEMENT, SOLVER,
                                    NGPSTIFF, NGPMASS >       DriverDynamic_;

public:
    typedef typename DriverDynamic_::Element                  Element;
    typedef typename DriverDynamic_::Node                     Node;

    typedef typename DriverDynamic_::Mesh                     Mesh;

    typedef typename DriverDynamic_::Integrand2               Integrand2;
    typedef typename DriverDynamic_::MatrixComputerMass       MatrixComputerMass;
    typedef typename DriverDynamic_::MatrixComputerStiff      MatrixComputerStiff;

    typedef typename DriverDynamic_::Solver                   Solver;

public:
    //--------------------------------------------------------------------------
    /// @name Coarse level methods
    //@{

    /// Default constructor
    Newmark( std::istream & smf, const double beta, const double gamma )
        : DriverDynamic_( smf ),
          beta_( beta ),
          gamma_( gamma )
    { }

    //@}

    //--------------------------------------------------------------------------
    /// @name Medium level methods
    //@{

    /// Compute and store dynamic force residual (ie negative RHS) on nodes
    /// followed by assembling them into global residual (global RHS)
    virtual void computeAndAssembleForces()
    {
        this -> DriverDynamic_::computeAndAssembleForces();

        this -> computeForcesNewmark( beta_, gamma_, timeStepSize_ );
    }

    /// Compute dynamic element stiffness (ie LHS is RHS differentiated WRT 
    /// displacements) and assemble them into global matrix (global LHS)
    virtual void computeAndAssembleStiffness()
    {
        this -> computeAndAssembleMatricesNewmark( beta_, gamma_, timeStepSize_ );
    }

    /// Update displacements, velocities and accelerations
    virtual void updateSolution()
    {
        this -> updateNewmark( beta_, gamma_, timeStepSize_ );

        this -> DriverDynamic_::updateSolution();
    }

    //@}

    //--------------------------------------------------------------------------
    /// @name Fine level methods
    //@{
    
    /// Compute complete matrix (LHS) containing stiffness and mass matrices
    void computeAndAssembleMatricesNewmark( const double beta,
                                            const double gamma,
                                            const double stepSize )
    {
        Integrand2 massIntegrand( &Element::massIntegrand );
        Integrand2 stiffnessIntegrand( &Element::stiffnessIntegrand );
    
        MatrixComputerMass mass( quadratureMass_, massIntegrand );
        MatrixComputerStiff stiff( quadrature_, stiffnessIntegrand );

        corlib::NewmarkMatrixAssembler<Solver,Element,MatrixComputerMass,MatrixComputerStiff>
            nma( beta, gamma, stepSize, mass, stiff );

        mesh_.iterateOverElements( boost::bind( nma, _1, solver_ ) );
    }

    /// Compute inertia forces and store them on nodes
    void computeForcesNewmark( const double beta,
                               const double gamma,
                               const double stepSize )
    {
        Integrand2 massIntegrand( &Element::massIntegrand );
        MatrixComputerMass mass( quadratureMass_, massIntegrand );
        corlib::NewmarkInertiaForceAssembler<Solver,Element,MatrixComputerMass>
            nifa( beta, gamma, stepSize, mass );
        mesh_.iterateOverElements( boost::bind( nifa, _1, solver_ ) );
    }

    /// Update nodal velocities and accelerations
    void updateNewmark( const double beta,
                        const double gamma,
                        const double stepSize )
    {
        corlib::NewmarkSolutionDistributor<Node> nsd( beta, gamma, stepSize );
        mesh_.iterateOverNodes( nsd );
    }

    //@}

protected:
    /// Shell mesh
    using DriverDynamic_::mesh_;
    /// Linear algebraic solver
    using DriverDynamic_::solver_;
    ///
    using DriverDynamic_::quadrature_;
    ///
    using DriverDynamic_::quadratureMass_;

    /// control parameters
    using DriverDynamic_::timeStepSize_;
    double         beta_;
    double         gamma_;
};

//==============================================================================
/// Driver for dynamic solid computations using backward-Euler time stepping
/// \tparam ELEMENT  Solid element
/// \tparam SOLVER   Linear algebraic solver
/// \tparam NGPSTIFF Number of Gauss quadrature points per element
///                  (defaults to 0: the number is deduced from traits)
/// \tparam NGPMASS  Number of Gauss quadrature points per element for mass matrix
///                  (defaults to 0: the number is deduced from traits)
template< typename ELEMENT, typename SOLVER,
          unsigned NGPSTIFF = 0, unsigned NGPMASS = 0 >
class solid::driver::BackwardEuler
    : public solid::driver::Dynamic< ELEMENT, SOLVER, NGPSTIFF, NGPMASS >
{
protected:
    typedef solid::driver::Dynamic< ELEMENT, SOLVER,
                                    NGPSTIFF, NGPMASS >       DriverDynamic_;

public:
    typedef typename DriverDynamic_::Element                  Element;
    typedef typename DriverDynamic_::Node                     Node;

    typedef typename DriverDynamic_::Mesh                     Mesh;

    typedef typename DriverDynamic_::Integrand2               Integrand2;
    typedef typename DriverDynamic_::MatrixComputerMass       MatrixComputerMass;
    typedef typename DriverDynamic_::MatrixComputerStiff      MatrixComputerStiff;

    typedef typename DriverDynamic_::Solver                   Solver;

public:
    //--------------------------------------------------------------------------
    /// @name Coarse level methods
    //@{

    /// Default constructor
    BackwardEuler( std::istream & smf ) : DriverDynamic_( smf ) { }

    //@}

    //--------------------------------------------------------------------------
    /// @name Medium level methods
    //@{

    /// Compute and store dynamic force residual (ie negative RHS) on nodes
    /// followed by assembling them into global residual (global RHS)
    virtual void computeAndAssembleForces()
    {
        this -> DriverDynamic_::computeAndAssembleForces();

        this -> computeForcesBE( timeStepSize_ );
    }

    /// Compute dynamic element stiffness (ie LHS is RHS differentiated WRT 
    /// displacements) and assemble them into global matrix (global LHS)
    virtual void computeAndAssembleStiffness()
    {
        this -> computeAndAssembleMatricesBE( timeStepSize_ );
    }

    /// Update displacements, velocities and accelerations
    virtual void updateSolution()
    {
        this -> updateBE( timeStepSize_ );

        this -> DriverDynamic_::updateSolution();
    }

    //@}

    //--------------------------------------------------------------------------
    /// @name Fine level methods
    //@{
    
    /// Compute complete matrix (LHS) containing stiffness and mass matrices
    void computeAndAssembleMatricesBE( const double stepSize )
    {
        Integrand2 massIntegrand( &Element::massIntegrand );
        Integrand2 stiffnessIntegrand( &Element::stiffnessIntegrand );
    
        MatrixComputerMass mass( quadratureMass_, massIntegrand );
        MatrixComputerStiff stiff( quadrature_, stiffnessIntegrand );

        corlib::BackwardEulerMatrixAssembler<Solver,Element,MatrixComputerMass,MatrixComputerStiff>
            bema( stepSize, mass, stiff );

        mesh_.iterateOverElements( boost::bind( bema, _1, solver_ ) );
    }

    /// Compute inertia forces and store them on nodes
    void computeForcesBE( const double stepSize )
    {
        Integrand2 massIntegrand( &Element::massIntegrand );
        MatrixComputerMass mass( quadratureMass_, massIntegrand );
        corlib::BackwardEulerInertiaForceAssembler<Solver,Element,MatrixComputerMass>
            beifa( stepSize, mass );
        mesh_.iterateOverElements( boost::bind( beifa, _1, solver_ ) );
    }

    /// Update nodal velocities and accelerations
    void updateBE( const double stepSize )
    {
        corlib::BackwardEulerSolutionDistributor<Node> besd( stepSize );
        mesh_.iterateOverNodes( besd );
    }

    //@}

protected:
    /// Shell mesh
    using DriverDynamic_::mesh_;
    /// Linear algebraic solver
    using DriverDynamic_::solver_;
    /// Quadrature for stiffness
    using DriverDynamic_::quadrature_;
    /// Quadrature for mass
    using DriverDynamic_::quadratureMass_;

    /// control parameters
    using DriverDynamic_::timeStepSize_;
};

#endif
