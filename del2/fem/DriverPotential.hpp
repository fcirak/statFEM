// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//
//                 Computational Structural Mechanics Lab
//                         University of Cambridge
//
// This software is copyrighted and all rights are retained by the University of
// Cambridge. All use, disclosure, and/or reproduction of any part not expressly
// authorised by the University of Cambridge is prohibited. (C) 2010
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

//! @file DriverPotential.hpp

#ifndef del2_fem_driverpotential_h
#define del2_fem_driverpotential_h

//------------------------------------------------------------------------------
// headers
#include <iostream>
#include <string>
#include <fstream>

#include <boost/bind.hpp>

#include <Eigen/Core>
#include <corlib/eigenX.hpp>

#include <corlib/Shape.hpp>
#include <corlib/verify.hpp>

#include <corlib/PropertiesParser.hpp>
#include <corlib/Quadrature.hpp>
#include <corlib/Constraints.hpp>
#include <corlib/BodyForce.hpp>
#include <corlib/Integrator.hpp>
#include <corlib/MatrixAssembler.hpp>
#include <corlib/VectorAssembler.hpp>
#include <corlib/Distributor.hpp>
#include <corlib/VTUwriter.hpp>
#include <corlib/DofNumberer.hpp>

#include <del2/fem/PotentialQuadrature.hpp>

//------------------------------------------------------------------------------
// forward declaration
namespace del2 {
    namespace fem {

        //======================================================================
        template< typename MESH, typename SYSTEMSOLVER >
        class DriverPotential;
        namespace eigenX = corlib::eigenX;
    }
}

//==============================================================================
//! A driver class to provide covenience functions to perform a numerical
//! solution of the Laplace (or Poisson) equation.
//! 
//! \tparam  MESH          Finite element mesh (containing nodes and elements)
//! \tparam  SYSTEMSOLVER  Linear system and its solver
template< typename MESH, typename SYSTEMSOLVER >
class del2::fem::DriverPotential
{
public:
    typedef MESH                                                   Mesh;
    typedef SYSTEMSOLVER                                           SysSolve;

public:
    typedef typename Mesh::Node                                    Node;
    typedef typename Node::VecDim                                  VecDim;

    typedef typename Mesh::Element                                 Element;

public:
    static const unsigned dim           = Node::dim;
    static const unsigned dof           = Node::dof;
    static const corlib::shape myShape  = Element::myShape;
    static const unsigned numNodes      = Element::numNodes;

    static const unsigned numQuadPoints = PotentialQuadratureTraits< myShape, numNodes >::numPoints;

public:
    //! Quadrature rule used for integration on element domain
    typedef corlib::Quadrature< myShape, numQuadPoints >            Quad;

    typedef corlib::NodalConstraint                                 Constraint;
    typedef boost::function< Constraint( const VecDim ) >           DirichFun;
    typedef corlib::ConstraintsFromFunction< Node, DirichFun >      Dirichlet;

    typedef eigenX::VectorSd<dof>                                   VecDof;
    typedef boost::function< VecDof( const VecDim ) >               BodyForceFun;
    typedef corlib::BodyForce< Element, BodyForceFun >              BodyForce;

    typedef boost::function< void( Element*, const VecDim &,
                                   const double & ) >               Integrand;

public:
    //! Constructor
    //!
    //! \param[in]  inputFileName  input file name
    DriverPotential( const std::string & inputFileName );

    DriverPotential() :
        fileBaseName_( "" ),
        conductivity_( 0.0 ),
        mesh_( NULL ),
        dirichlet_( NULL ),
        bodyForce_( NULL ),
        numDofs_( 0 ),
        sysMat_( NULL ){ }

    //! Destructor
    virtual ~DriverPotential( );

    //! Set Dirichlet boundary conditions
    virtual void createDirichletConstraints( DirichFun dirichFun );

    //! Set body force function
    virtual void createBodyForce( BodyForceFun bodyForceFun );

    //! Number DOFs of Potential field
    //!
    //! \param[in]  seedDof  starting DOF index
    //! \return              Number of DOFs living on Potential mesh
    virtual unsigned numberDofs( const unsigned seedDof );

    //! Create/allocate global system of linear equations
    virtual void createSystem( );

    //! Destroy/deallocate global system of linear equations
    virtual void destroySystem( );

    //! Set global matrix and vectors to zero.
    //!
    //! Useful if several steps are computed and the pattern of the matrix
    //! can be kept.
    virtual void clearSystem( );

    //! Create residual
    virtual void assembleResidual( );

    //! Assemble tangent
    virtual void assembleJacobian( );

    //! Finish assembly
    virtual void assembleFinish( );

    //! Apply constraints
    virtual void dirichletConstrainSystem( );

    //! Solve system of linear equations
    virtual void solveSystem( );

    //! Update quantities after iterative solution
    virtual void updateIncrement( );

    //! Update quantities after iteration has converged, i.e.
    //! update at load or time step
    virtual void updateStep( );

    //! Provide base name of application
    virtual std::string fileBaseName( ) const { return fileBaseName_; }

    //! Write step to VTU file
    //!
    //! \param[in] vtuFileName  The file name to which the VTU data is written.
    virtual void writeVtu( const std::string & vtuFileName );

protected:
    //! file base name
    std::string               fileBaseName_;

    //! Conductivity material constant
    double                    conductivity_;

    //! The mesh (containing elements and nodes)
    Mesh *                    mesh_;

    //! Dirichlet constraints
    Dirichlet *               dirichlet_;

    //! Body force
    BodyForce *               bodyForce_;

    //! Number of degree-of-freedom of Potential field
    unsigned                  numDofs_;

    //! Quadrature object for numerical integration
    Quad                      quadrature_;

    //! Solver for global system of linear equations
    SysSolve *                sysMat_;
};

//------------------------------------------------------------------------------
// implementations
#include "DriverPotential.ipp"

#endif
