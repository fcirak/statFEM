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
#include <iostream>
#include <fstream>
#include <string>
#include <cstdlib>
#include <functional>

#include <boost/function.hpp>

#include <corlib/Mesh.hpp>
#include <corlib/Shapefun.hpp>
#include <corlib/Quadrature.hpp>
#include <corlib/DofNumberer.hpp>
#include <corlib/SystemSolveEigenSparse.hpp>
#include <corlib/Distributor.hpp>
#include <corlib/Integrator.hpp>
#include <corlib/Accessor.hpp>
#include <corlib/MatrixAssembler.hpp>
#include <corlib/VectorAssembler.hpp>
#include <corlib/BodyForce.hpp>
#include <corlib/Constraints.hpp>
#include <corlib/PropertiesParser.hpp>
#include <corlib/VTUwriter.hpp>
#include <corlib/linalg.hpp>
#include <corlib/fuzzyEqual.hpp>
#include <corlib/verify.hpp>
#include <corlib/ComputeElementMatrix.hpp>
#include <corlib/NodeBasic.hpp>
#include <corlib/ElementBasic.hpp>
#include <corlib/NormError.hpp>
#include <corlib/ShapefunTraits.hpp>
#include <corlib/eigenX.hpp>

#include <del2/fem/NodePotential.hpp>
#include <del2/fem/ElementPotential.hpp>

//#define NOWATCHDOG
namespace eigenX = corlib::eigenX;

//------------------------------------------------------------------------------
template<unsigned DIM>
eigenX::VectorSd<1>
analytical( const eigenX::VectorSd<DIM> X )
{
    eigenX::VectorSd<1> u;

    //const double r2 = (X[0] - 1.1)*(X[0]-1.1) + (X[1] - 1.1)*(X[1] - 1.1);
    //u[0] = std::log( r2 );
    //u[0] = - std::exp( X[0] ) + (std::exp(1.)-1)*X[0] + 1.;
    //u[0] = 0.5*X[0]*(1.-X[0]);

    u[0] = std::sin( 6. * X[0] ) * std::sin( 8. * X[1] );
    return u;
}

//------------------------------------------------------------------------------
template<unsigned DIM>
eigenX::VectorSd<1>
force( const eigenX::VectorSd<DIM> X )
{
    eigenX::VectorSd<1> F = 100. * analytical<DIM>(X);
    //F[0] = 0.;
    //F[0] = std::exp( X[0] );
    //F[0] = 1.;
    return F;
}

//------------------------------------------------------------------------------
// homogeneous Dirichlet boundary conditions
template<unsigned DIM>
corlib::NodalConstraint
dirichletFun( const eigenX::VectorSd<DIM> X )
{
    corlib::NodalConstraint constraint;
    if ( corlib::fuzzyEqual( X[0], 0. ) or
         corlib::fuzzyEqual( X[1], 0. ) or
         corlib::fuzzyEqual( X[0], 1. ) or
         corlib::fuzzyEqual( X[1], 1. ) ) {

        eigenX::VectorSd<1> U = analytical<DIM>(X);
        constraint.setComponent( 0, U[0] );
    }
    return constraint;
}


//------------------------------------------------------------------------------
template<typename MESH, typename FUNC>
void writeMeshData( const std::string & vtuFile, MESH * mesh, FUNC func );

int main( int argc, char* argv[] )
{
    //! variables to be read from input file "./input.dat"
    std::string meshFile;
    double conductivity;

    //! instantiate and feed properties parser with the variables to read
    corlib::PropertiesParser *prop = new corlib::PropertiesParser;
    prop -> registerPropertiesVar( "meshFile",         meshFile   );
    prop -> registerPropertiesVar( "conductivity",     conductivity );

    //! read variables from the input.dat file
    std::ifstream inputFile( "./input.dat" );
    FTL_VERIFY( inputFile.is_open() );
    prop -> readValues( inputFile );

    //! delete properties parser and close the input file
    delete prop;
    inputFile.close( );

    //! indicate if we want to compute with triangles or quadrilaterals
    //! this has to, of course, correlate with the element type in meshfile
    const bool useTriangles = false;
    const corlib::shape elemShape = ( useTriangles ? corlib::TRIANGLE : corlib::QUADRILATERAL );

    //! choosing TRIANGLE or QUADRILATERAL determines the space dimension
    const unsigned  dim  = corlib::ShapeTraits<elemShape>::dim;

    //! instantiate quadratic TENSORLAGRANGE basis functions
    const unsigned degree = 2;
    const corlib::SfunClass sfunClass = ( useTriangles ? corlib::SIMPLEX : corlib::TENSORLAGRANGE );
    typedef corlib::ShapeFunType<degree,dim,sfunClass>::type Sfun;

    //! instantiate a quadrature object
    const unsigned numQuadPoints = ( useTriangles ? 6 : 9 );
    typedef corlib::Quadrature<elemShape,numQuadPoints>  Quad;

    //! instantiate the data structures for computing
    typedef corlib::NodeBasic<dim>                                           BasisNode;
    typedef del2::fem::NodePotential<BasisNode>                              Node;
    typedef corlib::ElementBasic<Node,Sfun>                                  BasisElement;
    typedef del2::fem::ElementPotential<BasisElement>                        Element;
    typedef corlib::Mesh<Element>                                            Mesh;
    typedef corlib::SystemSolveEigenSparse                                   SysSolve;

    //! read smf file and construct mesh
    std::ifstream smf( meshFile.c_str( ) );
    FTL_VERIFY( smf.is_open() );
    Mesh mesh( smf );
    smf.close( );

    //! cache material to elements
    mesh.iterateOverElements( boost::bind( &Element::setConductivity, _1, conductivity ) );

    //! Apply dirichlet boundary conditions
    typedef eigenX::VectorSd<dim>                                        VecDim;
    typedef boost::function< corlib::NodalConstraint( const VecDim ) >   VecDim2Constraint;
    typedef corlib::ConstraintsFromFunction<Node,VecDim2Constraint>      Dirichlet;
    Dirichlet dirichlet( boost::bind( dirichletFun<dim>, _1 ) );
    mesh.iterateOverNodes( dirichlet );

    //! number dofs
     unsigned numdofs = 0;
     mesh.iterateOverNodes( boost::bind( &Node::numberDOFs, _1, boost::ref( numdofs ) ) );
     std::cout << numdofs << " degrees of freedom" << std::endl;

    // alternative approach for dof numbering
    //  unsigned numdofs = mesh.iterateOverNodes( corlib::dofNumberFun( &Node::numberDOFs, 0 ) );

    //! system matrix object
    SysSolve * sysmat  = new SysSolve( numdofs );

    //! instantiate a quadrature object
    Quad quadrature;

    //! compute element stiffness matrices
    corlib::MatrixComputeAndAssembleFun<const Element,Quad,SysSolve>
        maf( &Element::stiffnessIntegrand,
             &Element::getDofIndices,
             &Element::getDofIndices,
             quadrature, sysmat );
    mesh.iterateOverElements( maf );

    //! Clear force arrays
    mesh.iterateOverNodes( boost::bind( &Node::clearForce, _1 ) );

    //! body force
    typedef eigenX::VectorSd<1>                                  Vec1;
    typedef boost::function< Vec1( const VecDim ) >              VecDim2Vec1;
    typedef corlib::BodyForce<Element,VecDim2Vec1>               BodyForce;
    BodyForce bodyForce( boost::bind( force<dim>, _1 ) );
    corlib::Integrator<Quad,BodyForce>  bodyForceIntegrator( quadrature, bodyForce );
    mesh.iterateOverElements( bodyForceIntegrator );

    //! assemble the  forces to the system matrices rhs
    mesh.iterateOverNodes( corlib::vectorAssemblerFun( &Node::getForce,
                                                       &Node::copyDofArray,
                                                       sysmat ) );

    //! collect and apply constraints
    sysmat -> finishAssembly();
    corlib::ConstraintFunctor< Node, SysSolve > constraint( sysmat );
    constraint = mesh.iterateOverNodes( constraint );
    constraint.applyConstraints();

    //! solve the system of equations
    sysmat -> solveSystem( );

    //! add nodal solutions to increment
    mesh.iterateOverNodes( corlib::distributorFun( sysmat,
                                                   &Node::copyDofArray,
                                                   &Node::setIncrement ) );

    //! clear the stiffness matrix and RHS
    sysmat -> clearMatrix( );
    sysmat -> clearRhs( );

    //! set potential equal to the increment
    mesh.iterateOverNodes( boost::bind( &Node::updatePotential, _1 ) );

    //! L2 error
    typedef boost::function<Vec1( const Element*, const VecDim &)> Approximate;
    typedef corlib::NormError<VecDim2Vec1, Approximate>         L2Error;
    Vec1 el2ErrorSquared; el2ErrorSquared.setZero( );
    L2Error el2Error( boost::bind( analytical<dim>, _1 ),
                      boost::bind( &Element::potential, _1, _2 ), el2ErrorSquared );
    corlib::Integrator<Quad,L2Error> errorIntegratorEl2( quadrature, el2Error );
    mesh.iterateOverElements( errorIntegratorEl2 );

    const double l2Error = std::sqrt( el2ErrorSquared[0] );

    std::cout << "L2-error: " << l2Error << std::endl;

    std::ofstream err( "error.dat" );
    err << "L2-error = " << l2Error << std::endl;
    err.close();

    //! write data to VTU file
    std::string vtuFile = "laplace.vtu";
    writeMeshData<Mesh,VecDim2Vec1>( vtuFile, &mesh, boost::bind( force<dim>, _1 ) );


    delete sysmat;

    return 0;
}

//------------------------------------------------------------------------------
template< typename MESH, typename FUNC>
void writeMeshData( const std::string & vtuFile, MESH * mesh, FUNC func )
{
    // set up the VTU writer and open a file
    std::ofstream vtu( vtuFile.c_str( ) );
    assert( vtu.is_open( ) );
    corlib::VTUwriter<MESH> vtuwriter( mesh, vtu );

    //! write mesh data
    vtuwriter.writeMesh();

    //-- write point data
    vtuwriter.openPointData();
    vtuwriter.writeNodalQuantity( corlib::accessorFun( &MESH::Node::getPotential, "Potential" ) );
    vtuwriter.closePointData( );

    //-- write cell data
#ifdef NOWATCHDOG
    vtuwriter.openCellData();

    typedef boost::function< typename MESH::Node::VecDof( const typename MESH::Element* ) > GetVecDof;

    GetVecDof elemResidual = boost::bind( &MESH::Element::template equationResidual<FUNC>, _1,
                                          corlib::ShapeTraits<MESH::Element::myShape>::centroid(), 
                                          func );

    corlib::Accessor<GetVecDof> accessResidual( elemResidual, "Residual" );
    vtuwriter.writeElementQuantity( accessResidual );

    typedef boost::function< typename MESH::Element::Vec3( const typename MESH::Element* ) >
        Elem2Vec3;
    typename MESH::Element::VecLDim xi; xi.setZero();
    Elem2Vec3 flux = boost::bind( &MESH::Element::internalFlux, _1, 
                                  corlib::ShapeTraits<MESH::Element::myShape>::centroid() );
    corlib::Accessor<Elem2Vec3> accessFlux( flux, "Flux" );
    vtuwriter.writeElementQuantity( accessFlux );

    vtuwriter.closeCellData( );
#endif


    // finish writing
    vtuwriter.finalize();
    vtu.close( );
}
