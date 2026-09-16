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
/** Compute potential flow around some object
 *      
 *                   p_n=0
 *         +---------------------+
 *         |     ___             |
 *         |    /   \            |
 *  p_n=-1 |    |   _\ <- Object | p_n=+1
 *         |    \__/             |
 *         |         p_n=0       |
 *         +---------------------+
 *        x=0                  x=x_max > 1.
 */
//------------------------------------------------------------------------------
#include <iostream>
#include <fstream>
#include <string>

#include <boost/function.hpp>

#include <Eigen/Core>
#include <corlib/eigenX.hpp>

#include <corlib/Mesh.hpp>
#include <corlib/Shapefun.hpp>
#include <corlib/Quadrature.hpp>
#include <corlib/DofNumberer.hpp>
#include <corlib/SystemSolveEigenSparse.hpp>
#include <corlib/Distributor.hpp>
#include <corlib/Accessor.hpp>
#include <corlib/MatrixAssembler.hpp>
#include <corlib/VectorAssembler.hpp>
#include <corlib/BodyForce.hpp>
#include <corlib/Constraints.hpp>
#include <corlib/PropertiesParser.hpp>
#include <corlib/VTUwriter.hpp>
#include <corlib/verify.hpp>
#include <corlib/MeshBoundary.hpp>
#include <corlib/NodeBasic.hpp>
#include <corlib/ElementBasic.hpp>
#include <corlib/Integrator.hpp>

#include <del2/fem/NodePotential.hpp>
#include <del2/fem/ElementPotential.hpp>

namespace eigenX = corlib::eigenX;


//------------------------------------------------------------------------------
//! Set Neumann datum to +/- 1 at the right/left extremal boundaries
template<unsigned DIM>
eigenX::VectorSd<1>
normalDerivative( const eigenX::VectorSd<DIM> & X )
{

    eigenX::VectorSd<1> res;
    if ( X[0] > 1.0 ) res[0] =  1.; // see picture above
    else              res[0] = -1.;
    return res;
}

//------------------------------------------------------------------------------
//! Helper to write the result file
template<typename MESH>
void writeMeshData( const std::string & vtuFile, MESH * mesh );


//------------------------------------------------------------------------------
int main( int argc, char* argv[] )
{ 
    //! Variable attributes (depending on element shape)
    const corlib::shape elemShape  = corlib::TRIANGLE;
    const unsigned      numNodesPE = 3;
    const unsigned numGaussPoints  = 3;
    const unsigned numGaussSurf    = 2;
    
    //! Fixed attributes
    const unsigned      dim        = corlib::ShapeTraits<elemShape>::dim;
    const unsigned      dof        = 1;

    //! variables to be read from input file
    std::string meshFile, boundaryFile;
    const double conductivity = 1.;

    //! feed properties parser with the variables to read
    corlib::PropertiesParser *prop = new corlib::PropertiesParser;
    prop -> registerPropertiesVar( "meshFile",         meshFile   );
    prop -> registerPropertiesVar( "boundaryFile",     boundaryFile );

    //! read variables from the input.dat file
    std::ifstream inputFile( "./input.dat" );
    FTL_VERIFY( inputFile.is_open() );
    prop -> readValues( inputFile );
    delete prop;
    inputFile.close( );

    //! set some typedefs
    typedef corlib::Shapefun<  elemShape, numNodesPE>                Sfun;
    typedef corlib::NodeBasic<dim>                                   BasisNode;
    typedef del2::fem::NodePotential<BasisNode>                      Node;
    typedef corlib::ElementBasic<Node,Sfun>                          BasisElement;
    typedef del2::fem::ElementPotential<BasisElement>                Element;
    typedef corlib::Mesh<Element>                                    Mesh;
    typedef corlib::SystemSolveEigenSparse                           SysSolve;
    typedef corlib::Quadrature<elemShape, numGaussPoints>            Quad;

    //! Derive boundary mesh attributes
    typedef corlib::ShapeFunType<Sfun::degree,dim-1,Sfun::sfc>::type SurfSFun;
    typedef corlib::ElementBasic<Node,SurfSFun>                      BasisSurfElement;
    typedef del2::fem::ElementPotential<BasisSurfElement>            SurfElement;

    //! Attributes for the boundary force computation 
    typedef eigenX::VectorSd<dof>                                    VecDof;
    typedef eigenX::VectorSd<dim>                                    VecDim;
    typedef boost::function< VecDof( const VecDim ) >                VecDimToVecDof;
    typedef corlib::Quadrature<SurfSFun::myShape,numGaussSurf>       SurfaceQuad;
    typedef corlib::BodyForce<SurfElement,VecDimToVecDof>            SurfaceForce;

    //! read smf file and construct mesh
    std::ifstream smf( meshFile.c_str( ) );
    FTL_VERIFY( smf.is_open() );
    Mesh mesh( smf );
    smf.close( );

    //! read mesh boundary
    std::ifstream bdry( boundaryFile.c_str() );
    FTL_VERIFY( bdry.is_open() );
    corlib::MeshBoundary<SurfElement> neumannBoundary( bdry, mesh );
    bdry.close();

    //! cache material to elements
    mesh.iterateOverElements( std::bind( &Element::setConductivity,
                                         std::placeholders::_1, conductivity ) );

    //! number dofs
    unsigned numdofs = mesh.iterateOverNodes( corlib::dofNumberFun( &Node::numberDOFs, 0 ) );
    std::cout << numdofs << " degrees of freedom" << std :: endl;

    //! system matrix object 
    SysSolve * sysmat  = new SysSolve( numdofs );

    //! quadrature object
    Quad quadrature;

    //! compute element stiffness matrices   
    corlib::MatrixComputeAndAssembleFun<const Element,Quad,SysSolve>
        mafStiffness( &Element::stiffnessIntegrand, 
                      &Element::getDofIndices,
                      &Element::getDofIndices,
                      quadrature, sysmat );
    mesh.iterateOverElements( mafStiffness );

    //! surface force computation
    SurfaceQuad  surfaceQuadrature;
    SurfaceForce surfaceForce( std::bind( normalDerivative<dim>, std::placeholders::_1 ) );
    corlib::Integrator<SurfaceQuad,SurfaceForce> surfaceForceIntegrator( surfaceQuadrature, surfaceForce );
    mesh.iterateOverNodes(    std::bind( &Node::clearForce, std::placeholders::_1 ) );
    neumannBoundary.iterateOverElements( surfaceForceIntegrator );

    //! assemble the  forces to the system matrices rhs
    mesh.iterateOverNodes( corlib::vectorAssemblerFun( &Node::getForce, 
                                                       &Node::copyDofArray,
                                                       sysmat ) );
    
    //! collect and apply constraints
    sysmat -> finishAssembly();

    //! Fix a DOF, since this is a pure Neumann problem
    sysmat -> fixDOF( 0 );

    //! Solve it
    sysmat -> solveSystem( );

    //! add nodal solutions to increment 
    mesh.iterateOverNodes( corlib::distributorFun( sysmat, 
                                                   &Node::copyDofArray, 
                                                   &Node::setIncrement ) );

    //! clear the stiffness matrix and RHS
    sysmat -> clearMatrix( );
    sysmat -> clearRhs( );

    //! Update nodal solution
    mesh.iterateOverNodes( std::bind( &Node::updatePotential, std::placeholders::_1 ) );

    //! write data to VTU file
    std::string vtuFile = "flow.vtu";
    writeMeshData( vtuFile, &mesh );

    //! Clear memory
    delete sysmat;

    return 0;
}

//------------------------------------------------------------------------------
template< typename MESH>
void writeMeshData( const std::string & vtuFile, MESH * mesh )
{
    // set up the VTU writer and open a file
    std :: ofstream vtu( vtuFile.c_str( ) );
    assert( vtu.is_open( ) );
    corlib::VTUwriter<MESH> vtuwriter( mesh, vtu );

    //! write mesh data
    vtuwriter.writeMesh();
    
    //-- write point data
    vtuwriter.openPointData();
    vtuwriter.writeNodalQuantity( corlib::accessorFun( &MESH::Node::getPotential, "Potential" ) );
    vtuwriter.closePointData( );

    // finish writing
    vtuwriter.finalize();
    vtu.close( );
}
