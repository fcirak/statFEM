// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//
//                 Computational Structural Mechanics Lab
//                         University of Cambridge
//
// This software is copyrighted and all rights are retained by the University of
// Cambridge. All use, disclosure, and/or reproduction of any part not expressly
// authorised by the University of Cambridge is prohibited. (C) 2010
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

//! @file laplaceRHSinfer.cpp

#include <iostream>
#include <fstream>
#include <string>
#include <map>
#include <cstdlib>
#include <functional>
#include <utility>

//！boost lib headers
#include <boost/function.hpp>

//！corlib headers
#include <corlib/Mesh.hpp>
#include <corlib/Shapefun.hpp>
#include <corlib/Quadrature.hpp>
#include <corlib/DofNumberer.hpp>
//#include <corlib/SystemSolve.hpp>
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

//！local headers
#include <del2/fem/NodePotential.hpp>
#include <del2/fem/ElementPotential.hpp>

#include "fem.hpp"
#include "inference.hpp"

namespace eigenX = corlib::eigenX;

//------------------------------------------------------------------------------
int main( int argc, char* argv[] )
{
    //! variables to be read from input file "./input.dat"
    std::string meshFile, allParamFile, infInputFile;

    //! instantiate and feed properties parser with the variables to read
    corlib::PropertiesParser *prop = new corlib::PropertiesParser;
    prop->registerPropertiesVar( "meshFile", meshFile );
    prop->registerPropertiesVar( "allParamFile", allParamFile );
    prop->registerPropertiesVar( "infInputFile", infInputFile );

    //! read variables from the input.dat file
    std::ifstream inputFile( "./input.dat" );
    FTL_VERIFY( inputFile.is_open( ) );
    prop -> readValues( inputFile );

    //! delete properties parser and close the input file
    delete prop;
    inputFile.close( );

    //! indicate if we want to compute with line, triangles or quadrilaterals
    //! this has to, of course, correlate with the element type in meshfile
    const corlib::shape elemShape = corlib::LINE;
    const unsigned dim  = corlib::ShapeTraits<elemShape>::dim;
    const corlib::SfunClass sfunClass = corlib::SIMPLEX;
    const unsigned numQuadPoints = 3;
    constexpr unsigned DEGREE = 1;
    typedef corlib::ShapeFunType<DEGREE,dim,sfunClass>::type Sfun;

    //! instantiate a quadrature object
    typedef corlib::Quadrature<elemShape,numQuadPoints>  Quad;
    Quad quadrature;

    //! instantiate the data structures for computing
    typedef corlib::NodeBasic<dim>                        BasisNode;
    typedef del2::fem::NodePotential<BasisNode>           Node;
    typedef corlib::ElementBasic<Node,Sfun>               BasisElement;
    typedef del2::fem::ElementPotential<BasisElement>     Element;
    typedef corlib::Mesh<Element>                         Mesh;
    typedef corlib::SystemSolveEigenSparse                SysSolve;

    //! read smf file and construct mesh
    std::ifstream smf( meshFile.c_str( ) );
    FTL_VERIFY( smf.is_open( ) );
    Mesh mesh( smf );
    smf.close( );

    //! number dofs
    unsigned numDofs = 0;
    mesh.iterateOverNodes( boost::bind( &Node::numberDOFs, _1, boost::ref( numDofs ) ) );

    //! get the parameters
    std::map<std::string, double> allParamMap;
    apps::readParamToMap( allParamFile, allParamMap );

    //! cache material to elements
    double conductivity = allParamMap[ "conductivity" ];
    mesh.iterateOverElements( boost::bind( &Element::setConductivity, _1, conductivity ) );

    //! system matrix object
    SysSolve * sysmat  = new SysSolve( numDofs );

    //! call the finite element module
    fem( mesh, quadrature, sysmat );
    sysmat -> solveSystem( );

    //! call the inference module
    inference( numDofs, mesh, quadrature, sysmat, infInputFile, allParamMap );

    //! clear the stiffness matrix and RHS
    sysmat -> clearMatrix( );
    sysmat -> clearRhs( );

    delete sysmat;

    return 0;
}
