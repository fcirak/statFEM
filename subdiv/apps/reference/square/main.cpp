//
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//
//                 Computational Structural Mechanics Lab
//                         University of Cambridge
//
// This software is copyrighted and all rights are retained by the CSMLab. All
// use, disclosure, and/or reproduction of any part not expressly authorized by
// F Cirak is prohibited. (C) 2010.
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//

//! @author Kosala Bandara
//! @date   2010

#include <iostream>
#include <fstream>
#include <string>

#include <corlib/Shape.hpp>
#include <corlib/verify.hpp>

#include <subdiv/surf/Vertex.hpp>
#include <subdiv/surf/Edge.hpp>
#include <subdiv/surf/Facet.hpp>
#include <subdiv/surf/FacetTree.hpp>

#include <subdiv/surf/Subdivision.hpp>
#include <subdiv/surf/Mesh.hpp>
#include <subdiv/surf/io.hpp>

#include "Parameters.hpp"

//==============================================================================
// everything starts here
int main() {

    //--------------------------------------------------------------------------
#if LOOP_
    const subdiv::surf::method subdivMethod = subdiv::surf::LOOP;
#elif CATMULL_CLARK_
    const subdiv::surf::method subdivMethod = subdiv::surf::CATMULL_CLARK;
#endif
    const corlib::shape myShape = subdiv::surf::SubdivisionShape< subdivMethod >::myShape;

    typedef subdiv::surf::Vertex                                       Vertex; 
    typedef subdiv::surf::Facet< myShape, Vertex >                     Facet;
    typedef subdiv::surf::Edge< Vertex, Facet >                        Edge;
    typedef subdiv::surf::FacetTree< myShape, Vertex >                 FacetTree;
    typedef subdiv::surf::Mesh< Vertex, Edge, Facet, FacetTree >       Mesh;

    typedef subdiv::surf::Subdivision< subdivMethod, FacetTree >       Subdivision;
    //--------------------------------------------------------------------------
    // read input data using Properties Parser
#if LOOP_
    const std::string inputData = "inputTria.dat";
#elif CATMULL_CLARK_
    const std::string inputData = "inputQuad.dat";
#endif

    std::ifstream inp( inputData.c_str() );
    FTL_VERIFY( inp.is_open() );
    std::cout << " Reading input data from: " << inputData << std::endl;
    app::Parameters parameters( inp );
    inp.close();

    //read input mesh and tags
    std::ifstream inputStream( parameters.inputMeshFileName.c_str() );
    FTL_VERIFY( inputStream.is_open() );
    std::ifstream secondaryInputStream( parameters.secondaryInputFileName.c_str() );
    FTL_VERIFY( secondaryInputStream.is_open() );

    Mesh mesh( inputStream );
    mesh.readTags( secondaryInputStream );
    mesh.buildFacetTopology( );  // OK
    mesh.implantFacetTrees( );

    // some debug output
    //mesh.iterateOverVertices( boost::bind( &subdiv::surf::write< 0 >, boost::ref( std::cout ), _1 ) );
    //mesh.iterateOverFacets( boost::bind( &subdiv::surf::write< myShape, Vertex >, boost::ref( std::cout ), _1 ) );

    // apply three levels of subdivsion
    Subdivision subdiv;
    mesh.subdivide< Subdivision >( &subdiv );
    mesh.subdivide< Subdivision >( &subdiv );
    mesh.subdivide< Subdivision >( &subdiv );

    //write output mesh 
    std::cout << " Writing subdivided mesh to: " << parameters.outputMeshFileName << std::endl;

    const std::string makeOutputFolder = "mkdir -p " + parameters.outputFolder;
    FTL_VERIFY( !system( makeOutputFolder.c_str( ) ) );
    std::ofstream outputStream( parameters.outputMeshFileName.c_str() );
    mesh.writeMeshWithIndex( outputStream );

    return 0;
}

