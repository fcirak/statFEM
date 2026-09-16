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

//! @author Burkhard Bornemann, Kosala Bandara
//! @date   2010

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <functional>

// includes from corlib
#include <corlib/PropertiesParser.hpp>
#include <corlib/verify.hpp>
#include <corlib/VTUwriter.hpp>

// includes from subdiv::surf
#include <subdiv/surf/Vertex.hpp>
#include <subdiv/surf/Edge.hpp>
#include <subdiv/surf/Facet.hpp>
#include <subdiv/surf/FacetTree.hpp>
#include <subdiv/surf/Subdivision.hpp>
#include <subdiv/surf/Mesh.hpp>

// local includes
#include "Parameters.hpp"
#include "TaggingFun.hpp"

//==============================================================================
// everything starts here
int main(int argc, char* args[] ){

    //--------------------------------------------------------------------------
#if LOOP_
    const subdiv::surf::method subdivMethod = subdiv::surf::LOOP;
#elif CATMULL_CLARK_
    const subdiv::surf::method subdivMethod = subdiv::surf::CATMULL_CLARK;
#endif
    // simplex shape
    const corlib::shape myShape = subdiv::surf::SubdivisionShape< subdivMethod >::myShape;

    typedef subdiv::surf::Vertex                                       Vertex; 
    typedef subdiv::surf::Facet< myShape, Vertex >                     Facet;
    typedef subdiv::surf::Edge< Vertex, Facet >                        Edge;
    typedef subdiv::surf::FacetTree< myShape, Vertex >                 FacetTree;
    typedef subdiv::surf::Mesh< Vertex, Edge, Facet, FacetTree >       Mesh;

    typedef subdiv::surf::Subdivision< subdivMethod, FacetTree >       Subdivision;

    // read input data using Properties Parser
    const std::string inputData = "input.dat";
    std::ifstream inp( inputData.c_str() );
    app::Parameters parameters( inp );
    inp.close();

    //read input mesh and tags
    std::ifstream inputStream( parameters.inputMeshFileName.c_str() );
    FTL_VERIFY_DESCRIPTIVE( inputStream.is_open(), 
                            "Could not open file %s\n",
                            parameters.inputMeshFileName.c_str() );
    std::ifstream secondaryInputStream( parameters.secondaryInputFileName.c_str() );
    FTL_VERIFY_DESCRIPTIVE( secondaryInputStream.is_open(),
                            "Could not open file %s\n",
                            parameters.secondaryInputFileName.c_str() );

    // Surface
    Mesh sMesh( inputStream );
    sMesh.readTags( secondaryInputStream );  // prior to buildFacetTopo
    sMesh.buildFacetTopology( );  // OK
    sMesh.iterateOverVertices( &app::taggingVertex< Vertex > );
    sMesh.iterateOverFacets( &app::taggingEdge< Facet > );
    // make the tree roots
    sMesh.implantFacetTrees( );

    // Subdivision scheme
    Subdivision subdiv;

    // Make animation folder
    const std::string makeOutputFolder = "mkdir -p " + parameters.outputFolder;
    FTL_VERIFY( !system( makeOutputFolder.c_str( ) ) );

    // Make complete subdivision
    const unsigned numSubdivisions = 1;
    for ( unsigned s = 0; s <= numSubdivisions; ++s ) {
        // subdivide
        std::cout << s << ". subdivision" << std::endl;
        if ( s > 0 ) {
            sMesh.subdivide< Subdivision >( &subdiv );
        }

        // write resulting file
        const std::string outputSmfFileName = parameters.outputBaseFileName +
            "." + std::to_string( s ) + ".smf";
        const std::string outputTgFileName = parameters.outputBaseFileName +
            "." + std::to_string( s ) + ".tg";
        std::ofstream outputSmfStream( outputSmfFileName.c_str() );
        std::ofstream outputTgStream( outputTgFileName.c_str() );
        sMesh.writeMeshWithIndex( outputSmfStream, &outputTgStream );
        outputSmfStream.close( );
        outputTgStream.close( );

        // determine number of non-manifold edges
        unsigned nonManifEdges = 0;
        for ( Mesh::FacetIterator f = sMesh.fBegin( ); f != sMesh.fEnd( ); ++f ) {
            FacetTree * ft = sMesh.facetTree( *f );
            typedef std::vector< FacetTree * >                    VecFTPtr;
            typedef VecFTPtr::iterator                            VecFTPtrIter;
            VecFTPtr leafs;  ft->addLeafs( std::back_inserter( leafs ) );
            for ( VecFTPtrIter l = leafs.begin(); l != leafs.end(); ++l ) {
                for ( unsigned n = 0; n < Facet::numNeighbors; ++n )
                    if ( not (*l)->isManifold( n ) ) nonManifEdges++;
            }
        }
        std::cout << "nonManifEdges=" << nonManifEdges << std::endl;
    }
    
    // done
    return 0;
}
    
