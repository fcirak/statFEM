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

#ifndef subdiv_surf_wrapper_h
#define subdiv_surf_wrapper_h

#include <subdiv/surf/Vertex.hpp>
#include <subdiv/surf/Edge.hpp>
#include <subdiv/surf/Facet.hpp>
#include <subdiv/surf/FacetTree.hpp>
#include <subdiv/surf/Subdivision.hpp>
#include <subdiv/surf/Mesh.hpp>
#include <subdiv/surf/ShapeFunSubdivision.hpp>

//==============================================================================
namespace subdiv {
    namespace surf {

        template< subdiv::surf::method SUBDMETHOD >
        class Wrapper;

    }
}

//==============================================================================
/// Helper to encapsulate the subdivision surface
///
/// \tparam SUBDMETHOD   Subdivision method
template< subdiv::surf::method SUBDMETHOD >
class subdiv::surf::Wrapper
{
public:
    static const subdiv::surf::method subdivMethod   = SUBDMETHOD;
    static const corlib::shape        myShape        = subdiv::surf::SubdivisionShape<subdivMethod>::myShape;
    
    typedef subdiv::surf::Vertex                                       Vertex;
    typedef subdiv::surf::Facet< myShape, Vertex >                     Facet;
    typedef subdiv::surf::Edge< Vertex, Facet >                        Edge;
    typedef subdiv::surf::FacetTree< myShape, Vertex >                 FacetTree;
    typedef subdiv::surf::Mesh< Vertex, Edge, Facet, FacetTree >       Mesh;
    typedef subdiv::surf::Subdivision< subdivMethod, FacetTree >       Subdivision;
    typedef subdiv::surf::ShapeFunSubdivision< Facet, Subdivision >    ShapeFun;
    
public:
    /// Constructor
    Wrapper( std::istream & smf, std::istream * tg ) 
    {
        this -> read( smf, tg );
    }

    /// Empty c'tor to read mesh later
    Wrapper() { }

    /// Read and build the mesh
    void read( std::istream & smf, std::istream * tg )
    {
        //Vertex::setNumAddCoords( numAdVariables );
        mesh_.readSmf( smf );
        if ( tg ) mesh_.readTags( *tg );
        mesh_.buildFacetTopology( );
    }

    /// Access mesh
    Mesh & accessMesh( ) { return mesh_; }

    /// Subdivision scheme
    Subdivision & accessSubdivision( ) { return subdiv_; }

protected:
    /// Surface mesh
    Mesh            mesh_;
    /// Subdivision scheme
    Subdivision     subdiv_;
};

#endif
