// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//
//                 Computational Structural Mechanics Lab
//                         University of Cambridge
//
// This software is copyrighted and all rights are retained by the University of
// Cambridge. All use, disclosure, and/or reproduction of any part not expressly
// authorised by the University of Cambridge is prohibited. (C) 2010
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

//! @file helpers.hpp

#include <corlib/Accessor.hpp>
#include <corlib/ElementBasic.hpp>
#include <corlib/NodeBasic.hpp>

#include <boost/function.hpp>
#include <boost/bind.hpp>

//------------------------------------------------------------------------------
void readFromInput( std::string & meshFile, 
                    std::string & materialFile, 
                    std::string & constraintsFile, 
                    unsigned    & numberOfSubsteps )
{
    // feed properties parser with the variables to read
    corlib::PropertiesParser *prop = new corlib::PropertiesParser;
    prop -> registerPropertiesVar( "meshFile",         meshFile   );
    prop -> registerPropertiesVar( "materialFile",     materialFile );
    prop -> registerPropertiesVar( "constraintsFile",  constraintsFile );
    prop -> registerPropertiesVar( "numberOfSubsteps", numberOfSubsteps );

    // read variables from the input.dat file
    std::ifstream inputFile( "./input.dat" );
    FTL_VERIFY( inputFile.is_open() );
    prop -> readValues( inputFile );
    delete prop;
    inputFile.close( );
}

//------------------------------------------------------------------------------
// write solution data to a file
template< typename MESH>
void writeMeshData( const std::string & vtuFile, MESH * mesh )
{
    // set up the VTU writer and open a file
    std::ofstream vtu( vtuFile.c_str( ) );
    corlib::VTUwriter<MESH> vtuwriter( mesh, vtu );

    // write mesh data
    vtuwriter.writeMesh();
    
    //-- write point data
    vtuwriter.openPointData();
    vtuwriter.writeNodalQuantity( corlib::accessorFun( &MESH::Node::giveDisplacements, "Displacements" ) );
    vtuwriter.closePointData( );
    
    // -- write element data
    vtuwriter.openCellData( );
    typedef boost::function< typename MESH::Element::Mat3x3( const typename MESH::Element * ) > Mat3x3Fun;
    Mat3x3Fun giveStress = boost::bind( &MESH::Element::firstPiolaKirchhoff, _1,
                                        corlib::ShapeTraits< MESH::Element::myShape >::centroid( ) );
    corlib::Accessor< Mat3x3Fun > accessStress( giveStress, "Stress" );
    vtuwriter.writeElementQuantity( accessStress );
    vtuwriter.closeCellData( );

    // finish writing
    vtuwriter.finalize();
    vtu.close( );
}

//------------------------------------------------------------------------------
template<corlib::shape SHAPE, unsigned numNodesPE, unsigned numGaussP>
struct Attributes
{
    static const unsigned dim = corlib::ShapeTraits<SHAPE>::dim;

    typedef corlib::Quadrature<SHAPE, numGaussP>                          Quad;
    typedef corlib::Shapefun<  SHAPE, numNodesPE>                         Sfun;
    typedef corlib::NodeBasic<dim>                                        BasisNode;
    typedef solid::fem::NodeStatic<BasisNode>                             Node;
    typedef solid::material::MaterialBase                                 Material;
    typedef corlib::ElementBasic<Node,Sfun>                               BasisElement;
    typedef solid::fem::ElementStatic<BasisElement,Material>              Element;
    typedef corlib::Mesh<Element>                                         Mesh;
    typedef boost::function<void( Element*, 
                                  const typename Element::VecLDim &, 
                                  const double & ) >                      Integrand;
    typedef corlib::Integrator<Quad,Integrand,Element>                    Integrator;
};
