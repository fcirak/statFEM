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

#ifndef gshell_apps_write_h
#define gshell_apps_write_h

//==============================================================================
namespace app {

    template< typename NODE >
    typename NODE::VecDim giveCartesianShear( const NODE * node );

    template< typename MESH >
    void writeMeshData( const std::string & vtuFile, MESH * mesh );

    template< typename MESH, typename DATAWRITER >
    void writeMeshCellData( const std::string & vtuFile, MESH * mesh,
                            DATAWRITER * writer);

}

//==============================================================================
template< typename NODE >
typename NODE::VecDim app::giveCartesianShear( const NODE * node )
{
    typedef typename NODE::VecDim VecDim;
        
    const VecDim shear = 
        node->giveDisplacementsRange( Eigen::VectorXi::LinSpaced(NODE::dim, 0, NODE::dim) );
    VecDim cartShear;  cartShear.setZero();
    for ( unsigned d = 0; d < (NODE::dof-NODE::dim); ++d ) {
        cartShear += shear( d ) * node->getTangent( gshell::fem::CURRENT, d );
    }
    return cartShear;
}

//==============================================================================
template< typename MESH >
void app::writeMeshData( const std::string & vtuFile, MESH * mesh )
{
    typedef typename MESH::NodeType    Node;
    typedef typename MESH::ElementType Element;

    // set up the VTU writer and open a file
    std::ofstream vtu( vtuFile.c_str( ) );
    FTL_VERIFY_DESCRIPTIVE( vtu.is_open( ),
                            "Failed to open file %s\n", vtuFile.c_str( ) );
    corlib::VTUwriter< MESH > vtuwriter( mesh, vtu );

    // write mesh data
    vtuwriter.writeMesh();
    
    //--------------------------------------------------------------------------
    // write point data
    vtuwriter.openPointData();
    // write displacements
    
    const unsigned dim = Node::dim;
        
    typedef std::function< typename Node::VecDim( const Node * )> VecDimFun;
    VecDimFun giveDispl = std::bind( &Node::giveDisplacementsRange, std::placeholders::_1,
                                     Eigen::VectorXi::LinSpaced(dim, 0, dim) );
    corlib::Accessor< VecDimFun > accessDispl( giveDispl, "Displacements" );
    vtuwriter.writeNodalQuantity( accessDispl );
    // write shears
//#ifdef KINEMATICS_SHEARFLEXIBLE
//    //VecDimFun giveShear = std::bind( &Node::giveDisplacementsRange, std::placeholders::_1,
//    //                                  Eigen::VectorXi::LinSpaced(dim, 0, dim) );
//    VecDimFun giveShear = std::bind( &app::giveCartesianShear<Node>, std::placeholders::_1 );
//    corlib::Accessor< VecDimFun > accessShear( giveShear, "Shears" );
//    vtuwriter.writeNodalQuantity( accessShear );
//#endif
    // write equivalent forces
    VecDimFun giveForce = std::bind( &Node::giveForceRange, std::placeholders::_1,
                                     Eigen::VectorXi::LinSpaced(dim, 0, dim) );
    corlib::Accessor< VecDimFun > accessForce( giveForce, "Force" );
    vtuwriter.writeNodalQuantity( accessForce );
    // write 1st tangent
    const gshell::fem::config conf = gshell::fem::CURRENT;
    VecDimFun giveTangent0 = std::bind( &Node::getTangent, std::placeholders::_1, conf, 0 );
    corlib::Accessor< VecDimFun > accessTangent0( giveTangent0, "Tangent0" );
    vtuwriter.writeNodalQuantity( accessTangent0 );
    // write 2nd tangent
    VecDimFun giveTangent1 = std::bind( &Node::getTangent, std::placeholders::_1, conf, 1 );
    corlib::Accessor< VecDimFun > accessTangent1( giveTangent1, "Tangent1" );
    vtuwriter.writeNodalQuantity( accessTangent1 );
    // write normal
    VecDimFun giveNormal = std::bind( &Node::getNormal, std::placeholders::_1, conf );
    corlib::Accessor< VecDimFun > accessNormal( giveNormal, "Normal" );
    vtuwriter.writeNodalQuantity( accessNormal );
    // finish of point data
    vtuwriter.closePointData( );

    //--------------------------------------------------------------------------
    // write element data
    vtuwriter.openCellData( );
    // element area
    typedef std::function< double( const Element * )> DoubleFun;
    DoubleFun giveArea = std::bind( &Element::giveArea, std::placeholders::_1 );
    corlib::Accessor< DoubleFun > accessArea( giveArea, "Area" );
    vtuwriter.writeElementQuantity( accessArea );
    // membrane strains
    typedef std::function< typename Element::Mat3x3( const Element * )> Mat3x3Fun;
    Mat3x3Fun giveMembStrain = std::bind( &Element::giveStrainResultant, std::placeholders::_1,
                                            Element::INPLANE, Element::GREENLAGRANGE,
                                            corlib::ShapeTraits< Element::myShape >::centroid( ) );
    corlib::Accessor< Mat3x3Fun > accessMembStrain( giveMembStrain, "MembraneStrain" );
    vtuwriter.writeElementQuantity( accessMembStrain );
    // bending strains
    Mat3x3Fun giveBendStrain = std::bind( &Element::giveStrainResultant, std::placeholders::_1,
                                            Element::OUTOFPLANE, Element::GREENLAGRANGE,
                                            corlib::ShapeTraits< Element::myShape >::centroid( ) );
    corlib::Accessor< Mat3x3Fun > accessBendStrain( giveBendStrain, "BendingStrain" );
    vtuwriter.writeElementQuantity( accessBendStrain );
    // membrane stresss
    Mat3x3Fun giveMembStress = std::bind( &Element::giveStressResultant, std::placeholders::_1,
                                            Element::INPLANE, Element::PIOLAKIRCHHOFF2,
                                            corlib::ShapeTraits< Element::myShape >::centroid( ) );
    corlib::Accessor< Mat3x3Fun > accessMembStress( giveMembStress, "MembraneStress" );
    vtuwriter.writeElementQuantity( accessMembStress );
    // bending stresss
    Mat3x3Fun giveBendStress = std::bind( &Element::giveStressResultant, std::placeholders::_1,
                                            Element::OUTOFPLANE, Element::PIOLAKIRCHHOFF2,
                                            corlib::ShapeTraits< Element::myShape >::centroid( ) );
    corlib::Accessor< Mat3x3Fun > accessBendStress( giveBendStress, "BendingStress" );
    vtuwriter.writeElementQuantity( accessBendStress );
    // finish of element data 
    vtuwriter.closeCellData( );

    // finish writing
    vtuwriter.finalize();
    vtu.close( );
}

//==============================================================================
template< typename MESH, typename DATAWRITER >
void app::writeMeshCellData( const std::string & vtuFile, MESH * mesh, 
                             DATAWRITER * writer )
{
    // set up the VTU writer and open a file
    std::ofstream vtu( vtuFile.c_str( ) );
    FTL_VERIFY_DESCRIPTIVE( vtu.is_open( ),
                            "Failed to create file %s\n", vtuFile.c_str( ) );
    corlib::VTUwriter<MESH> vtuwriter( mesh, vtu );


    // write mesh data
    vtuwriter.writeMesh();
  
    // -- write element data
    vtuwriter.openCellData( );
    vtuwriter.writeElementQuantity( *writer );
    vtuwriter.closeCellData( );

    // finish writing
    vtuwriter.finalize();
    vtu.close( );
}

#endif
