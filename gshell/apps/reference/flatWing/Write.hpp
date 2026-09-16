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
//! @date   2011

#ifndef gshell_app_write_h 
#define gshell_app_write_h

#include <iostream>
#include <string>
#include <algorithm>
#include <ctime>

#include <corlib/Accessor.hpp>

//==============================================================================
// declarations
namespace app{

    namespace eigenX = corlib::eigenX;

    template< typename MESH >
    void writeMeshData( const std::string & vtuFile, MESH * mesh );

    template< typename MESH >
    void writeHistory( std::ostream & histFile, const double time,
                       const unsigned timeStep, MESH & mesh );
    
    std::string secondsToPosixTime( const double & sec );

}

//==============================================================================
// definitions

//------------------------------------------------------------------------------
template< typename MESH >
void app::writeMeshData( const std::string & vtuFile, MESH * mesh )
{
    
    typedef typename MESH::NodeType    Node;
    typedef typename MESH::ElementType Element;

    // set up the VTU writer and open a file
    std::ofstream vtu( vtuFile.c_str() );
    FTL_VERIFY_DESCRIPTIVE( vtu.is_open(), "Failed to open file %s\n", vtuFile.c_str() );
    corlib::VTUwriter< MESH > vtuwriter( mesh, vtu );

    // write mesh data
    vtuwriter.writeMesh();
    
    //--------------------------------------------------------------------------
    // write point data
    vtuwriter.openPointData();
    
    const unsigned dim = Node::dim;
        
    // write displacements
    typedef boost::function< typename Node::VecDim( const Node * )> VecDimFun;
    VecDimFun giveDispl = std::bind( &Node::giveDisplacementsRange, std::placeholders::_1,
                                     Eigen::VectorXi::LinSpaced(dim, 0, dim) );
    corlib::Accessor< VecDimFun > accessDispl( giveDispl, "Displacements" );
    vtuwriter.writeNodalQuantity( accessDispl );
    // write velocities
#ifndef DYNAMICS_STATIC
    VecDimFun giveVeloc = std::bind( &Node::giveVelocitiesRange, std::placeholders::_1,
                                     Eigen::VectorXi::LinSpaced(dim, 0, dim) );
    corlib::Accessor< VecDimFun > accessVeloc( giveVeloc, "Velocities" );
    vtuwriter.writeNodalQuantity( accessVeloc );
#endif
    // write accelerations
#ifndef DYNAMICS_STATIC
    VecDimFun giveAccel = std::bind( &Node::giveAccelerationsRange, std::placeholders::_1,
                                     Eigen::VectorXi::LinSpaced(dim, 0, dim) );
    corlib::Accessor< VecDimFun > accessAccel( giveAccel, "Accelerations" );
    vtuwriter.writeNodalQuantity( accessAccel );
#endif
//#define WRITE_EQUIVALENTFORCES
#ifdef WRITE_EQUIVALENTFORCES
     // write element (or so-called equivalent) forces
     VecDimFun giveForce = std::bind( &Node::giveForceRange, std::placeholders::_1,
                                      Eigen::VectorXi::LinSpaced(dim, 0, dim) );
     corlib::Accessor< VecDimFun > accessForce( giveForce, "Force" );
     vtuwriter.writeNodalQuantity( accessForce );
#endif
//#define WRITE_TANGENTSANDNORMAL
#ifdef WRITE_TANGENTSANDNORMAL
     const gshell::fem::config conf = gshell::fem::CURRENT;
     // write 1st tangent
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
#endif
     // finish point data
    vtuwriter.closePointData();

    //--------------------------------------------------------------------------
    // write element data
    vtuwriter.openCellData();
//     // element area
//     typedef boost::function< double( const Element * )> DoubleFun;
//     DoubleFun giveArea = boost::bind( &Element::giveArea, _1 );
//     corlib::Accessor< DoubleFun > accessArea( giveArea, "Area" );
//     vtuwriter.writeElementQuantity( accessArea );
//     // membrane strains
    typedef boost::function< typename Element::Mat3x3( const Element * )> Mat3x3Fun;
//     Mat3x3Fun giveMembStrain = boost::bind( &Element::giveStrainResultant, _1,
//                                             Element::INPLANE, Element::GREENLAGRANGE,
//                                             corlib::ShapeTraits< Element::myShape >::centroid() );
//     corlib::Accessor< Mat3x3Fun > accessMembStrain( giveMembStrain, "MembraneStrain" );
//     vtuwriter.writeElementQuantity( accessMembStrain );
//     // bending strains
//     Mat3x3Fun giveBendStrain = boost::bind( &Element::giveStrainResultant, _1,
//                                             Element::OUTOFPLANE, Element::GREENLAGRANGE,
//                                             corlib::ShapeTraits< Element::myShape >::centroid() );
//     corlib::Accessor< Mat3x3Fun > accessBendStrain( giveBendStrain, "BendingStrain" );
//     vtuwriter.writeElementQuantity( accessBendStrain );
    // membrane stresss
    Mat3x3Fun giveMembStress = std::bind( &Element::giveStressResultant, std::placeholders::_1,
                                           Element::INPLANE, Element::PIOLAKIRCHHOFF2,
                                           corlib::ShapeTraits< Element::myShape >::centroid() );
    corlib::Accessor< Mat3x3Fun > accessMembStress( giveMembStress, "MembraneStress" );
    vtuwriter.writeElementQuantity( accessMembStress );
    // bending stresss
    Mat3x3Fun giveBendStress = std::bind( &Element::giveStressResultant, std::placeholders::_1,
                                           Element::OUTOFPLANE, Element::PIOLAKIRCHHOFF2,
                                           corlib::ShapeTraits< Element::myShape >::centroid() );
    corlib::Accessor< Mat3x3Fun > accessBendStress( giveBendStress, "BendingStress" );
    vtuwriter.writeElementQuantity( accessBendStress );
    // finish element data
    vtuwriter.closeCellData();

    // finish writing
    vtuwriter.finalize();
    vtu.close();
}

//------------------------------------------------------------------------------
//! Write average displacement, velocity and acceleration of body to file
template< typename MESH >
void app::writeHistory( std::ostream & histFile,
                        const double time,
                        const unsigned timeStep,
                        MESH & mesh )
{
    // convenience
    namespace eigenX = corlib::eigenX;
    typedef typename MESH::NodeType        Node;
    typedef typename Node::VecDim          VecDim;

    // header
    if ( timeStep == 0 ) {
        histFile << "# time timeStep dispX dispY dispZ velX velY velZ" << std::endl;
    }

    // number of nodes
    const double numNodes = static_cast< double >( mesh.numNodes() );

    // compute average displacements
    VecDim dis; dis.setZero();
    dis = std::accumulate( mesh.nodesBegin(), mesh.nodesEnd(), dis,
                           std::bind( std::plus< VecDim >(), std::placeholders::_1,
                                        std::bind( &Node::giveDisplacements, std::placeholders::_2 ) ) );
    dis /= numNodes;

    // compute average velocities
    VecDim vel = Eigen::VectorXd::Zero( Node::dim );
    vel = std::accumulate( mesh.nodesBegin(), mesh.nodesEnd(), vel,
                           std::bind( std::plus< VecDim >(), std::placeholders::_1,
                                        std::bind( &Node::giveVelocities, std::placeholders::_2 ) ) );
    vel /= numNodes;

    // write stuff
    histFile << time << " " << timeStep << " ";
    histFile << dis.transpose() << " ";
    histFile << vel.transpose() << " ";
    histFile << std::endl;
    
    return;
}

//-----------------------------------------------------------------------------
// helper function to convert seconds to HH:MM:SS
std::string app::secondsToPosixTime( const double & sec )
{
    std::time_t pt = 1000.*sec;
    return std::to_string( pt );
}

#endif
