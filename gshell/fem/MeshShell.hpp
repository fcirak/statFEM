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

#ifndef gshell_fem_meshshell_h
#define gshell_fem_meshshell_h

//------------------------------------------------------------------------------
#include <vector>
#include <string>
#include <iostream>
#include <algorithm>

#include <subdiv/surf/ShapeFunSubdivision.hpp>

#include <corlib/Mesh.hpp>

#include <gshell/fem/LinkRotational.hpp>

//------------------------------------------------------------------------------
// declarations
namespace gshell {
    namespace fem {

        template< typename GELEMENT, typename SHAPEFUN > class MeshShell;
        
    }
}

//------------------------------------------------------------------------------
/** \brief
 *  Container for the FE nodes and gshell elements, extends corlib::Mesh< GELEMENT >
 *
 *  \details
 *  This mesh object constructs a g-shell mesh using a readily subdivided mesh made of
 *  shells::fem::ElementShell and the respective nodes. The support nodes of 
 *  elements of the subdivided mesh is copied the g-shell mesh at hand.
 *
 *  \tparam GNODE     the type of node to be used
 *  \tparam GELEMENT  the type of element to be used
 */
template< typename GELEMENT, typename SHAPEFUN >
class gshell::fem::MeshShell :
    public corlib::Mesh< GELEMENT >
{
public:
    typedef GELEMENT                                           ElementType;
    typedef SHAPEFUN                                           ShapeFun;
    typedef typename ElementType::Node                         NodeType;

    typedef gshell::fem::LinkRotationalExt< ElementType >      LinkExt;
    typedef gshell::fem::LinkRotationalInt< ElementType >      LinkInt;
                                                          
private:                                                  
    typedef corlib::Mesh< ElementType >                        MeshBasic_;

    typedef std::vector< LinkExt * >                           VecLinkExt_;
    typedef std::vector< LinkInt * >                           VecLinkInt_;

public:
    typedef typename MeshBasic_::NodeConstIterator             NodeConstIterator;
    typedef typename MeshBasic_::ElementConstIterator          ElementConstIterator;

    typedef typename VecLinkExt_::iterator                     LinkExtIter;
    typedef typename VecLinkInt_::iterator                     LinkIntIter;

public:
    //! @name Constructors and destructor
    //!{

    //! Constructor with donator surface mesh and separate subdivision scheme
    template< typename SMESH, typename QUADRATURE >
    MeshShell( SMESH & sMesh, const QUADRATURE & quadrature )
        : MeshBasic_() 
    {
        this->initialise( sMesh );
        this->populateShapeFunctionCache( sMesh, quadrature );
    }


    //! Empty constructor to initiate mesh later
    MeshShell() : MeshBasic_() { }

    //! Initiate mesh using surface mesh and subdivision rule
    template< typename SMESH >
    void initialise( SMESH & sMesh );

    //! Populate the shape function cache by evaluating at quadrature points
    template< typename SMESH, typename QUADRATURE >
    void populateShapeFunctionCache( SMESH & sMesh,
                                     const QUADRATURE & quadrature );

    //! Destructor which clears the containers
    virtual ~MeshShell( );

    //!}

    //! @name Methods to handle links (ie clamping and creases)
    //!{

    //! Read links contained in stream
    //!
    //! Link file format:
    //!\beginverbatim
    //!     <numLinksExternal>  <numLinksInternal>
    //!     <nodeID> <elementID>
    //!     <nodeID> <elementID>
    //!     ...        ...
    //!     <nodeID> <elementID_0> <elementID_1> 
    //!     <nodeID> <elementID_0> <elementID_1>
    //!     ...        ...         ...
    //!\endverbatim
    //!
    //! \tparam          MESH   Mesh type
    //!
    //! \param[in]       mesh   Input mesh
    //! \param[in,out]   is     Input stream
    void readLinks( std::istream & is );

    //! Remove all links
    void clearLinks();

    //! Insert external Link
    void addLinkExt( const LinkExt & linkExt );

    //! Insert internal Link
    void addLinkInt( const LinkInt & linkInt );

    //! Total number of links
    unsigned numLinks() const
    {
        return ( linksExt_.size( ) + linksInt_.size( ) );
    }

    //! Access to begin iterator of external links
    LinkExtIter lExtBegin( ) { return linksExt_.begin( ); }
    //! Access to end iterator of external links
    LinkExtIter lExtEnd( ) { return linksExt_.end( ); }

    //! Access to begin iterator of internal links
    LinkIntIter lIntBegin( ) { return linksInt_.begin( ); }
    //! Access to end iterator of internal links
    LinkIntIter lIntEnd( ) { return linksInt_.end( ); }

    //! Iterate over external elements
    template< typename OP >
    OP iterateOverLinksExt( OP op );

    //! Iterate over internal elements
    template< typename OP >
    OP iterateOverLinksInt( OP op );

    //!}

protected:
    //! List of nodes
    using MeshBasic_::nodes_;
    //! List of elements
    using MeshBasic_::elements_;

    /// Vector of external links (ie clamping)
    VecLinkExt_     linksExt_;
    /// Vector of internal links (ie creases, t-intersections, ...)
    VecLinkInt_     linksInt_;

};
//------------------------------------------------------------------------------
#include "MeshShell.ipp"


#endif
