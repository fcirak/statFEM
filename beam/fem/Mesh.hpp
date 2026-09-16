// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//
//                 Computational Structural Mechanics Lab
//                         University of Cambridge
//
// This software is copyrighted and all rights are retained by the University of
// Cambridge. All use, disclosure, and/or reproduction of any part not expressly
// authorised by the University of Cambridge is prohibited. (C) 2010
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

//! @file   Mesh.hpp

#ifndef beam_fem_meshbeam_h
#define beam_fem_meshbeam_h
//------------------------------------------------------------------------------
#include <vector>
#include <array>
#include <string>
#include <iostream>
#include <algorithm>

#include <corlib/Mesh.hpp>
#include <beam/fem/LinkBasic.hpp>

//------------------------------------------------------------------------------
// declarations
namespace beam {
    namespace fem {

        template<typename ELEMENT> class Mesh;
        
    }
}

//------------------------------------------------------------------------------
/** \brief
 *  Container for the FE nodes and beam elements, extends corlib::Mesh< NODE, ELEMENT >
 *
 *  \details
 *  The base class of the object reads the input stream and generates the 
 *  (beam) elements. The base class stores containers with pointers to the nodes 
 *  and the elements.  Special about this derived class is its identifying
 *  of the neighbours of the beam elements.  The latter can find their set of
 *  interpolating nodes looking at their neighbours and knowing their
 *  B-spline shape funtion degree.
 *
 *  \tparam NODE     the type of node to be used
 *  \tparam ELEMENT  the type of element to be used
 */
template<typename ELEMENT >
class beam::fem::Mesh :
    public corlib::Mesh< ELEMENT >
{
public:
    typedef ELEMENT                                            Element;
    typedef typename Element::Node                             Node;
    typedef LinkBasic< Node >                                  Link;
                                                          
protected:                                                
    typedef std::vector<Node*>                                 NodeVec_;
    typedef std::vector<Element*>                              ElementVec_;
    typedef std::vector<Link*>                                 LinkVec_;

    typedef typename std::array<Element*, 2>                   VecEle2_;
    typedef typename std::map< Element*, VecEle2_ >            MapEleToEle2_;
    typedef typename MapEleToEle2_::iterator                   MapEleToEle2Iter_;
    typedef typename MapEleToEle2_::const_iterator             MapEleToEle2ConstIter_;

public:                                                   
    typedef typename NodeVec_::const_iterator                  NodeConstIterator;
    typedef typename ElementVec_::const_iterator               ElementConstIterator;
    typedef typename LinkVec_::const_iterator                  LinkConstIterator;

    //! @name Constructors and destructor
    //@{
    //! Constructor with SMF file
    Mesh( std::istream& smf );

    //! Destructor which clears the containers
    virtual ~Mesh( );
    //@}

    //! mesh topology type
    bool hasClosedTopology() const { return closedTopology_; }

    //! Give node at provided node ID
    //!
    //! \param[in]   nodeId   node ID as found in input files
    Node* nodeFind( const unsigned nodeId ) const;

    //! Give node which is stored previous to given node
    //!
    //! \param[in]   theNode   pointer to reference node
    Node* nodeFindPrevious( const Node * theNode ) const;

    //! Give node which is stored previous to given node
    //!
    //! \param[in]   theNode   pointer to reference node
    Node* nodeFindNext( const Node * theNode ) const;

    //! Determine if node is ghost node
    bool nodeIsGhost( const Node * theNode ) const
    {
        NodeConstIterator n = std::find( ghostNodes_.begin(), ghostNodes_.end(), theNode );
        if ( n != ghostNodes_.end() )
            return true;
        else
            return false;
    }

    //! Number of ghost nodes
    unsigned numNodesGhost() const { return ghostNodes_.size( ); }

    //! Convenience function to get left element out of map holding
    //! the neighbours of the elements
    //!
    //! \param[in]  ele   element of which left neighbour is sought
    //! \return           left neighbour or NULL
    Element* leftNeighbourElement( Element* ele ) const;

    //! Convenience function to get right element out of map holding
    //! the neighbours of the elements
    //!
    //! \param[in]  ele   element of which right neighbour is sought
    //! \return           right neighbour or NULL
    Element* rightNeighbourElement( Element* ele ) const;

    //! @name Functions to deal with links
    //@{

    //! Set list of links due to supports
    //!
    //! \tparam      SREADER         Supports reader class
    //! \param[in]   suppReader      Support input reader
    template< typename SREADER >
    void setLinksOfSupports( const SREADER& suppReader );

    //! Iterate over elements
    //! Apply a given functor OP to all links making use of std::for_each
    //!
    //! \tparam         OP  Type of the element operation 
    //! \param[in,out]  op  Specific operation applied to all element
    //! \retval         OP  for_each returns the operator
    template< typename OP >
    OP iterateOverLinks( OP op ); 

    //! Return number of links
    unsigned numLinks() const { return links_.size(); }

    //@}

    const Link * getLinkPointer( const unsigned numLink ) const
    {
        return links_.at( numLink );
    }

    Link * getLinkPointer( const unsigned numLink )
    {
        return links_.at( numLink );
    }

    void clearLinks() {
        this -> iterateOverLinks( corlib::deleteFunctor( ) );
        links_.clear( );
    }

    //@name Output
    //@{

    //! Return SMF-like stream 
    std::ostream & writeSmf( std::ostream & os,
                             bool withSmfHead = true ) const;

    //@}

private:
    //! Create mesh which is topologically a ring (closed curve)
    void createMeshClosed_();

    //! Create mesh which is topologically an open curve
    void createMeshOpen_();

protected:
    //! List of nodes <i>including</i> ghost nodes
    using corlib::Mesh< Element >::nodes_;
    //! List of pointers to ghost nodes
    NodeVec_        ghostNodes_;

protected:
    //! a closed loop mesh
    bool            closedTopology_;

    //! Map holding the neighbours of the elements
    MapEleToEle2_   mapEleToFaceEle_;

    //! List of ghost elements only needed to find ghost nodes
    ElementVec_     ghostElements_;

    //! List of links
    LinkVec_        links_;

};
//------------------------------------------------------------------------------
#include "Mesh.ipp"


#endif
