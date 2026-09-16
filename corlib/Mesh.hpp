// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//
//                 Computational Structural Mechanics Lab
//                         University of Cambridge
//
// This software is copyrighted and all rights are retained by the University of
// Cambridge. All use, disclosure, and/or reproduction of any part not expressly
// authorised by the University of Cambridge is prohibited. (C) 2010
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

//! @file Mesh.hpp

#ifndef corlib_mesh_h
#define corlib_mesh_h

//------------------------------------------------------------------------------
//! system includes
#include <vector>
#include <iostream>
#include <algorithm>
#include <functional>

//! boost includes
#include <boost/utility.hpp>
#include <boost/function.hpp>

//------------------------------------------------------------------------------
namespace corlib{
    template<typename ELEMENT> class Mesh;
}

//------------------------------------------------------------------------------
/** \brief container for the FE nodes and elements
 *  \details Has containers with pointers to the nodes and the elements
 * \tparam ELEMENT  the type of element to be used
 */
template<typename ELEMENT>
class corlib::Mesh
    : public boost::noncopyable
{
public:
    typedef ELEMENT                                       Element;
    typedef typename Element::Node                        Node;


protected:
    typedef std::vector<Node*>                            NodeVec_;
    typedef std::vector<Element*>                         ElementVec_;

public:
    typedef typename NodeVec_::const_iterator             NodeConstIterator;
    typedef typename ElementVec_::const_iterator          ElementConstIterator;

    //! @name Constructors and destructor
    //@{
    //! Constructor with SMF file
    Mesh( std::istream & smf );

    //! Read SMF file function
    void read( std::istream & smf );

    //! Constructor with node and element iterators
    template< typename IT1, typename IT2 >
    Mesh( IT1 vb, IT1 ve, IT2 eb, IT2 ee );

    //! Constructor with iterators with element pointers
    //! This constructor is used for duplicating meshes
    template< typename ITEP >
    Mesh( ITEP itepb, ITEP itepe,
          boost::function<bool(Node*,Node*)> nodeComp = std::less<Node*>() );

    //! Empty constructor
    Mesh( ) : isCopy_( false ) // set flag that says the mesh is just
                               // a copy of other mesh and should not deallocate pointers
    { }

    //! Destructor which clears the containers
    virtual ~Mesh( );
    //@}

public:
    //! @name Generic application of functors to node and element containers
    //@{
    //! Iterate over nodes
    template<typename OP>  OP iterateOverNodes(    OP op );
    //! Iterate over elements
    template<typename OP>  OP iterateOverElements( OP op );
    //! Iterate over nodes with predicate
    template<typename OP, typename PRED>
    OP iterateOverNodesWithPredicate( OP op, PRED p );
    //! Iterate over elements with predicate
    template<typename OP, typename PRED>
    OP iterateOverElementsWithPredicate( OP op, PRED p );
    //@}

    //! @name Access to node and element containers
    //@{
    //! Give const access to begin of node container
    NodeConstIterator nodesBegin( ) const { return nodes_.begin( ); }
    //! Give const access to end of node container
    NodeConstIterator nodesEnd(   ) const { return nodes_.end( ); }
    //! Give const access to begin of element container
    ElementConstIterator elementsBegin( ) const { return elements_.begin( ); }
    //! Give const access to end of element container
    ElementConstIterator elementsEnd(   ) const { return elements_.end( ); }
    //! Give const access to specific node
    const Node* getNodePointer( const unsigned numNode ) const { return nodes_.at(numNode);}
    //! Give const access to specific element
    const Element* getElementPointer( const unsigned numElem ) const { return elements_.at(numElem); }
    //! Give access to a node pointer
    Node* getNodePointer( const unsigned numNode ) { return nodes_.at(numNode); }
    //! Grant non-const access to specific element
    Element* getElementPointer( const unsigned numElem ) { return elements_.at(numElem); }
    //@}

    //! @name Size queries
    //@{
    //! Return number of nodes
    unsigned numNodes( ) const { return nodes_.size( ); }
    //! Return number of elements
    unsigned numElements( ) const { return elements_.size( ); }
    //@}

    //! Write an smf file from the mesh
    std::ostream & writeSmf( std::ostream & smf );

protected:
    NodeVec_        nodes_;       //!< Vector with node pointers
    ElementVec_     elements_;    //!< Vector with element pointers
private:
    const bool      isCopy_;    //!< Set to true if mesh only mirrors another mesh
                                //!< having pointers to data allocated elsewhere
};
//------------------------------------------------------------------------------
#include "Mesh.ipp"


#endif
