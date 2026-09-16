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

#include <corlib/Shape.hpp>
#include <corlib/SmfHead.hpp>
#include <corlib/ElementBasic.hpp>
#include <corlib/NodeBasic.hpp>
#include <corlib/misc.hpp>
#include <functional>
#include <beam/fem/Supports.hpp>

//------------------------------------------------------------------------------
/** Constructor given a stream to a simple mesh file. This file must have the
 *  following format:
 *
 *    NN   NE \n
 *    x1  y1  z1  \n
 *    ....        \n
 *    xNN yNN zNN \n
 *    v1_1  v1_2 ... v1_NVE   \n
 *    ...                     \n
 *    vNE_1 vNE_2 ... vNE_NVE \n
 *  
 * where NN is the number of nodes, NE the number of elements, and NVE refers to
 * the number of vertices per element. Note that the latter number is fixed at 
 * compile time. Note that the latter number NVE differs from the number of 
 * nodes per element, needed by the interpolation appraoch.
 *
 * \param[in]  smf   Input stream
 */
template< typename ELEMENT >
beam::fem::Mesh< ELEMENT >::Mesh( std::istream& smf ) :
    corlib::Mesh< ELEMENT >( smf )
{
    // element range
    const ElementConstIterator eb = this->elementsBegin();
    const ElementConstIterator ee = this->elementsEnd();

    // loop elements and collect nodes
    typedef typename std::multimap< Node*, Element* > NodeToElementMap;
    NodeToElementMap nodesToElems;
    for ( ElementConstIterator e=eb; e!=ee; ++e ) {
        typename ELEMENT::NodeIterator nIter = (*e)->nodesBegin();
        for ( ; nIter != (*e)->nodesEnd(); ++nIter ) {
            Node* np = *nIter;
            nodesToElems.insert( std::pair< Node*, Element* >( np, *e ) );
        }
    }

    // loop elements and store neighbours
    typedef typename NodeToElementMap::iterator NodeToElementMapIt;
    closedTopology_ = true;
    for ( ElementConstIterator e=eb; e!=ee; ++e ) {
        VecEle2_ faceEle = {{ NULL, NULL }};
        typename ELEMENT::NodeIterator nIter = (*e)->nodesBegin();
        for ( ; nIter != (*e)->nodesEnd(); ++nIter ) {
            Node* np = *nIter;
            const NodeToElementMapIt elemOfNodeLow = nodesToElems.lower_bound( np );
            const NodeToElementMapIt elemOfNodeUpp = nodesToElems.upper_bound( np );
            bool haveNeighbour = false;
            for ( NodeToElementMapIt eon=elemOfNodeLow; eon!=elemOfNodeUpp; ++eon ) {
                if ( (*eon).second != *e ) {
                    const unsigned n = std::distance( (*e)->nodesBegin(), nIter );                    
                    faceEle[ n ] = (*eon).second;
                    haveNeighbour = true;
                    break;
                }
            }
            if ( not haveNeighbour ) {
                closedTopology_ = false;
            }
        }
        //std::cout << *e << " : " << faceEle << std::endl;
        mapEleToFaceEle_.insert( std::make_pair( *e, faceEle ) );
    }

    // create the mesh
    if ( closedTopology_ ) {
        std::cout << "Detected closed mesh" << std::endl;
        this->createMeshClosed_();
    }
    else {
        std::cout << "Detected open mesh" << std::endl;
        this->createMeshOpen_();
    }

    return;
}

//------------------------------------------------------------------------------
//! Create mesh which is topologically a ring (closed curve)
template< typename ELEMENT >
void beam::fem::Mesh< ELEMENT >::createMeshClosed_()
{
    // element range
    const ElementConstIterator eb = this->elementsBegin();
    const ElementConstIterator ee = this->elementsEnd();

    // loop elements to set their nodes
    for ( ElementConstIterator e=eb; e!=ee; ++e ) {
        // number of interpolating nodes of element
        const unsigned numNodesSN = Element::numNodesSN;

        // find left most element relevant for interpolation
        const unsigned degree = numNodesSN - 1;
        //const unsigned offSet = degree % 2;

        // local number of left vertex
        // local ID   0     1     2     3           d (degree)
        //            +.....+.....+-----+.....+.....+
        // global ID              k    k+1
        const unsigned leftVertexLocId = degree / 2;  // k
        const unsigned rightVertexLocId = leftVertexLocId + 1;  // k+1

        // initialise node pointers
        for ( unsigned n=0; n<numNodesSN; ++n )
            (*e)->setNode( n, NULL );

        // find interpolating nodes on left side
        Element* leftMost = *e;
        for ( int n=leftVertexLocId; n>=0; --n ) {
            (*e)->setNode( n, *(leftMost->nodesBegin()) );
            if ( this->leftNeighbourElement( leftMost ) != NULL )
                leftMost = this->leftNeighbourElement( leftMost );
            else
                break;
        }

        // find interpolating nodes on right side
        Element* rightMost = *e;
        for ( unsigned n=rightVertexLocId; n<numNodesSN; ++n ) {
            (*e)->setNode( n, *(rightMost->nodesBegin()+1) );
            if ( this->rightNeighbourElement( rightMost ) != NULL )
                rightMost = this->rightNeighbourElement( rightMost );
            else
                break;
        }
    }

    return;
}

//------------------------------------------------------------------------------
//! Create mesh which is topologically an open curve
template< typename ELEMENT >
void beam::fem::Mesh< ELEMENT >::createMeshOpen_()
{
    // element range
    const ElementConstIterator eb = this->elementsBegin();
    const ElementConstIterator ee = this->elementsEnd();

    // find element without left neighbour
    for ( ElementConstIterator e=eb; e!=ee; ++e ) {
        if ( this->leftNeighbourElement( *e ) == NULL ) {
            // Create co-ordinate of ghost node
            // Ghost nodes depending on spline order/degree.
            //
            // Example 7th order spline
            //          <===== ghost nodes 'o' created in this direction
            //    o.....o.....o.....+=====+-----+-----+-----+
            //                      k    k+1
            const unsigned numGhostVertices = ELEMENT::numNodesSN / 2 - 1;
            const typename Node::VecDim leftCoor = (*((*e)->nodesBegin()))->giveCoordinates();
            Element* rightElement = (*e);
            Element* mirrorRightElement = (*e);
            for ( int g=numGhostVertices; g>0; --g ) {
                // Creation direction
                // 
                // Example 7th degree spline:
                //               <===== created in this direction
                //       o.....o.....o.....+=====+-----+-----+-----+
                //             :     : 
                //             o.....o.....+-----+=====+-----+-----+-----+
                //                   :
                //                   o.....+-----+-----+=====+-----+-----+-----+

                // Create ghost element
                Element* ghostElement = new Element( );
                ghostElements_.push_back( ghostElement );
                VecEle2_ faceEle;
                faceEle[ 0 ] = NULL;
                faceEle[ 1 ] = rightElement;
                mapEleToFaceEle_.insert( std::make_pair( ghostElement, faceEle ) );

                // Attach ghost element
                faceEle[ 0 ] = ghostElement;
                faceEle[ 1 ] = this->rightNeighbourElement( rightElement );
                mapEleToFaceEle_.find( rightElement )->second = faceEle;

                // Create ghost node
                Node* ghostNode = new Node;
                ghostNode -> setId( nodes_.size() );
                nodes_.push_back( ghostNode );  // store node
                ghostNodes_.push_back( ghostNode );  // store node
                const typename Node::VecDim rightCoor = 
                    (*(mirrorRightElement->nodesBegin()+1))->giveCoordinates();
                ghostNode->setCoordinates( 2.0*leftCoor - rightCoor );
//                ghostNode->setCoordinates( 2.0*leftCoor - 3.0*rightCoor );
                mirrorRightElement = this->rightNeighbourElement( mirrorRightElement );

                // Attach ghost node to ghost element
                *(ghostElement->nodesBegin()) = ghostNode;
                *(ghostElement->nodesBegin()+1) = *(rightElement->nodesBegin());

                // advance to left
                rightElement = this->leftNeighbourElement( rightElement );
            }
        }
    }

    // find element without right neighbour
    for ( ElementConstIterator e=eb; e!=ee; ++e ) {
        if ( this->rightNeighbourElement( *e ) == NULL ) {
            // create co-ordinate of ghost node
            const unsigned numGhostVertices = ELEMENT::numNodesSN / 2 - 1;
            const typename Node::VecDim rightCoor = 
                (*((*e)->nodesBegin()+1))->giveCoordinates();
            Element* leftElement = (*e);
            Element* mirrorLeftElement = (*e);
            for ( unsigned g=0; g<numGhostVertices; ++g ) {

                // Create ghost element
                Element* ghostElement = new Element( );
                ghostElements_.push_back( ghostElement );
                VecEle2_ faceEle;
                faceEle[ 0 ] = leftElement;
                faceEle[ 1 ] = NULL;
                mapEleToFaceEle_.insert( std::make_pair( ghostElement, faceEle ) );

                // Attach ghost element
                faceEle[ 0 ] = this->leftNeighbourElement( leftElement );
                faceEle[ 1 ] = ghostElement;
                mapEleToFaceEle_.find( leftElement )->second = faceEle;

                // Create ghost node
                Node* ghostNode = new Node;
                ghostNode -> setId( nodes_.size() );
                nodes_.push_back( ghostNode );  // store node
                ghostNodes_.push_back( ghostNode );  // store node
                const typename Node::VecDim leftCoor = 
                    (*(mirrorLeftElement->nodesBegin()))->giveCoordinates();
                ghostNode->setCoordinates( 2.0*rightCoor - leftCoor );
//                ghostNode->setCoordinates( 4.0*rightCoor - 3.0*leftCoor );
                mirrorLeftElement = this->leftNeighbourElement( mirrorLeftElement );

                // Attach ghost node to ghost element
                *(ghostElement->nodesBegin()) = *(leftElement->nodesBegin()+1);
                *(ghostElement->nodesBegin()+1) = ghostNode;

                // advance to right
                leftElement = this->rightNeighbourElement( leftElement );
            }
        }
    }

    // loop elements to set their nodes
    for ( ElementConstIterator e=eb; e!=ee; ++e ) {
        // number of interpolating nodes of element
        const unsigned numNodesSN = Element::numNodesSN;

        // find left most element relevant for interpolation
        const unsigned degree = numNodesSN - 1;
        //const unsigned offSet = degree % 2;

        // local number of left vertex
        // local ID   0     1     2     3           d (degree)
        //            +.....+.....+-----+.....+.....+
        // global ID              k    k+1
        const unsigned leftVertexLocId = degree / 2;  // k
        const unsigned rightVertexLocId = leftVertexLocId + 1;  // k+1

        // initialise node pointers
        for ( unsigned n=0; n<numNodesSN; ++n )
            (*e)->setNode( n, NULL );

        // find interpolating nodes on left side
        Element* leftMost = *e;
        for ( int n=leftVertexLocId; n>=0; --n ) {
            (*e)->setNode( n, *(leftMost->nodesBegin()) );
            if ( this->leftNeighbourElement( leftMost ) != NULL )
                leftMost = this->leftNeighbourElement( leftMost );
            else
                break;
        }

        // find interpolating nodes on right side
        Element* rightMost = *e;
        for ( unsigned n=rightVertexLocId; n<numNodesSN; ++n ) {
            (*e)->setNode( n, *(rightMost->nodesBegin()+1) );
            if ( this->rightNeighbourElement( rightMost ) != NULL )
                rightMost = this->rightNeighbourElement( rightMost );
            else
                break;
        }
    }

#if 0
    for ( ElementConstIterator e=eb; e!=ee; ++e ) {
        std::cout << "element=" << *e
                  << " : left=" << this->leftNeighbourElement( *e )
                  << ", right=" << this->rightNeighbourElement( *e )
                  << std::endl;
    }
#endif

    // ghost element range
#if 0
    const ElementConstIterator geb = ghostElements_.begin();
    const ElementConstIterator gee = ghostElements_.end();
    for ( ElementConstIterator e=geb; e!=gee; ++e ) {
        std::cout << "ghostElement=" << *e
                  << " : left=" << this->leftNeighbourElement( *e )
                  << ", right=" << this->rightNeighbourElement( *e )
                  << std::endl;
    }
#endif

    return;
}

//------------------------------------------------------------------------------
//! Free all dynamically allocated memory 
template< typename ELEMENT >
beam::fem::Mesh< ELEMENT >::~Mesh()
{
    this->iterateOverLinks( corlib::deleteFunctor( ) );
    links_.clear( );

    // destroy ghosts
    mapEleToFaceEle_.clear();
    std::for_each( ghostElements_.begin(), ghostElements_.end(), corlib::deleteFunctor() );
    ghostElements_.clear();
    ghostNodes_.clear();


    return;
}

//------------------------------------------------------------------------------
//! Give node at provided node ID
template< typename ELEMENT >
typename ELEMENT::Node* beam::fem::Mesh< ELEMENT >::nodeFind( const unsigned nodeId ) const
{
    Node* node = NULL;
    for ( NodeConstIterator n=this->nodesBegin(); n!=this->nodesEnd(); ++n ) {
        if ( (*n)->giveId() == nodeId ) {
            node = *n;
            break;
        }
    }
    return node;
}

//------------------------------------------------------------------------------
//! Give node which is stored previous to given node
template< typename ELEMENT >
typename ELEMENT::Node* beam::fem::Mesh< ELEMENT >::nodeFindPrevious( const Node* theNode ) const
{
    // element range
    const ElementConstIterator eb = this->elementsBegin();
    const ElementConstIterator ee = this->elementsEnd();

    // initialise searched for node
    Node* previousNode = NULL;

    // loop elements and find element with left vertex being the node in question
    for ( ElementConstIterator e=eb; e!=ee; ++e ) {
        if ( *((*e)->nodesBegin()) == theNode ) {
            previousNode = *(this->leftNeighbourElement( *e )->nodesBegin());
            break;
        }
        if ( *((*e)->nodesBegin()+1) == theNode ) {
            previousNode = *((*e)->nodesBegin());
            break;
        }
    }

    // deliver
    return previousNode;
}

//------------------------------------------------------------------------------
//! Give node which is stored previous to given node
template< typename ELEMENT >
typename ELEMENT::Node* beam::fem::Mesh< ELEMENT >::nodeFindNext( const Node* theNode ) const
{
    // element range
    const ElementConstIterator eb = this->elementsBegin();
    const ElementConstIterator ee = this->elementsEnd();

    // initialise searched for node
    Node* nextNode = NULL;

    // loop elements and find element with left vertex being the node in question
    for ( ElementConstIterator e=eb; e!=ee; ++e ) {
        if ( *((*e)->nodesBegin()) == theNode ) {
            nextNode = *((*e)->nodesBegin()+1);
            break;
        }
        if ( *((*e)->nodesBegin()+1) == theNode ) {
            nextNode = *(this->rightNeighbourElement( *e )->nodesBegin()+1);
            break;
        }
    }

    // deliver
    return nextNode;
}

//------------------------------------------------------------------------------
template< typename ELEMENT >
ELEMENT* beam::fem::Mesh< ELEMENT >::leftNeighbourElement( Element* ele ) const
{
    // store resulting pointer to element
    Element* neighbourEle = NULL;

    // element in map?
    const MapEleToEle2ConstIter_ eleMap = mapEleToFaceEle_.find( ele );
    if ( eleMap != mapEleToFaceEle_.end() ) {
        neighbourEle = (eleMap->second)[ 0 ];
    }

    // done
    return neighbourEle;
}

//------------------------------------------------------------------------------
template< typename ELEMENT >
ELEMENT* beam::fem::Mesh< ELEMENT >::rightNeighbourElement( Element* ele ) const
{
    // store resulting pointer to element
    Element* neighbourEle = NULL;

    // element in map?
    const MapEleToEle2ConstIter_ eleMap = mapEleToFaceEle_.find( ele );
    if ( eleMap != mapEleToFaceEle_.end() ) {
        neighbourEle = (eleMap->second)[ 1 ];
    }

    // done
    return neighbourEle;
}

//------------------------------------------------------------------------------
//! Set list of links
template< typename ELEMENT >
template< typename SREADER >
void beam::fem::Mesh< ELEMENT >::setLinksOfSupports( const SREADER& suppReader )
{
    typedef typename SREADER::SupportConstIter SupportConstIter;
    const SupportConstIter lb = suppReader.begin();
    const SupportConstIter le = suppReader.end();
    for ( SupportConstIter el=lb; el!=le; ++el ) {
        const unsigned nodeId = el->first;
        Node * theNode = this->nodeFind( nodeId );
        FTL_VERIFY_DESCRIPTIVE( theNode != NULL,
                                "Node ID not found; cannot apply support" );
        const unsigned degree = Element::numNodesSN - 1;
        FTL_STATIC_ASSERT_MSG( ( degree % 2 == 1 ),
                               "Only odd spline degrees are allowed" );
        // store #theNode and previous and next node
        std::array< const Node *, degree > nodes;
        Node * prevNode = theNode;
        for ( unsigned p=0; p<(degree/2); ++p ) {
            prevNode = this->nodeFindPrevious( prevNode );
            nodes[ p ] = prevNode;
        }
        FTL_VERIFY_DESCRIPTIVE( prevNode != NULL,
                                "Node ID not found; cannot apply support" );
        nodes[ degree/2 ] = theNode;
        Node * nextNode = theNode;
        for ( unsigned p=(degree/2+1); p<degree; ++p ) {
            nextNode = this->nodeFindNext( nextNode );
            nodes[ p ] = nextNode;
        }
        FTL_VERIFY_DESCRIPTIVE( nextNode != NULL,
                                "Node ID not found; cannot apply support" );

        beam::fem::SupportsAtNode< Link, degree > supportsAtNode;
        supportsAtNode.set( links_, (el->second), nodes );
    }
}


//------------------------------------------------------------------------------
/** Apply a given functor to all links making use of std::for_each
 *  \tparam OP         Type of the element operation 
 *  \param[in,out] op  Specific operation applied to all element
 *  \retval OP         for_each returns the operator
 */
template< typename ELEMENT >
template< typename OP >
OP beam::fem::Mesh< ELEMENT >::iterateOverLinks( OP op )
{
    typename LinkVec_::iterator begin = links_.begin( );
    typename LinkVec_::iterator end   = links_.end( );
    return std::for_each( begin, end, op );
}


//------------------------------------------------------------------------------
//! Return SMF-like stream
template< typename ELEMENT >
std::ostream & beam::fem::Mesh< ELEMENT >::writeSmf( std::ostream & os,
                                                           bool withSmfHead ) const
{
    if ( withSmfHead ) {
        corlib::SmfHead smfHead;
        smfHead.write( convertShapeEnumToString( Element::myShape ),
                       Element::numNodes,
                       os );
    }

    os << ( this->numNodes() - this->numNodesGhost() ) << " "
       << this->numElements() << std::endl;

    const NodeConstIterator nb = this->nodesBegin();
    const NodeConstIterator ne = this->nodesEnd();
    for ( NodeConstIterator n=nb; n!=ne; ++n ) {
        if ( not this->nodeIsGhost( *n ) ) {
            //(*n)->writeId( os );
            (*n)->writeCoordinates( os );
        }
    }
        
    const ElementConstIterator eb = this->elementsBegin();
    const ElementConstIterator ee = this->elementsEnd();
    for ( ElementConstIterator e=eb; e!=ee; ++e ) {
        (*e)->writeNodeIndices( os );
    }

    return os;       
}
