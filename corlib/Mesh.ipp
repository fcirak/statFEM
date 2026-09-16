// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//
//                 Computational Structural Mechanics Lab
//                         University of Cambridge
//
// This software is copyrighted and all rights are retained by the University of
// Cambridge. All use, disclosure, and/or reproduction of any part not expressly
// authorised by the University of Cambridge is prohibited. (C) 2010
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

//! @file Mesh.ipp

//------------------------------------------------------------------------------
//! system includes
#include <set>
#include <limits>
#include <iterator>

//! boost includes
#include <boost/tuple/tuple.hpp>
#include <boost/bind.hpp>

//! corlib includes
#include <corlib/SmfHead.hpp>
#include <corlib/misc.hpp>
#include <corlib/verify.hpp>
#include <corlib/EvaluateField.hpp>

//------------------------------------------------------------------------------
/** Constructor given a stream to a simple mesh file, see #read() for details
 *
 * \param smf Input stream
 */
template<typename ELEMENT>
corlib::Mesh<ELEMENT>::Mesh( std::istream & smf )
    : isCopy_( false ) // set flag that says the mesh is just 
                         // a copy of other mesh and should not dealocate pointers
{
    // read the file
    this -> read( smf );
}

//------------------------------------------------------------------------------
/** Add nodes and elements given a stream to a simple mesh file.
 *
 *  This file must have the following format:
 *  <pre>
 *    NN   NE  
 *    x1  y1  z1   
 *    ....         
 *    xNN yNN zNN  
 *    v1_1  v1_2  ... v1_NPE    
 *    ...                      
 *    vNE_1 vNE_2 ... vNE_NPE  
 *  </pre>
 * where NN is the number of nodes, NE the number of elements, and NPE refers to
 * the number of nodes per element. Note that the latter number is fixed at 
 * compile time.
 *
 * \param smf Input stream
 */
template<typename ELEMENT>
void corlib::Mesh<ELEMENT>::read( std::istream & smf )
{
    FTL_VERIFY( smf.good() );
    
    // validate input file
    const corlib::shape elementShape  = Element::myShape;
    const unsigned numNodesPerElement = Element::numNodes;
    SmfHead smfHead;
    smfHead.readValidated( smf, elementShape, numNodesPerElement );

    // start reading data
    unsigned numNodes, numElements;
    //! read number of nodes and elements
    smf >> numNodes >> numElements;
    smf.ignore( std::numeric_limits<std::streamsize>::max(), '\n' );

    //! read nodes
    nodes_.reserve( numNodes );
    for ( unsigned n = 0; n < numNodes; n ++ ) {
        Node* theNode = new Node;
        // Important: 'readSelf' might set a new ID for parallel computations
        theNode -> setId( n );
        theNode -> readSelf( smf );
        nodes_.push_back( theNode );
    }

    //! read elements
    elements_.reserve( numElements );
    for ( unsigned e = 0; e < numElements; e ++ ) {
        Element* theElement = new Element;
        // split elements by ends of lines
        theElement -> readSelf( smf, nodes_ );
        smf.ignore( std::numeric_limits<std::streamsize>::max(), '\n' );
        elements_.push_back( theElement );
    }

    FTL_VERIFY( smf.good() );

    return;
}

//------------------------------------------------------------------------------
/** Constructor using iterators of a container with coordinates and a container
 *  with a connectivity.
 *  \tparam IT1    Type of iterator to coordinate container
 *  \tparam IT2    Type of iterator to connectivity container
 *  \param[in]  vb Begin of coordinate container
 *  \param[in]  ve End   of coordinate container
 *  \param[in]  eb Begin of connectivity container
 *  \param[in]  ee End   of connectivity container
 */
template<typename ELEMENT>
template< typename IT1, typename IT2 > 
corlib::Mesh<ELEMENT>::Mesh( IT1 vb, IT1 ve, IT2 eb, IT2 ee )
    : isCopy_( false ) // set flag that says the mesh is just 
                       // a copy of other mesh and should not dealocate pointers
{
    FTL_VERIFY( ( vb != ve ) and ( eb != ee ) );

    // set up nodes
    const unsigned numNodes = 
        static_cast<unsigned>( std::distance( vb, ve ) );
    nodes_.reserve( numNodes );

    unsigned nodeId = 0;
    for( ; vb != ve; ++ vb ) {
        Node* theNode = new Node;
        theNode -> setId( nodeId++ );
        theNode -> setCoordinates( *vb );
        nodes_.push_back( theNode );
    }

    // set up elements
    const unsigned numElements = 
        static_cast<unsigned>( std::distance( eb, ee ) );
    elements_.reserve( numElements );

    for( ; eb != ee; ++ eb ) {
        Element* theElement = new Element;

        typename Element::NodeIterator first = theElement -> nodesBegin();
        typename Element::NodeIterator  last = theElement -> nodesEnd();
        for ( unsigned v = 0; first != last; ++first, v ++ ) {
            const unsigned vIndex = (*eb)[ v ];
            *first = nodes_[ vIndex ];
        }
        elements_.push_back( theElement );
    }
}

//------------------------------------------------------------------------------
/** Constructor with iterators with element pointers
 *  This constructor is used for duplicating meshes.
 *  \tparam     ITEP      Type of iterator to element pointer container
 *  \param[in]  itepb     Begin of element pointers container
 *  \param[in]  itepe     End of element pointers container
 *  \param[in]  nodeComp  Optional comparator how to order nodes in mesh's node vector
 */
template<typename ELEMENT>
template< typename ITEP >
corlib::Mesh<ELEMENT>::Mesh( ITEP itepb, ITEP itepe,
                             boost::function<bool(Node*,Node*)> nodeComp )
    : isCopy_( true ) // set flag that says the mesh is just 
                      // a copy of other mesh and should not dealocate pointers
{
    FTL_VERIFY( itepb != itepe );

    // set up elements
    unsigned numElements = std::distance( itepb, itepe );
    elements_.resize( numElements );
    std::copy( itepb, itepe, elements_.begin( ) );

    // create set of pointers of really used nodes
    typedef boost::function<bool(Node*,Node*)> NodeComp;
    std::set<Node*,NodeComp> nodesSet( nodeComp );

    ElementConstIterator iteb = elements_.begin( );
    ElementConstIterator itee = elements_.end( );
    for ( ; iteb != itee; ++iteb ) {
        // insert the nodes of the element to the set
        typename Element::NodeConstIterator itnb = (*iteb) -> nodesBegin();
        typename Element::NodeConstIterator itne = (*iteb) -> nodesEnd();

        nodesSet.insert ( itnb, itne );
    }

    // set up nodes
    unsigned numNodes = nodesSet.size( );
    nodes_.resize( numNodes );
    std::copy( nodesSet.begin( ), nodesSet.end( ), nodes_.begin( ) );
}

//------------------------------------------------------------------------------
/** Free all dynamically allocated memory                                     */
template<typename ELEMENT>
corlib::Mesh<ELEMENT>::~Mesh( )
{
    // if #isCopy_==true only a shallow destruction occurs
    // a deep removal is carried out for #isCopy_==false (the default)

    if ( not isCopy_ ) {
        this -> iterateOverElements( corlib::deleteFunctor( ) );
    }
    elements_.clear( );
    if ( not isCopy_ ) {
        this -> iterateOverNodes( corlib::deleteFunctor( ) );
    }
    nodes_.clear( );

    return;
}

//------------------------------------------------------------------------------
/** Apply a given functor to all nodes making use of std::for_each
 *  \tparam OP         Type of the node operation 
 *  \param[in,out] op  Specific operation applied to all nodes
 *  \retval OP         for_each returns the operator
 */
template<typename ELEMENT>
template< typename OP >
OP corlib::Mesh<ELEMENT>::iterateOverNodes( OP op ) 
{
    return std::for_each( nodes_.begin(), nodes_.end(), op );
}

//------------------------------------------------------------------------------
/** Apply a given functor to all elements making use of std::for_each
 *  \tparam OP         Type of the element operation 
 *  \param[in,out] op  Specific operation applied to all element
 *  \retval OP         for_each returns the operator
 */
template<typename ELEMENT>
template< typename OP >
OP corlib::Mesh<ELEMENT>::iterateOverElements( OP op )
{
    return std::for_each( elements_.begin(), elements_.end(), op );
}

//------------------------------------------------------------------------------
/** Apply a given functor only to those nodes which have a given predicate.
 *  Given the predicate 'p', the functor 'op' will only be applied to those
 *  nodes for which holds 'p(X) = true'.
 *  \tparam OP          Type of functor to be applied to nodes
 *  \tparam PRED        Type of predicate functor
 *  \param[in,out] op   The functor to be applied
 *  \param[in] p        The predicate
 *  \retval op          Returns the given functor
 */
template<typename ELEMENT>
template<typename OP, typename PRED> 
OP corlib::Mesh<ELEMENT>::iterateOverNodesWithPredicate( OP op, PRED p )
{
    return corlib::for_each_if( nodes_.begin(), nodes_.end(), op, p );
}

//------------------------------------------------------------------------------
/** Apply a given functor only to those elements which have a given predicate.
 *  Given the predicate 'p', the functor 'op' will only be applied to those
 *  elements for which holds 'p(X) = true'.
 *  \tparam OP          Type of functor to be applied to elements
 *  \tparam PRED        Type of predicate functor
 *  \param[in,out] op   The functor to be applied
 *  \param[in] p        The predicate
 *  \retval op          Returns the given functor
 */
template<typename ELEMENT>
template<typename OP, typename PRED> 
OP corlib::Mesh<ELEMENT>::iterateOverElementsWithPredicate( OP op, PRED p )
{
    return corlib::for_each_if( elements_.begin(), elements_.end(), op, p );
}


//------------------------------------------------------------------------------
/** Write an SMF file containing the mesh vertex coordinates and the mesh
 *  topology. 
 *  \param[in,out] smf  Stream to write to
 */
template<typename ELEMENT>
std::ostream & corlib::Mesh<ELEMENT>::writeSmf( std::ostream & smf ) 
{
    FTL_VERIFY( smf.good() );

    // write appropriate header
    corlib::SmfHead smfHead;
    smfHead.setElementShape(     Element::myShape );
    smfHead.setElementNumPoints( Element::numNodes );
    smfHead.write( smf );
    // number of nodes and elements
    smf << this -> numNodes() << "  " << this -> numElements() << "\n";
    //--------------------------------------------------------------------------
    // write node coordinates (always 3D)
    for ( NodeConstIterator first = this -> nodesBegin();
          first != this -> nodesEnd(); ++ first ) {
        typename Node::VecDim x = (*first) -> giveCoordinates();
        for ( unsigned d = 0; d < Node::dim; d ++ )
            smf << x[d] << "  ";
        for ( unsigned d = Node::dim; d < 3; d ++ )
            smf << "0  ";
        smf << "\n";
    }
    //--------------------------------------------------------------------------
    // write element connectivity
    typedef std::ostream_iterator<unsigned>                            Out;
    typedef boost::function<unsigned (const typename ELEMENT::Node*) > GetId;
    // iterate over all elements
    for( ElementConstIterator first = this -> elementsBegin();
         first != this -> elementsEnd(); ++ first ) {
        // functor to insert id into ostream iterator
        corlib::detail_::GetNodalDatum<GetId,Out> 
            getId( boost::bind( &ELEMENT::Node::giveId, _1), 
                   std::ostream_iterator<unsigned>( smf, "  " ) );
        // apply to all nodes
        (*first) -> iterateOverNodes( getId );
        smf << "\n";
    }

    return smf;
}
