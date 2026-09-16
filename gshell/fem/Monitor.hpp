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

#ifndef gshell_fem_monitor_h
#define gshell_fem_monitor_h

//==============================================================================
// declarations
namespace gshell {
    namespace fem {

        template< typename NODE, typename ELEMENT >
        class Monitor;

    }
}

//==============================================================================
// definitions

//------------------------------------------------------------------------------
//! Monitor interpolated (limit) quantities, such as displacement, on
//! surface locations. The locations are chosen by providing a reference
//! co-ordinate of a control node in a file. The monitored quantities
//! are written to a file #out_.
//!
//! The input file looks like
//!\verbatim
//! numberOfMonitoredNodes
//! coordX1 coordY1 coordZ1
//! coordX2 coordY2 coordZ2
//! ...     ...     ...
//!\endverbatim
//!
//! Since the interpolation is attached to the elements, an element is found
//! of which the monitored node is a vertex. Its local element-wise
//! index in the element is stored in #monPoints_ together with the
//! element pointer.
//!
//! \tparam   NODE          Type of node
//! \tparam   ELEMENT       Type of element
template< typename NODE, typename ELEMENT >
class gshell::fem::Monitor :
    public std::unary_function< ELEMENT *, void >
{
public:
    typedef NODE                                             Node;
    typedef ELEMENT                                          Element;

    typedef typename Node::VecDim                            VecDim;
    typedef typename Element::VecLDim                        VecLDim;

    //! Read operation function type
    typedef std::function< VecDim( const Element *, const unsigned ) > ReadFun;

    //! NodeId and a local element-wise node ID
    typedef std::pair< Node *, unsigned >                    NodeAndLIndex;
    //! Map from elements to monitored local element-wise node ID
    typedef std::multimap< Element *, NodeAndLIndex >        MonPoints;
    typedef typename MonPoints::const_iterator               MonPointsConstIter;

public:

    //! Constructor
    //!
    //! \tparam     MESH             Type of mesh
    //! \param[in]  mesh              A mesh
    //! \param[in]  readFun           Element function to read monitored quantity
    //! \param[in]  inputFileStream   Stream to input file
    //! \param[in]  outputFileStream  Stream to output file
    //! \param[in]  tolerance         Tolerance to find control nodes based on co-ordinates
    //! \param[in]  offset            Small offset to determine element-wise parametric co-ordinate
    template< typename MESH >
    Monitor( const MESH & mesh,
             ReadFun readFun,
             std::ifstream & inputFileStream,
             std::ofstream & outputFileStream,
             const double & tolerance );

    //! Constructor
    //!
    //! \tparam     MESH             Type of mesh
    //! \param[in]  mesh              A mesh
    //! \param[in]  readFun           Element function to read monitored quantity
    //! \param[in]  inputFileName     Name to input file
    //! \param[in]  outputFileName    Name to output file
    //! \param[in]  tolerance         Tolerance to find control nodes based on co-ordinates
    //! \param[in]  offset            Small offset to determine element-wise parametric co-ordinate
    template< typename MESH >
    Monitor( const MESH & mesh,
             const ReadFun readFun,
             const std::string & inputFileName,
             const std::string & outputFileName,
             const double & tolerance );

    //! Copy constructor
    //!
    //! \param[in]  old               Monitor to copy
    Monitor( const Monitor< NODE, ELEMENT > & old );

    //! Destructor
    ~Monitor( );

    //! Print quantity for elements at stored local co-ordinates
    void operator()( Element * elem,
                     const unsigned step,
                     const double time );

    //! Write monitored nodes
    std::ostream & write( std::ostream & os ) const;

private:
    //! Set-up
    //!
    //! \tparam     MESH             Type of mesh
    //! \param[in]  mesh              A mesh
    //! \param[in]  inputFileStream   Stream to input file
    template< typename MESH >
    void setup_( const MESH & mesh,
                 std::ifstream & inputFileStream );

    //! write mapping of co-ordinate to node Id
    void writeMap_( const std::pair< Element *, NodeAndLIndex > & mp ) const;

    //! write header for data
    void writeHead_( ) const
    {
        *out_ << "#" << std::endl;
        *out_ << "# step time nodeId uX uY uZ" << std::endl;
        return;
    }

private:
    //! Tolerance for identifying nodes
    double tolerance_;

    //! Output file name
    std::string outputFileName_;
    //! Output stream
    std::ofstream * out_;

    //! Function to get quantity at monitored points
    ReadFun readFun_;

    //! points to monitor
    MonPoints monPoints_;
};

//------------------------------------------------------------------------------
template< typename NODE, typename ELEMENT >
template< typename MESH >
gshell::fem::Monitor< NODE, ELEMENT >::Monitor(
    const MESH & mesh,
    ReadFun readFun,
    std::ifstream & inputFileStream,
    std::ofstream & outputFileStream,
    const double & tolerance
    ) :
    tolerance_( tolerance ),
    outputFileName_( "" ),
    out_( &outputFileStream ),
    readFun_( readFun )
{
    this->setup_( mesh, inputFileStream );
    return;
}

//------------------------------------------------------------------------------
template< typename NODE, typename ELEMENT >
template< typename MESH >
gshell::fem::Monitor< NODE, ELEMENT >::Monitor(
    const MESH & mesh,
    const ReadFun readFun,
    const std::string & inputFileName,
    const std::string & outputFileName,
    const double & tolerance
    ) :
    tolerance_( tolerance ),
    outputFileName_( outputFileName ),
    out_( NULL ),
    readFun_( readFun )
{
    out_ = new std::ofstream( outputFileName.c_str( ) );
    FTL_VERIFY_DESCRIPTIVE( out_->is_open( ),
                            "Could not open file %s\n", inputFileName.c_str( ) );
    std::ifstream inputFileStream( inputFileName.c_str( ) );
    FTL_VERIFY_DESCRIPTIVE( inputFileStream.is_open( ),
                            "Could not open file %s\n", outputFileName.c_str( ) );

    this->setup_( mesh, inputFileStream );

    inputFileStream.close( );
    return;
}

//------------------------------------------------------------------------------
template< typename NODE, typename ELEMENT >
gshell::fem::Monitor< NODE, ELEMENT >::Monitor(
    const Monitor< NODE, ELEMENT > & old
    ) :
    tolerance_( old.tolerance_ ),
    outputFileName_( old.outputFileName_ ),
    out_( old.out_ ),
    readFun_( old.readFun_ ),
    monPoints_( old.monPoints_ )
{
    
    if ( outputFileName_ != "" ) {
        out_ = new std::ofstream( outputFileName_.c_str( ), std::ios_base::app );
    }
}

//------------------------------------------------------------------------------
template< typename NODE, typename ELEMENT >
gshell::fem::Monitor< NODE, ELEMENT >::~Monitor( )
{
    if ( outputFileName_ != "" and out_ ) {
        out_->close( );
        delete out_;
        out_ = NULL;
    }
    return;
}

//------------------------------------------------------------------------------
template< typename NODE, typename ELEMENT >
template< typename MESH >
void gshell::fem::Monitor< NODE, ELEMENT >::setup_(
    const MESH & mesh,
    std::ifstream & inputFileStream
    )
{
    typedef MESH                                   Mesh;
    typedef typename Mesh::NodeConstIterator       NodeConstIterator;
    typedef typename Mesh::ElementConstIterator    ElementConstIterator;
    static_assert(( std::is_same<Node,typename Mesh::NodeType>::value ));
    static_assert(( std::is_same<Element,typename Mesh::ElementType>::value ));


    // checks
    FTL_VERIFY_DESCRIPTIVE( tolerance_ > 0.0, "Tolerance must be strictly positive\n" );
    static_assert( ELEMENT::myShape == corlib::TRIANGLE or
                         ELEMENT::myShape == corlib::QUADRILATERAL );

    // get number of monitored nodes
    unsigned numMonitored = 0;
    inputFileStream >> numMonitored;

    // get monitored nodes, determine attached elements and local co-ordinates
    for ( unsigned m=0; m<numMonitored; ++m ) {
        // node co-ordinate
        VecDim monCoord;
        for ( unsigned d=0; d<Node::dim; ++d )
            inputFileStream >> monCoord[ d ];

        // loop nodes to find node to monitor
        Node * monNode = NULL;
        const NodeConstIterator nBegin = mesh.nodesBegin( );
        const NodeConstIterator nEnd = mesh.nodesEnd( );
        for ( NodeConstIterator n=nBegin; n!=nEnd and not monNode; ++n ) {
            const VecDim coord = (*n)->giveCoordinates( );
            const double diff = ( coord - monCoord ).norm( );
            if ( diff < tolerance_ ) monNode = *n;
        }
        //FTL_VERIFY_DESCRIPTIVE( monNode,
        //                        "A node with co-ordinate (%g,%g,%g)"
        //                        " could not be found using a tolerance of %g\n",
        //                        monCoord( 0 ), monCoord( 1 ), monCoord( 2 ),
        //                        tolerance_ );

        // point with tolerance failed, look for closest
        if ( not monNode ) {
            monNode = *nBegin;
            double monDist = ( monCoord - monNode->giveCoordinates( ) ).norm( );
            for ( NodeConstIterator n = ++(mesh.nodesBegin( )); n != nEnd; ++n ) {
                double dist = ( monCoord - (*n)->giveCoordinates( ) ).norm( );
                if ( dist < monDist ) {
                    monNode = *n;
                    monDist = dist;
                }
            }
        }

        // node ID
        const unsigned monNodeId = monNode->giveId( );
        std::cout << "Monitor : (x,y,z)=" << monCoord << " on nodeId=" << monNodeId << std::endl;

        // loop elements in mesh
        const ElementConstIterator eBegin = mesh.elementsBegin( );
        const ElementConstIterator eEnd = mesh.elementsEnd( );
        for ( ElementConstIterator e = eBegin; e != eEnd; ++e ) {

            // collect node IDs of element vertices
            std::vector< unsigned > vertexIndices( Element::numVertices );
            (*e)->giveNodeIndices( vertexIndices.begin( ) );
            // look for node ID
            const std::vector< unsigned >::iterator ni =
                std::find( vertexIndices.begin( ), vertexIndices.end( ), monNodeId );
            // go on if element has the vertex in question
            if ( ni != vertexIndices.end( ) ) {
                // store
                const unsigned localIndex = std::distance( vertexIndices.begin(), ni );
                NodeAndLIndex p = std::make_pair( monNode, localIndex );
                std::pair< Element *, NodeAndLIndex > mp = std::make_pair( *e, p );
                monPoints_.insert( mp );

                // write mapping of co-ordinate to node Id
                this->writeMap_( mp );
                
                // quit loop
                break;
            }

        }
    }
    // check if all nodes were found
    FTL_VERIFY_DESCRIPTIVE( monPoints_.size( ) == numMonitored, 
                            "%d != %d\n", monPoints_.size( ), numMonitored );

    // write header for data
    this->writeHead_( );

    return;
}

//------------------------------------------------------------------------------
template< typename NODE, typename ELEMENT >
void gshell::fem::Monitor< NODE, ELEMENT >::writeMap_(
    const std::pair< Element *, NodeAndLIndex > & mp
    ) const
{
    const Element * monElement = mp.first;
    const Node * monNode = mp.second.first;
    const unsigned localIndex = mp.second.second;
    
    *out_ << "#";
    *out_ << " coord=";
    for ( unsigned d=0; d<Node::dim; ++d )
        *out_ << " " << monNode->giveCoordinates( )( d );
    *out_ << ", nodeId=" << monNode->giveId( );
    *out_ << ", elemNodes=";
    std::vector< unsigned > elemNodeIndices( Element::numVertices );
    monElement->giveNodeIndices( elemNodeIndices.begin( ) );
    std::copy( elemNodeIndices.begin( ), elemNodeIndices.end( ),
               std::ostream_iterator< unsigned >( *out_, " " ) );
    *out_ << ", localIndex=" << localIndex;
    *out_ << std::endl;
    return;
}

//------------------------------------------------------------------------------
//! Print quantity for elements at stored local co-ordinates
//!
//! \param[in]  elem       Element pointer
template< typename NODE, typename ELEMENT >
void gshell::fem::Monitor< NODE, ELEMENT >::operator()(
    Element * elem,
    const unsigned step,
    const double time
    )
{
    // find monitored points of element
    const std::pair< MonPointsConstIter, MonPointsConstIter > myPoints = 
        monPoints_.equal_range( elem );
    // loop element's monitored points
    if ( myPoints.first != myPoints.second ) {
        for ( MonPointsConstIter p = myPoints.first; p != myPoints.second; ++p ) {
            // adjust format
            *out_ << std::scientific << std::setprecision( 9 ) << std::right;
            FTL_VERIFY( elem == (p->first) );
            // write time step and time
            *out_ << std::setw( 4 ) << step << " " << std::setw( 16 ) << time;
            const Node * monNode = (p->second).first;
            // write node ID
            const unsigned nodeId = monNode->giveId( );
            *out_ << " " << std::setw( 5 ) << nodeId;
            // interpolate quantity
            const unsigned localIndex = (p->second).second;
            const VecDim quan = readFun_( elem, localIndex );
            for ( unsigned d=0; d<Element::dim; ++d )
                *out_ << " " << std::setw( 16 ) << quan( d );
            *out_ << std::endl;
        }
    }
    return;
}

//------------------------------------------------------------------------------
//! Write list of monitored nodes (convenience call)
//!
//! \param[in,out]  os     Output stream
//! \return                Output stream
template< typename NODE, typename ELEMENT >
std::ostream & gshell::fem::Monitor< NODE, ELEMENT >::write(
    std::ostream & os
    ) const
{
    const MonPointsConstIter mBegin = monPoints_.begin( );
    const MonPointsConstIter mEnd   = monPoints_.end( );
    for ( MonPointsConstIter m = mBegin; m != mEnd; ++m ) {
        os << "Coord=" << (m->second).first->giveCoordinates( );
        os << ", NodeId=" << (m->second).first->giveId( );
        os << ", Element=" << (m->first);
        os << ", LocalIndex=" << (m->second).second;
        os << std::endl;
    }
    return os;
}

#endif
