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

#ifndef gshell_fem_concentratedload_h
#define gshell_fem_concentratedload_h

#include <Eigen/Core>

//==============================================================================
// declarations
namespace gshell {
    namespace fem {

        template< typename NODE >
        class ConcentratedLoad;

    }
}

//==============================================================================
// definitions

//------------------------------------------------------------------------------
//! This object applies concentrated loads to nodes.
//!
//! The input file looks like
//!\verbatim
//! numberOfConcentratedLoads
//! coordX1 coordY1 coordZ1  loadX1 loadY1 loadZ1
//! coordX2 coordY2 coordZ2  loadX2 loadY2 loadZ2
//! ...     ...     ...      ...    ...    ...
//!\endverbatim
//!
//! \tparam   NODE          Type of node
template< typename NODE >
class gshell::fem::ConcentratedLoad :
    public std::unary_function< NODE *, void >
{
public:
    typedef NODE                                             Node;

    typedef typename Node::VecDim                            VecDim;
    typedef typename Node::VecDof                            VecDof;

    //! nodeId and a local co-ordinate
    typedef std::multimap< Node *, VecDim >                  LoadPoints;
    typedef typename LoadPoints::const_iterator              LoadPointsConstIter;

public:
    //! Empty constructor
    ConcentratedLoad() : tolerance_( -1. ), factor_( 0. ) { }

    //! Constructor
    //!
    //! \param[in]  nodesBegin        Node begin iterator
    //! \param[in]  nodesEnd          Node end iterator
    //! \param[in]  inputFileStream   Stream to input file
    //! \param[in]  tolerance         Tolerance to find control nodes based on co-ordinates
    template< typename NODEITER >
    ConcentratedLoad( const NODEITER & nodesBegin,
                      const NODEITER & nodesEnd,
                      std::ifstream & inputFileStream,
                      const double & tolerance,
                      const bool stickToTolerance = false ) :
        tolerance_( tolerance ), factor_( 0. )
    {
        this->read( nodesBegin, nodesEnd, inputFileStream, tolerance,
                    stickToTolerance );
    }

    //! Destructor
    ~ConcentratedLoad( ) { }

    //! Read file and collect nodes with concentrated loads
    //!
    //! \param[in]  nodesBegin        Node begin iterator
    //! \param[in]  nodesEnd          Node end iterator
    //! \param[in]  inputFileStream   Stream to input file
    //! \param[in]  tolerance         Tolerance to find control nodes based on co-ordinates
    template< typename NODEITER >
    void read( const NODEITER & nodesBegin,
               const NODEITER & nodesEnd,
               std::ifstream & inputFileStream,
               const double & tolerance,
               const bool stickToTolerance = false );

    //! Set load factor
    void setFactor( const double factor ) { factor_ = factor; return; }

    //! Apply load to nodal forces
    void operator()( Node * node ) const;

private:
    //! Tolerance for identifying nodes
    double tolerance_;
    //! points to apply load
    LoadPoints loadPoints_;
    //! LOad factor
    double factor_;
};

//------------------------------------------------------------------------------
// Read file and collect nodes with concentrated loads
template< typename NODE >
template< typename NODEITER >
void gshell::fem::ConcentratedLoad< NODE >::read(
    const NODEITER & nodesBegin,
    const NODEITER & nodesEnd,
    std::ifstream & inputFileStream,
    const double & tolerance,
    const bool stickToTolerance
    )
{
    typedef NODEITER NodeConstIterator;
    tolerance_ = tolerance;

    // get number of loaded nodes
    unsigned numConcentratedLoaded = 0;
    inputFileStream >> numConcentratedLoaded;

    // get loaded nodes, determine attached elements and local co-ordinates
    unsigned lostPoints = 0;
    for ( unsigned m=0; m<numConcentratedLoaded; ++m ) {
        // node co-ordinate
        VecDim loadCoord;
        for ( unsigned d=0; d<Node::dim; ++d )
            inputFileStream >> loadCoord[ d ];
        VecDim loadVec;
        for ( unsigned d=0; d<Node::dim; ++d )
            inputFileStream >> loadVec[ d ];

        // loop nodes to find node to apply load
        Node * loadNode = NULL;
        const NodeConstIterator nBegin = nodesBegin;
        const NodeConstIterator nEnd = nodesEnd;
        for ( NodeConstIterator n=nBegin; n!=nEnd and not loadNode; ++n ) {
            const VecDim coord = (*n)->giveCoordinates( );
            const double diff = ( coord - loadCoord ).norm( );
            if ( diff < tolerance_ ) {
                loadNode = *n;
            }
        }
        //FTL_VERIFY_DESCRIPTIVE( loadNode,
        //                        "A node with co-ordinate (%g,%g,%g)"
        //                        " could be found using a tolerance of %g\n",
        //                        loadCoord( 0 ), loadCoord( 1 ), loadCoord( 2 ),
        //                        tolerance_ );

        // point with tolerance failed, look for closest
        if ( not loadNode and not stickToTolerance ) {
            loadNode = *nBegin;
            double loadDist = ( loadCoord - loadNode->giveCoordinates( ) ).norm( );
            for ( NodeConstIterator n = nBegin; n != nEnd; ++n ) {
                double dist = ( loadCoord - (*n)->giveCoordinates( ) ).norm( );
                if ( dist < loadDist ) {
                    loadNode = *n;
                    loadDist = dist;
                }
            }
        }

        // store node
        if ( loadNode )
            loadPoints_.insert( std::make_pair( loadNode, loadVec ) );
        else
            lostPoints += 1;

    }
    // check if all nodes were found
    FTL_VERIFY_DESCRIPTIVE( loadPoints_.size( ) == numConcentratedLoaded - lostPoints,
                            "%d != %d\n", loadPoints_.size( ), numConcentratedLoaded );

    return;
}

//------------------------------------------------------------------------------
// Apply load to nodal forces
template< typename NODE >
void gshell::fem::ConcentratedLoad< NODE >::operator()(
    Node * node
    ) const
{
    // find concentrated loads of node
    const std::pair< LoadPointsConstIter, LoadPointsConstIter > myLoads = 
        loadPoints_.equal_range( node );
    // loop node's concentrated loads
    if ( myLoads.first != myLoads.second ) {
        for ( LoadPointsConstIter l=myLoads.first; l!=myLoads.second; ++l ) {
            FTL_VERIFY( node == (*l).first );
            // get load vector
            VecDof loadVec; loadVec.setZero( );
            loadVec.head( Node::dim ) = factor_ * (*l).second;
            // apply load vector
            node->addToForce( loadVec );   // negative RHS
        }
    }
    return;
}

#endif
