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

//! @file   BasicRCB.ipp
//! @author Matija Kecman, Fehmi Cirak
//! @date   2011
//! @todo   - remove the sizetypedefs by switching to iterators

//------------------------------------------------------------------------------
/** This function performs the Recursive Coordinate Bisection (RCB) Algorithm.
 *
 *  \param[in]  points     Points (just coordinates) to be partitioned
 *  \param[in]  numParts   Number of parts (subdomains) to create
 *  \param[out] point2part Vector containing the part number for each point
 */
void tools::input::partitioners::rcb::algo::
performRCB( const std::vector< ublas::bounded_vector<double,3> > & points,
            const int & numParts,
            std::vector<int> & point2part )
{
    typedef std::vector< ublas::bounded_vector<double,3> > VecVecD3;

    // construct RCBPoints
    RCBPointContainer RCBPoints;
    VecVecD3::const_iterator pIter = points.begin();
    for ( ; pIter != points.end(); ++pIter ) {
        const int pId = std::distance( points.begin(), pIter );
        RCBPoint * point = new RCBPoint( *pIter, pId );
        RCBPoints.push_back( point );
    }

    // Create RCB tree
    RCBNode * head = new RCBNode( RCBPoints, 0 );

    // trivial partition
    if ( numParts == 1 )
        std::cout << " Warning: Dividing into a single part ... Is that what you want? " << std::endl; 

    // loop over number of parts
    for ( int part = 1; part != numParts; ++part ) {

        // fill leaves vector with leaf nodes
        RCBNodeContainer leaves;
        algo::findLeaves( head, leaves );

        // find leaf to be partitioned
        RCBNode * leafToPartition = leaves[0];
        double   maxBox = 0.0; // maximal axis-aligned length
        unsigned maxDim = 0;   // axis of maximal length
        for ( unsigned i = 0; i < leaves.size(); ++ i ) {
            double   leafBox = 0.0;
            unsigned leafDim = 0;

            algo::getMaxBoxAndMaxDim( leaves[i], leafBox, leafDim );

            if ( leafBox > maxBox ) {
                maxBox = leafBox;
                maxDim = leafDim;
                leafToPartition = leaves[i];
            }
        }

        // partition that leaf
        algo::partitionNode( leafToPartition, maxBox, maxDim );

    }

    // extract leaves
    RCBNodeContainer leaves;
    algo::findLeaves( head, leaves );
    FTL_VERIFY_DESCRIPTIVE( static_cast<int>( leaves.size() ) == numParts,
                            "with leaves.size()=%d, numParts=%d\n", leaves.size(), numParts );

    //! Extract partitioning
    algo::extractPartitioning( leaves, point2part );

    // clean up
    delete head;  // use pointers on RCBPoints
    std::for_each( RCBPoints.begin(), RCBPoints.end(), corlib::deleteFunctor() );

    return;
}

//----------------------------------------------------------------------------
/** \brief Finds leaves in RCB tree
 */
void tools::input::partitioners::rcb::algo::findLeaves( RCBNode * node, 
                                                        RCBNodeContainer & leaves )
{
    if( leaf( node ) ){
        assert( node != NULL );
        leaves.push_back( node );
    }
    else {
        findLeaves( node->childL, leaves );
        findLeaves( node->childR, leaves );
    }

    return;
}

//------------------------------------------------------------------------------
/** \brief Leaf checker. Make this inline.
 */
bool tools::input::partitioners::rcb::algo::leaf( const RCBNode * node )
{
    return ( !node->childL );
}

//------------------------------------------------------------------------------
/** \brief Function to compute the *longest* axis and the dimension in which 
 *  this occurs of a partition defined by begin and end.
 */
void tools::input::partitioners::rcb::algo::getMaxBoxAndMaxDim( const RCBNode * node,
                                                                double & maxBox, 
                                                                unsigned & maxDim )
{
    // smallest/largest center component in each coordinate direction
    double lowerCoord[3];
    double upperCoord[3];

    const RCBPointContainer::const_iterator begin = (node->data).begin();
    const RCBPointContainer::const_iterator end   = (node->data).end();
    assert( begin != end );

    // initialise lowerCoord and upperCoord
    for ( unsigned i = 0; i < 3; ++i ) {
        lowerCoord[i] = (*begin)->center[i];
        upperCoord[i] = lowerCoord[i];
    }

    //! Iterate over data stored at node to position lowerCoord and upperCoord
    for ( RCBPointContainer::const_iterator iter = begin; iter != end; ++iter ) {
        for ( unsigned i = 0; i < 3; ++i ) {
            lowerCoord[i] = std::min( lowerCoord[i], (*iter)->center[i] );
            upperCoord[i] = std::max( upperCoord[i], (*iter)->center[i] );
        }
    }

    // find maxBox and maxDim
    maxDim = std::numeric_limits<unsigned>::max();
    for ( unsigned i = 0; i < 3; ++i ) {
        double box = upperCoord[i] - lowerCoord[i];
        if ( box > maxBox ) {
            maxBox = box;
            maxDim = i;
        }
    }
    assert( maxDim < 3 );

    return;
}

//------------------------------------------------------------------------------
/** \brief Partition the given Node
 */
void tools::input::partitioners::rcb::algo::partitionNode( RCBNode * parent, 
                                                           double & maxBox, 
                                                           unsigned & maxDim )
{
    // parent data
    const RCBPointContainer::const_iterator begin = (parent->data).begin();
    const RCBPointContainer::const_iterator end   = (parent->data).end();

    // threshold determines which child the parent data will belong to 
    double threshold = std::numeric_limits<double>::max();
    for ( RCBPointContainer::const_iterator iter = begin; iter != end; ++iter ) {
        threshold = std::min( threshold, ((*iter)->center)[maxDim] );
    }
    threshold += maxBox / 2.;

    // data for new child nodes
    RCBPointContainer dataChildL, dataChildR;

    // dividing parent RCB nodes according to coordinate direction maxDim
    for ( RCBPointContainer::const_iterator iter = begin; iter != end; ++iter ){
        if ( ((*iter)->center)[maxDim] < threshold )
            dataChildL.push_back( *iter );
        else
            dataChildR.push_back( *iter );
    }

    // create two new children
    parent->childL = new RCBNode( dataChildL, parent->level + 1 );
    parent->childR = new RCBNode( dataChildR, parent->level + 1 );

    return;
}

//------------------------------------------------------------------------------
/** This function generates a point to subdomain id mapping from a container of
 *  RCB nodes.
 *
 *  \param[in]  leaves     Container of RCB leaves from which to extract partitioning
 *  \param[out] point2part Vector containing the part number for each point
 */
void
tools::input::partitioners::rcb::algo::
extractPartitioning( const RCBNodeContainer & leaves,
                     std::vector<int> & pointToPart )
{
    // allocate and invalidate pointToPart vector
    unsigned numPoints = 0;
    for (RCBNodeConstIterator it = leaves.begin(); it != leaves.end(); ++it) {
        numPoints += (*it)->data.size();
    }
    pointToPart.resize( numPoints );
    std::fill(pointToPart.begin(), pointToPart.end(), -1);

    // generate point to part mapping
    int partId = 0;
    for (RCBNodeConstIterator itl = leaves.begin(); itl != leaves.end(); ++itl, ++partId) // parts
        for (RCBPointConstIterator itp = (*itl)->data.begin(); itp != (*itl)->data.end(); ++itp) // points
            pointToPart[(*itp)->pointId] = partId;

    // post-conditions: all points are allocated to a part
    assert(std::find(pointToPart.begin(), pointToPart.end(), -1) == pointToPart.end());

    return;
}
