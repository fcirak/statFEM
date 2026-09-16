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

//! @file   PartitionerRCB.hpp
//! @author Jakub Sistek, Matija Kecman, Fehmi Cirak
//! @date   5/2011

//------------------------------------------------------------------------------
/** \brief Constructor given the total number of nodes, total number of elements, 
 *  element shape, and number of nodes in one element.
 */
tools::input::partitioners::rcb::mesh::
PartitionerRCB::PartitionerRCB( const unsigned numNodes, 
                                const unsigned numElements,
                                const unsigned numNodesPerElement )
    : PartitionerBasic_( numNodes, numElements, numNodesPerElement )
{ }

//------------------------------------------------------------------------------
/** Given mesh is divided into required number of parts. Partitionings of 
 *  elements and nodes is stored in private data.
 */
void tools::input::partitioners::rcb::mesh::
PartitionerRCB::createPartition( std::vector< std::vector<unsigned> > & connectivity,
                                 const std::vector< ublas::bounded_vector<double,3> > & coords,
                                 const int numParts, const int overlapLayers )
{

    // make a note on the parts
    numParts_ = numParts;
    overlapLayers_ = overlapLayers;

    // preconditions:
    FTL_VERIFY_DESCRIPTIVE( connectivity.size() == static_cast<unsigned>(numElements_), 
                            "number of element mismatch" );
    //--------------------------------------------------------------------------
    // for convenience
    namespace algo =  tools::input::partitioners::rcb::algo;
    typedef algo::RCBNode            RCBNode;
    typedef algo::RCBNodeContainer   RCBNodeContainer;

    // compute the centres of the elements
    std::vector< ublas::bounded_vector<double,3> > elementCenters;
    tools::input::partitioners::utils::computeElementCenters( connectivity, coords, elementCenters );

    //! Do the RCB partitioning
    algo::performRCB( elementCenters, numParts, elem2part_ );

    // process partitioning
    this -> processPartitioning_( connectivity );

    return;
}
