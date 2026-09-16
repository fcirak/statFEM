// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//
//                 Computational Structural Mechanics Lab
//                         University of Cambridge
//
// This software is copyrighted and all rights are retained by the University of
// Cambridge. All use, disclosure, and/or reproduction of any part not expressly
// authorised by the University of Cambridge is prohibited. (C) 2010
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

//! @file Parallel.ipp

//------------------------------------------------------------------------------
template< typename DRIVERSEQ >
unsigned solid::driver::ParallelBasic<DRIVERSEQ>::generateDofs()
{
    //create global vector of dof starts
    std::map<int,int> numNodeDofsGlobal;
    std::pair<int,int> numsU = this -> generateDofStarts( commAll_, numNodeDofsGlobal );

    // find number of global dofs
    const unsigned numTotalNodes  = numsU.first;
    const unsigned numTotalDofs   = numsU.second;

    // number dofs with the given array
    this -> dofNumberFromStartsFun( numNodeDofsGlobal );
    //free memory needed
    numNodeDofsGlobal.clear();

    // inform user about basic parameters
    this -> accessOutStream_()
        << "number of processes " << this ->getNumProcessors() << std::endl
        << "total number of DOFs:          " << numTotalDofs << std::endl;
    this -> accessOutStream_().flush();

    return numTotalDofs;
}

//------------------------------------------------------------------------------
template< typename ISTREAM, typename OSTREAM >
void solid::driver::partitionWithRCB( const enum corlib::shape myShape,
                                      std::istream & smf,
                                      std::vector<ISTREAM *> & lists,
                                      const unsigned overlapLayers,
                                      const int numSub,
                                      const int iSub,
                                      OSTREAM & subSmf,
                                      std::vector<OSTREAM *> & subLists )
{
    std::vector< ublas::bounded_vector<double,3> > coordinates;
    std::vector< std::vector<unsigned> >           elements;
    enum corlib::shape                             eleShape;
        
    smf.seekg( 0 );  // rewind smf stream
    tools::input::smf2vtu::readSmfFile( smf, eleShape, coordinates, elements );
    FTL_VERIFY( eleShape == myShape );

    const unsigned numNodes           = coordinates.size();
    const unsigned numElements        = elements.size();
    const unsigned numNodesPerElement = elements[0].size();

    // initialize partitioner
    namespace partitioners = tools::input::partitioners;
    partitioners::rcb::mesh::PartitionerRCB partitioner( numNodes, numElements, 
                                                         numNodesPerElement );

    // partition the mesh by RCB
    partitioner.createPartition( elements, coordinates, numSub, overlapLayers );
    
    // get submesh
    std::vector< ublas::bounded_vector<double,3> > coordinatesSub;
    std::vector<unsigned>                          nodesGlobal4local;
    std::vector<unsigned>                          nodesOwner;

    std::vector< std::vector<unsigned> >           elementsSub;
    std::vector<unsigned>                          elementsGlobal4local;
    std::vector<unsigned>                          elementsOwner;

    // get data for the subdomain from partitioner 
    partitioner.getSubMesh( coordinates, elements, iSub, 
                            elementsSub, elementsGlobal4local, elementsOwner,
                            coordinatesSub, nodesGlobal4local, nodesOwner );

    // write the submesh file
    tools::input::smf2vtu::writeSmfFile( subSmf, eleShape,
                                         coordinatesSub, elementsSub, 
                                         &nodesGlobal4local, 
                                         &elementsGlobal4local, 
                                         &nodesOwner, 
                                         &elementsOwner );

    typename std::vector<ISTREAM *>::iterator ins  = lists.begin();
    typename std::vector<ISTREAM *>::iterator insE = lists.end();
    for ( ; ins != insE; ++ins ) {

        // get corresponding list
        const int listIndex = std::distance( lists.begin(), ins );

        FTL_VERIFY( (*ins)->good( ) );
        
        unsigned length;
        *(*ins) >> length;
        std::vector< unsigned > oneListGlobal;
        for ( unsigned i = 0; i < length; i++ ) {
            unsigned nodeId;
            *(*ins) >> nodeId;
            oneListGlobal.push_back( nodeId );
        }

        // update list to new node numbering
        std::vector< unsigned > oneListGlobalNew;
        partitioner.updateListNumbering( oneListGlobal, oneListGlobalNew );

        // make sure the list is sorted
        std::sort( oneListGlobalNew.begin(), oneListGlobalNew.end() );

        // debug
        //std::copy( nodeListGlobal.begin(), nodeListGlobal.end(), 
        //           std::ostream_iterator<unsigned>( std::cout, " " ) );
        //std::cout << " \n ";

        // prepare array for intersection
        std::vector<unsigned> locNodesIntersect( oneListGlobal.size() + nodesGlobal4local.size() ) ;

        // localize the list of nodes to subdomain by comparison of list ( *itList ) with nodesGlobal4local
        std::vector<unsigned>::iterator posEnd;
        posEnd = std::set_intersection( oneListGlobal.begin(), oneListGlobal.end(), 
                                        nodesGlobal4local.begin(), nodesGlobal4local.end(), 
                                        locNodesIntersect.begin() );

        // export the set of fixed nodes 
        FTL_VERIFY( (subLists[listIndex])->good( ) );
            
        *(subLists[listIndex]) << std::distance( locNodesIntersect.begin(), posEnd ) << std::endl;
        std::copy( locNodesIntersect.begin(), posEnd, 
                   std::ostream_iterator<unsigned>( *(subLists[listIndex]), " \n" ) );

    }



    // set-up sub-mesh based on psmf
    subSmf.seekg( 0 );

}
