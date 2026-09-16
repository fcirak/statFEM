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

//! @file   PartitionerBasic.hpp
//! @author Jakub Sistek
//! @date   3/2011

//------------------------------------------------------------------------------
/** Constructor given the total number of nodes, total number of elements, 
*  element shape, and number of nodes in one element.
*/
tools::input::partitioners::utils::
PartitionerBasic::PartitionerBasic( const unsigned numNodes, 
                                    const unsigned numElements,
                                    const unsigned numNodesPerElement )
  : numNodes_( numNodes ),
    numElements_( numElements ),
    numNodesPerElement_( numNodesPerElement ),
    elem2part_( numElements, -1 ),
    node2part_( numNodes, -1 ),
    isPartitioned_( false )
{ }

//------------------------------------------------------------------------------
/** Given mesh is divided into required number of parts. Partitionings of 
 *  elements and nodes is stored in private data.
 */
void tools::input::partitioners::utils::
PartitionerBasic::processPartitioning_( std::vector< std::vector<unsigned> > & connectivity )
{

#ifdef VERBOSE
    boost::timer::cpu_timer routineTimer;
    routineTimer.start();
#endif

    // Elements are partitioned, now derive partitioning of nodes
    // This has been inspired by METIS
    std::vector<unsigned>  nptr( numNodes_ + 1, 0 );

    // count elements connected to nodes
    for ( unsigned e = 0; e < numElements_; e ++ ) {
        // go throught nodes of this element
        for ( unsigned i = 0; i < numNodesPerElement_; i++ ) {
            nptr[ connectivity[e][i] ]++;  // += 1 per connected element
        }
    }

    // change counts to starts
    FTL_VERIFY( nptr[ numNodes_ ] == 0 );
    for ( unsigned i = 1; i < numNodes_; i++ )
        nptr[i] = nptr[i-1] + nptr[i];
    // shift it one back
    for ( unsigned i = numNodes_; i > 0; i-- )
        nptr[i] = nptr[i-1];
    // start from zero
    nptr[0] = 0;

    // create list containing entries of each element ID and recreate starts
    // nind = [ ... elem_0^k elem_1^k elem_2^k ... elem_0^j elem_1^j .... ]
    //              --------------------------     -----------------
    //                element indices               element indices
    //                attached to node k            attached to node j
    std::vector<unsigned>  nind( nptr[numNodes_] );
    for ( unsigned e = 0; e <  numElements_; e ++ ) {
        // go throught nodes of this element
        for ( unsigned i = 0; i < numNodesPerElement_; i++ ) {
            nind[ nptr[ connectivity[e][i] ]++ ] = e;  // inner ++ to be ready next occurrence of node
        }
    }
    for ( unsigned i = numNodes_; i > 0; i-- )  // rotate back as nptr list was shifted by accum. ++ operations
        nptr[i] = nptr[i-1];
    nptr[0] = 0;

#ifdef VERBOSE
    std::cout << " Time for creating list of elements indices for each node: " << routineTimer.elapsed().wall * 1e-9 << std::endl; 
    routineTimer.start();
#endif

    // assign nodes to subdomains based on element partitioning
    std::vector<unsigned> pwgts( numParts_, 0 );
    for ( unsigned i = 0; i < numNodes_; i++ ) {
        if ( nptr[ i+1 ] - nptr[ i ] > 0 ) {
            // only continue for connected nodes
            const unsigned me = elem2part_[ nind[ nptr[i] ] ];  // partition index based on first element attached to node
            // loop over number of elements attached to node
            bool inSamePart = true;
            for ( unsigned j = nptr[i] + 1; j < nptr[i + 1] and inSamePart; j++ ) {
                inSamePart = ( elem2part_[ nind[j] ] == me );  // check if element is in same partition
            }
            if ( inSamePart ) {  // if all elements connected to node are in same partition, then ...
                node2part_[i] = me;    // ... make node a node in partition
                pwgts[me]++;           // ... augment number of nodes in partition
            }
        }
    }

#ifdef VERBOSE
    std::cout << " Time for assigning interior nodes to subdomains: " << routineTimer.elapsed().wall * 1e-9 << std::endl; 
    routineTimer.start();
#endif

    // maximal inbalance 3 %
    const unsigned maxpwgt = static_cast<unsigned>( 1.03 * numNodes_ / numParts_ ); // intentional conversion to int

    for ( unsigned i = 0; i < numNodes_; i++ ) {
        if ( node2part_[i] == -1 ) { // Assign the boundary element (node2part_ was initialised with -1)
            if ( nptr[i+1] - nptr[i] > 0 ) {  // only continue for connected nodes

                // fill nbrind with partition indices attached to node
                // fill nbrwgt with how many elements of one partition are attached
                std::vector<unsigned> nbrind;  nbrind.reserve( 200 );
                std::vector<unsigned> nbrwgt;  nbrwgt.reserve( 200 );

                for ( unsigned j = nptr[i]; j < nptr[i+1]; j++ ) {  // loop over node's elements
                    const unsigned me = elem2part_[ nind[j] ];

                    const typename std::vector<unsigned>::iterator meIt =
                        std::find( nbrind.begin(), nbrind.end(), me );
                    if ( meIt != nbrind.end() ) {
                        nbrwgt[ meIt - nbrind.begin() ] += 1;
                    }
                    else {
                        nbrind.push_back( me );
                        nbrwgt.push_back( 1 );
                    }
                }

                // Try to assign it first to the domain with most things in common 
                const unsigned maxpos = std::distance( nbrwgt.begin(),
                                                     std::max_element( nbrwgt.begin(), nbrwgt.end() ) );

                if ( pwgts[nbrind[maxpos]] < maxpwgt ) {
                    node2part_[i] = nbrind[maxpos];
                }
                else {
                    // If that fails, assign it to a light domain
                    node2part_[i] = nbrind[0];
                    for ( unsigned j = 0; j < nbrind.size(); j++ ) {
                        if ( pwgts[nbrind[j]] < maxpwgt ) {
                            node2part_[i] = nbrind[j];
                            break;
                        }
                    }
                }

                pwgts[ node2part_[i] ]++;  // increment number of nodes in partition
            }
            else {
                node2part_[i] = -2; // mark node which is not present in any element
            }
        }
    }

#ifdef VERBOSE
    std::cout << " Time for assigning boundary nodes to subdomains: " << routineTimer.elapsed().wall * 1e-9 << std::endl; 
    routineTimer.start();
#endif

    // check node partitioning
    FTL_VERIFY_DESCRIPTIVE( std::find( node2part_.begin( ), node2part_.end( ), -1 ) == node2part_.end( ),
                            "error in partitioning of nodes" );
    // count number of nodes not connected to elements
    const unsigned numDangling = std::count( node2part_.begin( ), node2part_.end( ), -2 );
    if ( numDangling > 0 ) {
        std::cout << "WARNING: " << numDangling << " dangling nodes detected in the mesh, removing them." << std::endl;
    }

    const unsigned numNodesCorrected = numNodes_ - numDangling;

    // erase element at nodes list
    nptr.clear();
    std::vector<unsigned>().swap( nptr );
    nind.clear();
    std::vector<unsigned>().swap( nind );

    // erase neighbours
    pwgts.clear();
    std::vector<unsigned>().swap( pwgts );

    // map nodes into new numbering based on their partition
    nodesOld2new_.resize( numNodes_, -1 );
    nodesNew2old_.resize( numNodesCorrected, -1 );

    // create from node2part a list of assigned nodes for each subdomain
    std::vector<unsigned> numNodesInParts( numParts_ );
    for ( unsigned iNode  = 0; iNode < numNodes_; iNode++ ) {
        const int iPart = node2part_[iNode];  // partition of node
        if ( iPart != -2 ) {
            numNodesInParts[iPart]++;
        }
    }
    // reserve memory for the lists to avoid reallocations
    std::vector< std::vector<unsigned> > listNodesInParts( numParts_ );
    for ( unsigned iPart  = 0; iPart < numParts_; iPart++ ) {
        listNodesInParts[iPart].reserve( numNodesInParts[iPart] );
    }
    // now create the lists
    for ( unsigned iNode  = 0; iNode < numNodes_; iNode++ ) {
        const int iPart = node2part_[iNode];  // partition of node
        if ( iPart != -2 ) {
            listNodesInParts[iPart].push_back( iNode );
        }
    }
    // get rid of former node2part_
    node2part_.clear();
    std::vector<int>().swap( node2part_ );

    // where nodes of subdomains start
    sub2nodeStart_.resize( numParts_ + 1, 0 );
    // new index
    unsigned newInd = 0;  // holds step-by-step new global index of nodes
    for ( unsigned iPart = 0; iPart < numParts_; iPart++ ) {
        unsigned numNodesUnique = 0;  // nodes in partition
        for ( unsigned ilocNode  = 0; ilocNode < numNodesInParts[iPart]; ilocNode++ ) {
            unsigned iNode = listNodesInParts[iPart][ilocNode];
            nodesOld2new_[ iNode ]  = newInd;
            nodesNew2old_[ newInd ] = iNode;
            // increase counters
            newInd++;
            numNodesUnique++;
        }
        sub2nodeStart_[ iPart + 1 ] = sub2nodeStart_[ iPart ] + numNodesUnique;
    }
    // clear memory
    for ( unsigned iPart  = 0; iPart < numParts_; iPart++ ) {
        listNodesInParts[iPart].clear( );
    }
    listNodesInParts.clear( );
    std::vector< std::vector<unsigned> >().swap( listNodesInParts );
    numNodesInParts.clear( );
    std::vector<unsigned>().swap( numNodesInParts );

    // mark unconnected nodes
    std::remove_copy_if( node2part_.begin(), node2part_.end(),
                         nodesOld2new_.begin(), boost::bind( std::not_equal_to<unsigned>(), _1, -2 ) );

#ifdef VERBOSE
    std::cout << " Time for creating new numbering of nodes: " << routineTimer.elapsed().wall * 1e-9 << std::endl; 
    routineTimer.start();
#endif

    // checks
    FTL_VERIFY_DESCRIPTIVE( std::count( nodesOld2new_.begin( ),nodesOld2new_.end( ), -2 ) == (int)numDangling,
                            "error in number of dangling nodes in mapping" );
    FTL_VERIFY_DESCRIPTIVE( std::find( nodesOld2new_.begin( ), nodesOld2new_.end( ), -1 ) == nodesOld2new_.end( ),
                            "error in permutation, missing indices" );
    FTL_VERIFY_DESCRIPTIVE( std::find( nodesNew2old_.begin( ), nodesNew2old_.end( ), -1 ) == nodesNew2old_.end( ),
                            "error in permutation, missing indices " );

    // map elements into new numbering based on their partition
    elementsOld2new_.resize( numElements_, -1 );
    elementsNew2old_.resize( numElements_, -1 );

    // create from elem2part a list of assigned elements for each subdomain
    std::vector<unsigned> numElemsInParts( numParts_ );
    for ( unsigned iElem  = 0; iElem < numElements_; iElem++ ) {
        const unsigned iPart = elem2part_[iElem];  // partition of node
        numElemsInParts[iPart]++;
    }
    // reserve memory for the lists to avoid reallocations
    std::vector< std::vector<unsigned> > listElemsInParts( numParts_ );
    for ( unsigned iPart  = 0; iPart < numParts_; iPart++ ) {
        listElemsInParts[iPart].reserve( numElemsInParts[iPart] );
    }
    // now create the lists
    for ( unsigned iElem  = 0; iElem < numElements_; iElem++ ) {
        const unsigned iPart = elem2part_[iElem];  // partition of node
        listElemsInParts[iPart].push_back( iElem );
    }

    // where elements of subdomains start
    sub2elementStart_.resize( numParts_ + 1, 0 );
    // new index
    newInd = 0;  // holds step-by-step new global index of elements
    for ( unsigned iPart  = 0; iPart < numParts_; iPart++ ) {
        unsigned numElementsUnique = 0;
        for ( unsigned iLocEl = 0; iLocEl < numElemsInParts[iPart]; iLocEl++ ) {
            unsigned iEl = listElemsInParts[iPart][iLocEl];
            elementsOld2new_[ iEl ]    = newInd;
            elementsNew2old_[ newInd ] = iEl;
            // increase counters
            newInd++;
            numElementsUnique++;
        }
        sub2elementStart_[ iPart + 1 ] = sub2elementStart_[ iPart ] + numElementsUnique;
    }
    // clear memory
    for ( unsigned iPart  = 0; iPart < numParts_; iPart++ ) {
        listElemsInParts[iPart].clear( );
    }
    listElemsInParts.clear( );
    std::vector< std::vector<unsigned> >().swap( listElemsInParts );
    numElemsInParts.clear( );
    std::vector<unsigned>().swap( numElemsInParts );

#ifdef VERBOSE
    std::cout << " Time for creating new numbering of elements: " << routineTimer.elapsed().wall * 1e-9 << std::endl; 
    routineTimer.start();
#endif

    // check number
    FTL_VERIFY_DESCRIPTIVE( sub2elementStart_[numParts_] == numElements_,
                            "error in element partition: %d %d \n", 
                            sub2elementStart_[numParts_] - 1, numElements_ );

    // check mappings
    FTL_VERIFY_DESCRIPTIVE( std::find( elementsOld2new_.begin( ), elementsOld2new_.end( ), -1 ) == elementsOld2new_.end( ),
                            "error in permutation, missing indices" );
    FTL_VERIFY_DESCRIPTIVE( std::find( elementsNew2old_.begin( ), elementsNew2old_.end( ), -1 ) == elementsNew2old_.end( ),
                            "error in permutation, missing indices " );

    isPartitioned_ = true;

    // Determine overlaps
    buildCluster_ = ( overlapLayers_ > 0 );

#ifdef VERBOSE
    if ( buildCluster_ ) // derive partitions from local nodes
        std::cout << " Creating division based on overlapping clusters." <<  std::endl;
    else                 // derive partitions from local elements
        std::cout << " Creating non-overlapping division." <<  std::endl;
#endif

    // update connectivity to new node numbering
    PartitionerBasic::updateConnectivityNumbering_( connectivity, false );

#ifdef VERBOSE
    std::cout << "Time for updating connectivity: " << routineTimer.elapsed().wall * 1e-9 << std::endl;
    routineTimer.start();
#endif

    // transpose the updated connectivity - i.e. make list of elements at nodes
    // this is needed only for generating clusters
    if ( buildCluster_ ) {
        transposeConnectivity_( connectivity );
    }
#ifdef VERBOSE
    std::cout << "Time for transposing connectivity: " << routineTimer.elapsed().wall * 1e-9 << std::endl;
    routineTimer.start();
#endif

}

//------------------------------------------------------------------------------
/** Update the connectivity matrix to new node numbering. Must be called after createPartition.
 *  This is done in-place.
 *  If backward is not given, it is assumed that connectivity in original numbering is translated to new node numbering.
 *  If backward is TRUE, connectivity in NEW numbering is translated into original numbering.
 */
void tools::input::partitioners::utils::
PartitionerBasic::updateConnectivityNumbering_( std::vector< std::vector<unsigned> > & connectivity,
                                                bool backward )
{
    // check data
    FTL_VERIFY_DESCRIPTIVE( isPartitioned_ == true,
                            "New numbering is not ready. Call createPartition first." );
    FTL_VERIFY_DESCRIPTIVE( connectivity.size( ) == static_cast<unsigned>( numElements_ ),
                            "Dimension mismatch." );

//    unsigned nodeIndNew;
    for ( unsigned e = 0; e < static_cast<unsigned>( numElements_ ); e ++ ) {
        // go throught nodes of this element
            
        // use this function if we have it
        std::vector<unsigned> element;
        this -> updateListNumbering( connectivity[e], element, backward );
        connectivity[e] = element;
        //for ( unsigned i = 0; i < connectivity[e].size(); i++ ) {
        //    unsigned nodeInd = connectivity[e][i];

        //    // get new number for node
        //    if ( backward ) {
        //        nodeIndNew = nodesNew2old_ [ nodeInd ];
        //    }
        //    else {
        //        nodeIndNew = nodesOld2new_ [ nodeInd ];
        //    }
        //    connectivity[e][i] = nodeIndNew;
        //}
    }
}
//------------------------------------------------------------------------------
/** Update the list of global nodes to new node numbering. Used with boundary conditions given by lists lf nodes.
 *  This is done in-place.
 *  If backward is not given, it is assumed that a list in original numbering is translated to new node numbering.
 *  If backward is TRUE, list in NEW numbering is translated into original numbering.
 */
void tools::input::partitioners::utils::
PartitionerBasic::updateListNumbering( std::vector<unsigned> & listOld,
                                       std::vector<unsigned> & listNew,
                                       bool backward = false ) const
{
                                         

    // check data
    FTL_VERIFY_DESCRIPTIVE( isPartitioned_ == true,
                            "New numbering is not ready. Call createPartition first." );

    // go throught the list
    for ( unsigned i = 0; i < listOld.size(); i++ ) {
        unsigned nodeInd = listOld[i];

        int nodeIndNew;

        // get new number for node
        if ( backward ) {
            nodeIndNew = nodesNew2old_ [ nodeInd ];
        }
        else {
            nodeIndNew = nodesOld2new_ [ nodeInd ];
        }

        if ( nodeIndNew >= 0 ) {
            listNew.push_back( nodeIndNew );
        }
        else {
            std::cout << " Ignoring a node from list: "<< nodeInd << std::endl;
        }
    }
}
//------------------------------------------------------------------------------
/** Transpose the connectivity matrix, i.e. generate list of elements at nodes.
 */
void tools::input::partitioners::utils::
PartitionerBasic::transposeConnectivity_( const std::vector< std::vector<unsigned> > & connectivity ) {
    
    // check data
    FTL_VERIFY_DESCRIPTIVE( connectivity.size( ) == static_cast<unsigned>( numElements_ ),
                            "Dimension mismatch." );

    // prepare space for list of elements at nodes
    elems4nodes_.resize( numNodes_ );

    for ( unsigned e = 0; e < static_cast<unsigned>( numElements_ ); e ++ ) {
        // go throught nodes of this element
        for ( unsigned i = 0; i < connectivity[e].size(); i++ ) {
            unsigned nodeInd = connectivity[e][i];

            // mark this element 
            elems4nodes_[nodeInd].insert( e );
        }
    }
}

//------------------------------------------------------------------------------
/** For a given mesh and privately stored partitioning, generate a submesh of subdomain. 
 *  Connectivity table is assumed updated to new node numbering.
 */
void tools::input::partitioners::utils::
PartitionerBasic::getSubMesh( const std::vector< ublas::bounded_vector<double,3> > & coordinates, 
                              const std::vector< std::vector<unsigned> > & elements,
                              const unsigned subIndex, 
                              std::vector< std::vector<unsigned> > & elementsSub,
                              std::vector<unsigned> & elementsGlobal4local,
                              std::vector<unsigned> & elementsOwner,
                              std::vector< ublas::bounded_vector<double,3> > & coordinatesSub, 
                              std::vector<unsigned> & nodesGlobal4local,
                              std::vector<unsigned> & nodesOwner )
{

    //timing of this routine
    boost::timer::cpu_timer partTimer;

    // set of global element indices included in local part
    std::set<unsigned> elementsGlobal4localSet;

    // set of global node indices included in local part
    std::set<unsigned> nodesGlobal4localSet;

    // build cluster if asked to - i.e. make sure that all elements in which local nodes are contained 
    // are involved in the subdomain
    if ( buildCluster_ ) {
        FTL_VERIFY_DESCRIPTIVE( elems4nodes_.size() != 0,
                                "Computing clusters but missing info on elements at nodes. \n" );

        for ( unsigned indNode = sub2nodeStart_[ subIndex ]; 
              indNode < sub2nodeStart_[ subIndex + 1 ]; 
              indNode++ ) {
            nodesGlobal4localSet.insert( indNode );
        }
    }
    else {

        // Go through all elements and make note on nonoverlapping partition
        for ( unsigned e = 0; e < elements.size(); e++ ) {
            // If element belongs to this subdomain
            if ( elem2part_[e] == subIndex ) {

                // make a note to include this element
                elementsGlobal4localSet.insert( e );
            }
        }
        // deduce local nodes for elements included in the set
        for ( std::set<unsigned>::const_iterator iterE = elementsGlobal4localSet.begin(); 
              iterE != elementsGlobal4localSet.end();
              ++iterE ) {
            // initialize new element to store local numbering
            unsigned ie = *iterE;
            for ( unsigned n = 0; n < elements[ ie ].size( ); n ++ ) {
                const unsigned nodeNumGlob = elements[ ie ][n];
                // add this node to mapping
                nodesGlobal4localSet.insert( nodeNumGlob );
            }
        }
    }
#ifdef VERBOSE
    std::cout << "---getSubMesh - time for initial set: " << partTimer.elapsed().wall * 1e-9 << std::endl;
#endif
    partTimer.start();

    for ( unsigned iOverlap = 0; iOverlap < overlapLayers_; iOverlap++ ) {
        // add layer of elements to overlap
#ifdef VERBOSE
        std::cout << " Creating overlap layer " <<  iOverlap + 1 <<  std::endl;
#endif

        std::set<unsigned>::const_iterator iterN = nodesGlobal4localSet.begin();
        for ( ; iterN != nodesGlobal4localSet.end( ); ++iterN ) { 
            unsigned indNode = *iterN;

            for ( std::set<unsigned>::const_iterator iterE = elems4nodes_[indNode].begin(); 
                  iterE != elems4nodes_[indNode].end();
                  ++iterE ) {
                elementsGlobal4localSet.insert( *iterE );
            }
        }
#ifdef VERBOSE
        std::cout << "---getSubMesh - time for elements set: " << partTimer.elapsed().wall * 1e-9 << std::endl;
#endif
        partTimer.start();

        // deduce local nodes for elements included in the set
        for ( std::set<unsigned>::const_iterator iterE = elementsGlobal4localSet.begin(); 
              iterE != elementsGlobal4localSet.end();
              ++iterE ) {
            // initialize new element to store local numbering
            unsigned ie = *iterE;
            for ( unsigned n = 0; n < elements[ ie ].size( ); n ++ ) {
                const unsigned nodeNumGlob = elements[ ie ][n];
                // add this node to mapping
                nodesGlobal4localSet.insert( nodeNumGlob );
            }
        }
#ifdef VERBOSE
        std::cout << "---getSubMesh - time for nodes set: " << partTimer.elapsed().wall * 1e-9 << std::endl;
#endif
        partTimer.start();
    }

    // turn the global 4 local map of elements to vector
    elementsGlobal4local.resize( elementsGlobal4localSet.size( ) );
    std::copy(  elementsGlobal4localSet.begin( ), elementsGlobal4localSet.end( ),
                elementsGlobal4local.begin( ));
    elementsGlobal4localSet.clear();

    // turn the global 4 local map of nodes to vector
    nodesGlobal4local.resize( nodesGlobal4localSet.size( ) );
    std::copy(  nodesGlobal4localSet.begin( ), nodesGlobal4localSet.end( ),
                nodesGlobal4local.begin( ));
    nodesGlobal4localSet.clear( );

    // for the sake of speed of global to local map, convert the vector to map 
    std::map<unsigned,unsigned> nodesGlobal4localMap;
    for ( unsigned i = 0; i < nodesGlobal4local.size(); ++i ) {
        nodesGlobal4localMap.insert( std::make_pair( nodesGlobal4local[i], i ) );
    }

#ifdef VERBOSE
    std::cout << "---getSubMesh - time for converting sets to vectors: " << partTimer.elapsed().wall * 1e-9 << std::endl;
#endif
    partTimer.start();

    // extract local elements
    elementsSub.reserve( elementsGlobal4local.size() );
    elementsOwner.reserve( elementsGlobal4local.size() );
    for ( std::vector<unsigned>::const_iterator iterE = elementsGlobal4local.begin(); 
          iterE != elementsGlobal4local.end();
          ++iterE ) {

        unsigned elGlobIndex = *iterE;
        std::vector<unsigned> element = elements[ elGlobIndex ];
        // now element contains node numbers in original numbering, change it to subdomain
        for ( unsigned n = 0; n < element.size( ); n ++ ) {
            const unsigned nodeNumGlob = element[n];
            const unsigned nodeNumLoc  = nodesGlobal4localMap[nodeNumGlob];
            element[n] = nodeNumLoc;
        }

        // insert element in local numbering
        elementsSub.push_back( element );

        // insert subdomain that owns this element
        //elementsOwner.push_back( PartitionerBasic::getSubIndex_( elementsOld2new_[ elGlobIndex ], 
        //                                                         sub2elementStart_ ) );
        
        // OR
        
        // decide if element is active
        unsigned isActive; 
        if ( PartitionerBasic::getSubIndex_( elementsOld2new_[ elGlobIndex ], sub2elementStart_ ) == subIndex ) {
            isActive = 1;
        }
        else {
            isActive = 0;
        }
        // insert activity flag
        elementsOwner.push_back( isActive );
    }
    nodesGlobal4localMap.clear();
#ifdef VERBOSE
    std::cout << "---getSubMesh - time for extracting elements: " << partTimer.elapsed().wall * 1e-9 << std::endl;
#endif
    partTimer.start();


    // extract array of local coordinates
    coordinatesSub.reserve( nodesGlobal4local.size( ) );
    nodesOwner.reserve( nodesGlobal4local.size( ) );
    for ( std::vector<unsigned>::const_iterator it = nodesGlobal4local.begin( );
          it != nodesGlobal4local.end( ); ++it ) {
        unsigned nodeNumNew = *it;
        unsigned nodeNumOld = nodesNew2old_[ nodeNumNew ];
        coordinatesSub.push_back( coordinates[ nodeNumOld ] );

        // insert subdomain that owns this node
        //nodesOwner.push_back( PartitionerBasic::getSubIndex_( nodeNumNew, sub2nodeStart_ ) );
        
        // OR
        
        // decide if node is active
        unsigned isActive; 
        if ( getSubIndex_( nodeNumNew, sub2nodeStart_ ) == subIndex ) {
            isActive = 1;
        }
        else {
            isActive = 0;
        }
        // insert activity flag
        nodesOwner.push_back( isActive );
    }
#ifdef VERBOSE
    std::cout << "---getSubMesh - time for extracting nodes: " << partTimer.elapsed().wall * 1e-9 << std::endl;
#endif
    partTimer.start();


    return;
}

//------------------------------------------------------------------------------
/** Print the partitions inside the object.
 */
void tools::input::partitioners::utils::
PartitionerBasic::printPartitions( ) const
{
    // print partitioning of elements
    std::cout << " partition of elements, size " <<  elem2part_.size( ) << std::endl;
    for ( unsigned e = 0; e < elem2part_.size( ); e++ ) {
        std::cout << e << "  " << elem2part_[e] << std::endl;
    }
    std::cout << " partition of nodes, size " <<  node2part_.size( ) << std::endl;
    for ( unsigned n = 0; n < node2part_.size( ); n++ ) {
        std::cout << n << "  " << node2part_[n] << std::endl;
    }
    std::cout << " renumbering of nodes, size " <<  nodesOld2new_.size( ) << std::endl;
    for ( unsigned n = 0; n < nodesOld2new_.size( ); n++ ) {
        std::cout << n << "  " << nodesOld2new_[n] << std::endl;
    }
    std::cout << " starts of node indices for subdomains, size " <<  sub2nodeStart_.size( ) << std::endl;
    for ( unsigned n = 0; n < sub2nodeStart_.size( ); n++ ) {
        std::cout << n << "  " << sub2nodeStart_[n] << std::endl;
    }

    return;
}

//------------------------------------------------------------------------------
//! return part number in linear partition given by starts for an entry
unsigned tools::input::partitioners::utils::
PartitionerBasic::getSubIndex_( const unsigned index,
                                const std::vector <unsigned> & starts )
{
    
    FTL_VERIFY_DESCRIPTIVE( starts.size() > 0,
                            "Array starts not ready to be used by getSubIndex \n" );
        
    for ( unsigned iSub = 0; iSub < starts.size() - 1; iSub ++ ) {
        if ( starts[ iSub ] <= index and starts[ iSub + 1 ] > index )
            return iSub;
    }
        
    FTL_VERIFY_DESCRIPTIVE( false, "Error: Index %d out of range for getSubIndex. \n", index );
    exit(-1);
}
