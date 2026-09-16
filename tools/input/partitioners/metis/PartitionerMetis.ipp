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

//! @author Jakub Sistek
//! @date   3/2011


//------------------------------------------------------------------------------
/** Constructor given the total number of nodes, total number of elements, 
*  element shape, and number of nodes in one element.
*/
tools::input::partitioners::metis::PartitionerMetis::
PartitionerMetis( const unsigned numNodes, 
                  const unsigned numElements,
                  const corlib::shape eleShape,
                  const unsigned numNodesPerElement )
    : PartitionerBasic_( numNodes, numElements, numNodesPerElement )
{
    // determine type of elements
    if ( eleShape == corlib::TRIANGLE ) {
        elemTypeM_ = 1; // triangle 
        numNodesPerLinElementM_ = 
            corlib::ShapeTraits<corlib::TRIANGLE>::numVertices;
        nCommonNodes_ = 2;
    }
    else if ( eleShape == corlib::QUADRILATERAL ) {
        elemTypeM_ = 4; // quadrilateral
        numNodesPerLinElementM_ = 
            corlib::ShapeTraits<corlib::QUADRILATERAL>::numVertices;
        nCommonNodes_ = 2;
    }
    else if ( eleShape == corlib::TETRAHEDRON ) {
        elemTypeM_ = 2; // tetrahedron
        numNodesPerLinElementM_ = 
            corlib::ShapeTraits<corlib::TETRAHEDRON>::numVertices;
        nCommonNodes_ = 3;
    }
    else if ( eleShape == corlib::HEXAHEDRON ) {
        elemTypeM_ = 3; // hexahedron
        numNodesPerLinElementM_ = 
            corlib::ShapeTraits<corlib::HEXAHEDRON>::numVertices;
        nCommonNodes_ = 4;
    }
    else {
        std::string shapeString = corlib::convertShapeEnumToString( eleShape ); 
        FTL_VERIFY_DESCRIPTIVE( false,
                                "Unsupported type of element  %s \n", 
                                shapeString.c_str()  );
    }
}

//------------------------------------------------------------------------------
/** Given mesh is divided into required number of parts. Partitionings of 
 *  elements and nodes is stored in private data.
 */
void tools::input::partitioners::metis::
PartitionerMetis::createPartition( std::vector< std::vector<unsigned> > & connectivity,
                                   const int numParts, const unsigned overlapLayers )
{

#ifdef VERBOSE
    boost::timer routineTimer;
    routineTimer.restart();
#endif

    // make o note on the parts
    numParts_      = numParts;
    overlapLayers_ = overlapLayers;

    // check inputs
    FTL_VERIFY_DESCRIPTIVE( connectivity.size() == static_cast<unsigned>( numElements_ ), 
                            "number of element mismatch" );

    // remap nodes so that connectivity may be truncated to a linear mesh
    std::vector<METIS_Index_> lmNodesOld2new( numNodes_, -1 );

    // mark connectivity nodes in the first portion of all elements
    unsigned glbNode = 0;
    for ( int e = 0; e < numElements_; e ++ ) {
        // go throught nodes of this element
        for ( unsigned i = 0; i < static_cast<unsigned>( numNodesPerLinElementM_ ); i++ ) {
            unsigned nodeInd = connectivity[e][i];
            // mark this node 
            if ( lmNodesOld2new[ nodeInd ] == -1 ) {
                lmNodesOld2new[ nodeInd ] = glbNode;
                glbNode++;
            }
        }
    }
    // number of nodes in linear mesh
    METIS_Index_ numNodesLin = static_cast<METIS_Index_>( glbNode );
    // now finish with higher order nodes
    for ( unsigned in = 0; in < static_cast<unsigned>( numNodes_ ); in++ ) {
        if ( lmNodesOld2new[in] == -1 ) {
            // node not assigned yet, i.e. it must be higher order node
            lmNodesOld2new[ in ]      = glbNode;
            glbNode++;
        }
    }
    // check mappings
    FTL_VERIFY_DESCRIPTIVE( glbNode == static_cast<unsigned>( numNodes_ ), "error in permutation, did not achieve global number of indices" );
    FTL_VERIFY_DESCRIPTIVE( find( lmNodesOld2new.begin( ), lmNodesOld2new.end( ), -1 ) == lmNodesOld2new.end( ),
                            "error in permutation, missing indices" );

    // Prepare space for linearized connectivity table for METIS
    std::vector<METIS_Index_> nodeVec;
    nodeVec.reserve( numElements_ * numNodesPerLinElementM_ );       
    std::vector<METIS_Index_> nodeVecPoint;
    nodeVecPoint.reserve( numElements_ + 1 );       
    nodeVecPoint.push_back( 0 );
    // go through all elements and fill array for METIS
    for ( int e = 0; e < numElements_; e ++ ) {
        // go throught nodes of this element
        for ( unsigned i = 0; i < static_cast<unsigned>( numNodesPerLinElementM_ ); i++ ) {
            // get node number
            METIS_Index_ v = static_cast<METIS_Index_>( lmNodesOld2new[ connectivity[e][i] ]);
            // store node number representing an element to node edge
            nodeVec.push_back( v );
        }
        nodeVecPoint.push_back( numNodesPerLinElementM_ );
    }
    // change numbers of nodes to starts
    for ( int e = 0; e < numElements_; e ++ ) {
        nodeVecPoint[e+1] += nodeVecPoint[e];
    }

#ifdef VERBOSE
    std::cout << " Time for creating first order mesh for METIS: " << routineTimer.elapsed() << std::endl; 
    routineTimer.restart();
#endif

    // clear memory
    lmNodesOld2new.clear();
    std::vector<METIS_Index_>().swap( lmNodesOld2new );

    // debugging export of linear mesh
    //std::ofstream smf( "test.smf" );
    //// write header
    //smf << "! elementShape tetrahedron" << endl;
    //smf << "! elementNumPoints 4" << endl;
    //smf << numNodesLin << "  " << numElements_ << endl;
    //for ( int p = 0; p < numNodesLin; p++ ) {
    //    for ( int comp = 0; comp < coordinates[nodesNew2old[p]].size(); comp++ ) {
    //        smf << "  " << coordinates[nodesNew2old[p]][comp];
    //    }
    //    smf << endl;
    //}
    //for ( int e = 0; e < numElements_; e ++ ) {
    //    // go throught nodes of this element
    //    for ( unsigned i = 0; i < numNodesPerLinElementM_; i++ ) {
    //        smf << "  " << nodesNew2old[ connectivity[e][i] ];
    //    }
    //    smf << endl;
    //}

    // debug
    //std::cout << "Elements:" << std::endl; 
    //for ( int p = 0; p < nodeVec.size(); p++ ) {
    //    std::cout << p << "  " << nodeVec.at( p ) << std::endl;
    //}

    // partition mesh
    unsigned edgeCut;         // number of cut edges in final division
    if ( numParts_ > 1 ) {
        METIS_Index_ pnumflag = 0;
        METIS_Index_ numElementsM = static_cast<METIS_Index_>( numElements_ );

        // create dual graph of the mesh
#if (METIS_VER_MAJOR >= 5)
        METIS_Index_ * pxadj;
        METIS_Index_ * padjncy;
        METIS_MeshToDual( &numElementsM, &numNodesLin, &(nodeVecPoint[0]), &(nodeVec[0]), 
                          &nCommonNodes_, &pnumflag, &pxadj, &padjncy );
#else
        std::vector<METIS_Index_> xadj;
        std::vector<METIS_Index_> adjncy;

        xadj.resize( numElements_ + 1 );
        adjncy.resize( numNodesPerLinElementM_ * numElements_ ); 
        METIS_MeshToDual( &numElementsM, &numNodesLin, &(nodeVec[0]), &elemTypeM_, 
                          &pnumflag, &(xadj[0]), &(adjncy[0]) );
#endif
        nodeVecPoint.clear();
        std::vector<METIS_Index_>().swap( nodeVecPoint );
        nodeVec.clear();
        std::vector<METIS_Index_>().swap( nodeVec );

#ifdef VERBOSE
        std::cout << " Time for creating dual graph by METIS: " << routineTimer.elapsed() << std::endl; 
        routineTimer.restart();
#endif

        METIS_Index_ wgtflag=0;
        METIS_Index_ options[10];
        options[0] = 0;
        METIS_Index_ edgeCutM;
        METIS_Index_ numPartsM    = static_cast<METIS_Index_>( numParts_ );
        std::vector<METIS_Index_> elem2partM(numElements_);
        // call METIS to partition the graph
#if (METIS_VER_MAJOR >= 5)
        METIS_Index_ ncon = 1;
        METIS_PartGraphKway( &numElementsM, &ncon, pxadj, padjncy, NULL, NULL, NULL,
                             &numPartsM, NULL, NULL, NULL, &edgeCutM, &(elem2partM[0]) );
        int ierr;
        ierr = METIS_Free( pxadj );   FTL_VERIFY( ierr == METIS_OK );
        ierr = METIS_Free( padjncy ); FTL_VERIFY( ierr == METIS_OK );
#else
        METIS_PartGraphKway( &numElementsM, &(xadj[0]), &(adjncy[0]), NULL, NULL, 
                             &wgtflag, &pnumflag, &numPartsM, options, &edgeCutM, &(elem2partM[0]) );
        // erase graph
        xadj.clear( );
        std::vector<METIS_Index_>().swap( xadj );
        adjncy.clear( );
        std::vector<METIS_Index_>().swap( adjncy );
#endif
        // convert the result to unsigned integers
        std::copy( elem2partM.begin(), elem2partM.end(), elem2part_.begin() );

        // clear memory
        elem2partM.clear();
        std::vector<METIS_Index_>().swap( elem2partM );

        edgeCut = static_cast<unsigned>( edgeCutM );

    }
    else if ( numParts_ == 1 ) {
        // trivial part
        std::cout << " Warning: Dividing into a single part ... Is that what you want? " << std::endl; 

        // set all subdomains to zero
        elem2part_.assign( numElements_, 0 );
        edgeCut = 0;
    }
    else {
        FTL_VERIFY_DESCRIPTIVE( false, " Illegal number of subdomains: %d \n", numParts_ );
    }

    // verbose result
#ifdef VERBOSE
    std::cout << " Time for partitioning of graph by METIS: " << routineTimer.elapsed() << std::endl; 
    routineTimer.restart();
    std::cout << " Number of cut edges in resulting division: " << edgeCut << std::endl; 
#endif

    // process partitioning
    this->processPartitioning_( connectivity );

    return;
}

