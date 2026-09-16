// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//
//                 Computational Structural Mechanics Lab
//                         University of Cambridge
//
// This software is copyrighted and all rights are retained by the University of
// Cambridge. All use, disclosure, and/or reproduction of any part not expressly
// authorised by the University of Cambridge is prohibited. (C) 2010
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

//! @author Jakub Sistek

#ifndef corlib_generatedofstarts_h
#define corlib_generatedofstarts_h

//------------------------------------------------------------------------------

#include <mpi.h>
#include <vector>
#include <ostream>

#include <corlib/Mesh.hpp>

#include <boost/bind.hpp>
#include <boost/function.hpp>

//------------------------------------------------------------------------------
namespace corlib{

    namespace detail_ {
/** \brief Helper function changing counts to starts.
 */
        template <typename TYPE> 
        void counts2Starts ( std::vector<TYPE> & array ) {
            unsigned i;
            unsigned n = array.size()-1;
            for (i=1; i<n; i++) array[i] += array[i-1]; 
            for (i=n; i>0; i--) array[i] = array[i-1]; 
            array[0] = 0; 
        }

/** \brief Helper function for searching rank of processor in map of starts 
 *         if we have array of starts:   0 5 9 14 16
 *         of processors:                0 1 2 3  
 *         One knows, that rank 0 has nodes 0-4
 *                         rank 1 has nodes 5-8
 *                         rank 2 has nodes 9-13
 *                     and rank 3 has nodes 14-15
 *         The initial simple array was converted to map of the form (start,rank) to make search faster.
 */
        template <typename TYPE> 
        int proc4Id ( TYPE iD, std::map<TYPE,TYPE> & procStarts ) {
            typename std::map<TYPE,TYPE>::iterator posUp = ( procStarts.lower_bound( iD ) );
            typename std::map<TYPE,TYPE>::iterator posLow;
            FTL_VERIFY_DESCRIPTIVE( posUp != procStarts.end(), "I am behind the array for ID %d \n", iD );
            if ( posUp -> first != iD ) {
                // if keys do not match, the lower_bound returns first larger value
                FTL_VERIFY_DESCRIPTIVE( posUp != procStarts.begin(), "I am already in the beginning for node ID %d \n", iD );
                posLow = --posUp;
            }
            else {
                posLow = posUp;
            }

            int proc = posLow -> second;

            //std::cout << "Claiming that node " << iD << " belongs to proc " << proc << std::endl;
            return proc;
        }
    }

/** \brief Function for generating dof numbering for parallel applications.
 *  \details The function takes local meshes at each processor. It assumes that 
 *  nodes are able to return number of dofs they contain.
 *  It returns a local map with starting (global) dof number for each node. 
 *  Total number of nodes and degrees of freedom is returned in pair.
 *  \param[in] subMesh  subdomain mesh 
 *  \param[in] comm     MPI communicator
 *  \param[out] dofStartsMap  starts of degrees of freedom at nodes
 *  \param[in] queryNumDofs  functor to obtain number of degrees of freedom at node
 *  \param[in] queryActivity functor asking for activity of nodes
 */
    template< typename MESH >
    std::pair<int,int>
    generateDofStarts( MESH & subMesh, MPI_Comm comm, std::map<int,int> & dofStartsMap, 
                            typename boost::function< unsigned ( const typename MESH::Node * ) > queryNumDofs 
                                = boost::bind( &MESH::Node::giveNumDofs, _1 ),
                            typename boost::function< bool ( const typename MESH::Node * ) > queryActivity
                                = boost::bind( &MESH::Node::isDDActive, _1 ) ) {

        typedef std::map<int,int> DofStartsMap;

        // Orient in the communicator
        int nProc, rank;
        MPI_Comm_size( comm, &nProc );
        MPI_Comm_rank( comm, &rank  );

        // Collect numbers of local dof starts at active nodes
        std::vector<int> localDofsStarts;
        subMesh.iterateOverNodesWithPredicate( boost::bind( &std::vector<int>::push_back, &localDofsStarts, 
                                                            boost::bind( queryNumDofs, _1) ),
                                               boost::bind( queryActivity, _1 ) );
        // reserve space for one more entry for conversion to starts
        localDofsStarts.push_back(0);
        // change the array to starts
        detail_::counts2Starts ( localDofsStarts );

        // gather lists of global Ids of local and non-local nodes in each subdomain
        boost::function< unsigned ( const typename MESH::Node * ) > queryIds = boost::bind( &MESH::Node::giveId, _1 );
        std::vector<int> globalIdsOfLocalNodes;
        subMesh.iterateOverNodesWithPredicate( boost::bind( &std::vector<int>::push_back, &globalIdsOfLocalNodes,
                                                            boost::bind( queryIds, _1 ) ),
                                               boost::bind( queryActivity, _1 ) );
        std::vector<int> globalIdsOfNonLocalNodes;
        subMesh.iterateOverNodesWithPredicate( boost::bind( &std::vector<int>::push_back, &globalIdsOfNonLocalNodes, 
                                                            boost::bind( queryIds, _1 ) ),
                                               not boost::bind( queryActivity, _1 ) );

        // determine number of local unique nodes
        int numLocalNodes    = globalIdsOfLocalNodes.size();
        int numNonLocalNodes = globalIdsOfNonLocalNodes.size();
        int numLocalDofs     = localDofsStarts[ numLocalNodes ];

        int numTotalNodes; 
        MPI_Allreduce( &numLocalNodes, &numTotalNodes, 1, MPI_INT, MPI_SUM, comm );
        int numTotalDofs;
        MPI_Allreduce( &numLocalDofs,  &numTotalDofs,  1, MPI_INT, MPI_SUM, comm );

        // determine number of total unique nodes on all processors
        std::vector<int> procNodesStarts( nProc + 1, 0 );
        MPI_Allgather( &numLocalNodes,        1, MPI_INT,
                       &(procNodesStarts[0]), 1, MPI_INT, comm );
        // change the array to starts
        detail_::counts2Starts ( procNodesStarts );
        //std::cout << "procNodesStart" << std::endl;
        //std::copy( procNodesStarts.begin(), procNodesStarts.end(), std::ostream_iterator<int>(std::cout, " ") );

        // convert the vector to map
        std::map<int,int> procNodesStartsMap;
        for ( int i = 0; i < procNodesStarts.size( ); i++ ) {
            procNodesStartsMap.insert( std::make_pair( procNodesStarts[i], i ) );
        }


        // determine number of total unique dofs
        std::vector<int> procDofsStarts(  nProc + 1, 0 );
        MPI_Allgather( &numLocalDofs,        1, MPI_INT,
                       &(procDofsStarts[0]), 1, MPI_INT, comm );
        // change the array to starts
        detail_::counts2Starts ( procDofsStarts );

        // Create a map of local nodes [globalID, dofStart]
        for ( int i = 0; i < globalIdsOfLocalNodes.size(); i++ ) {
            int shiftedStart = localDofsStarts[i] + procDofsStarts[rank];
            dofStartsMap.insert( std::make_pair( globalIdsOfLocalNodes[i], shiftedStart ) );
        }

        // Now let us handle the nonlocal nodes...
        
        //  find owners of non-local nodes
        //                neighbour       IDs               dof starts
        typedef std::map< int, std::pair< std::vector<int>, std::vector<int> > > Neighbours;
        Neighbours neighbours;
        std::vector<int>::iterator it1  = globalIdsOfNonLocalNodes.begin();
        std::vector<int>::iterator it1e = globalIdsOfNonLocalNodes.end();
        for ( ; it1 != it1e; it1++ ) {
            int iD = *it1;
            int owner = detail_::proc4Id( iD, procNodesStartsMap );
            Neighbours::iterator pos = neighbours.find( owner );
            if ( pos != neighbours.end( ) ) {
                // the rank already exists in neighbours, just add a new entry into its list of nodes
                (( pos -> second ).first).push_back( iD );
            }
            else {
                // introduce a new owner with a vector
                std::vector<int> list( 1, iD );
                std::vector<int> dofStarts( 1, 0 );
                neighbours.insert( std::make_pair( owner, std::make_pair( list, dofStarts ) ) );
            }
        }

        // find number of my neighbours
        int numNeighbours = neighbours.size( );
        // and maximum of it
        int numNeighboursMax;
        MPI_Allreduce( &numNeighbours, &numNeighboursMax, 1, MPI_INT, MPI_MAX, comm );

        // create the 2dimensional array for communication
        std::vector<int> whoWantsWhoTmp( numNeighboursMax, 0 );
        typename Neighbours::iterator it2   = neighbours.begin();
        typename Neighbours::iterator it2E  = neighbours.end();
        int index = 0;
        for ( ; it2 != it2E; ++it2 ) {
            int owner = it2 -> first;
            whoWantsWhoTmp[index] = owner;
            index++;
        }

        // find global connectivity pattern
        int lWhoWantsWho = numNeighboursMax * nProc;
        std::vector<int> whoWantsWho( lWhoWantsWho );
        MPI_Allgather( &(whoWantsWhoTmp[0]), numNeighboursMax, MPI_INT,
                       &(whoWantsWho[0]),    numNeighboursMax, MPI_INT, comm );
        std::vector<int> numNeighboursAll( nProc );
        MPI_Allgather( &numNeighbours,         1, MPI_INT,
                       &(numNeighboursAll[0]), 1, MPI_INT, comm );

        // who wants to speaks to me
        std::vector<int> whoWantsMyData;
        for ( int iProc = 0; iProc < nProc; iProc++ ) {
            for ( int j = 0; j < numNeighboursAll[iProc]; j++ ) {
                int wanted = whoWantsWho[iProc*numNeighboursMax + j];

                if ( wanted == rank ) {
                    // register the rank

                    int index;
                    whoWantsMyData.push_back( iProc );
                }
            }
        }

        // Send requests to data owners 
        std::vector<MPI_Request> requests;
        typename Neighbours::iterator it4   = neighbours.begin();
        typename Neighbours::iterator it4E  = neighbours.end();
        std::vector<int> lengths( neighbours.size() );
        int ind = 0;
        for ( ; it4 != it4E; ++it4 ) {
            int owner = it4 -> first;
            int length = ((it4 -> second).first).size();
            lengths[ind] = length;

            // how many is coming
            MPI_Request oneRequest1;
            MPI_Isend( &(lengths[ind]), 1, MPI_INT, owner, rank, comm, &oneRequest1 );
            requests.push_back( oneRequest1 );

            // send IDs
            MPI_Request oneRequest2;
            MPI_Isend( &(((it4 -> second).first)[0]), length, MPI_INT, owner, rank, comm, &oneRequest2 );
            requests.push_back( oneRequest2 );

            // receive dof starts
            MPI_Request oneRequest3;
            ((it4 -> second).second).resize(length);
            MPI_Irecv( &(((it4 -> second).second)[0]), length, MPI_INT, owner, rank, comm, &oneRequest3 );
            requests.push_back( oneRequest3 );

            ind++;
        }

        // satisfy others' needs
        typename std::vector<int>::iterator it5   = whoWantsMyData.begin();
        typename std::vector<int>::iterator it5E  = whoWantsMyData.end();
        for ( ; it5 != it5E; it5++ ) {
            int sender = *it5;

            // receive length of demanded data
            int length;
            MPI_Status oneStatus1;
            MPI_Recv( &length, 1, MPI_INT, sender, sender, comm, &oneStatus1 );

            // receive array of demanded indices
            std::vector<int> desiredIds( length );
            MPI_Status oneStatus2;
            MPI_Recv( &(desiredIds[0]), length, MPI_INT, sender, sender, comm, &oneStatus2 );

            // go through the array of IDs and fill the global dof number
            std::vector<int> desiredDofStarts( length );
            for ( int i = 0; i < desiredIds.size(); i++ ) {
                int iD = desiredIds[i];

                // it should be in my map
                typename DofStartsMap::const_iterator pos = dofStartsMap.find( iD );
                if ( pos == dofStartsMap.end() ) {
                    std::cout << " DesiredIds "; 

                    std::copy( desiredIds.begin(), desiredIds.end(), std::ostream_iterator<int>( std::cout, " " ) );
                    std::cout << std::endl;
                }
                FTL_VERIFY_DESCRIPTIVE( pos != dofStartsMap.end(),
                                        "Something went wrong! Cannot find node ID %d although I should have it. \n ", iD );
                int globalDofStart = pos -> second;

                desiredDofStarts[i] = globalDofStart;
            }

            MPI_Send( &(desiredDofStarts[0]), length, MPI_INT, sender, sender, comm );
        }

        // wait until communication finishes
        int nRequests = requests.size();
        std::vector< MPI_Status > statuses( nRequests );
        MPI_Waitall( nRequests, &(requests[0]), &(statuses[0]) );

        MPI_Barrier( comm );
        
        // go through my portion of non-local nodes and add their dof starts to my map
        typename Neighbours::iterator it6   = neighbours.begin();
        typename Neighbours::iterator it6E  = neighbours.end();
        for ( ; it6 != it6E; ++it6 ) {
            int owner = it6 -> first;
            int length = ((it6 -> second).first).size();
            for( int i = 0; i < length; i++ ) {
                int iD           = ((it6 -> second).first)[i];
                int shiftedStart = ((it6 -> second).second)[i];
                dofStartsMap.insert( std::make_pair( iD, shiftedStart ) );
            }
        }


        // check the map
        //typename DofStartsMap::iterator it7  = dofStartsMap.begin();
        //typename DofStartsMap::iterator it7E = dofStartsMap.end();
        //for ( int irank = 0; irank < nProc; irank++ ) {
        //    if (rank == irank) {
        //        std::cout << "resulting map on " << rank << std::endl;
        //        for ( ; it7 != it7E; ++it7 ) {
        //            std::cout << it7 -> first << "  " << it7 -> second << std::endl;
        //        }
        //    }
        //    MPI_Barrier( comm );
        //}

        return std::make_pair(numTotalNodes, numTotalDofs);

    }

}

#endif
