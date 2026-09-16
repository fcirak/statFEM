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

#ifndef tools_input_partitioners_utils_partitionerbasic_h
#define tools_input_partitioners_utils_partitionerbasic_h

//------------------------------------------------------------------------------
#include <vector>
#include <string>
#include <set>
#include <map>
#include <cassert>

#include <boost/timer/timer.hpp>
#include <boost/type_traits.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/bind.hpp>

#include <corlib/verify.hpp>
#include <corlib/misc.hpp>

//------------------------------------------------------------------------------
namespace tools{
    namespace input{
        namespace partitioners{
            namespace utils{

                class PartitionerBasic;

                namespace ublas = boost::numeric::ublas;

            }
        }
    }
}



//------------------------------------------------------------------------------
/** \brief Partitioning of a mesh (abstract base class for PartitionerMetis, PartitionerRcb)
 *
 *  \details This class provides common methods to act on an already partitioned mesh.
 *           The partitioner (METIS, RCB) provides a non-overlapping division of 
 *           elements in the elem2part_ vector. 
 *           Subsequently, in this class, nodes are renumbered 
 *           to achieve optimal overlap of local nodes with local elements.
 */
class tools::input::partitioners::utils::PartitionerBasic
{
public:

    //! Initialize a partition
    PartitionerBasic( const unsigned numNodes, const unsigned numElements,
                      const unsigned numNodesPerElement );

    //! Retrieve status
    bool createdPartitions() const { return isPartitioned_; }

protected:
    //!
    void processPartitioning_( std::vector< std::vector<unsigned> > & connectivity );

private:
    //! Update the connectivity matrix to new node numbering
    void updateConnectivityNumbering_( std::vector< std::vector<unsigned> > & connectivity,
                                       bool backward = false );

    //! Transpose the connectivity matrix, i.e. generate list of elements at nodes.
    void transposeConnectivity_( const std::vector< std::vector<unsigned> > & connectivity );

public:
    //! Update the numbering of given vector of global nodes
    void updateListNumbering( std::vector<unsigned> & listOld,
                              std::vector<unsigned> & listNew,
                              bool backward ) const;

    //! Create submesh using privately stored partitioning
    void getSubMesh( const std::vector< ublas::bounded_vector<double,3> > & coordinates, 
                     const std::vector< std::vector<unsigned> > & elements,
                     const unsigned subIndex,
                     std::vector< std::vector<unsigned> > & elementsSub,
                     std::vector<unsigned> & elementsGlobal4local,
                     std::vector<unsigned> & elementsOwner,
                     std::vector< ublas::bounded_vector<double,3> > & coordinatesSub, 
                     std::vector<unsigned> & nodesGlobal4local,
                     std::vector<unsigned> & nodesOwner );

    

    // debugging
    void printPartitions( ) const;


protected:
    //! return part number in linear partition given by starts for an entry
    static
    unsigned getSubIndex_( const unsigned index,
                           const std::vector<unsigned> & starts );

protected: 
    // mesh constants
    unsigned numNodes_;                       //!< total number of nodes
    unsigned numElements_;                    //!< total number of elements
    unsigned numNodesPerElement_;             //!< number of nodes in elements - constant

    // partitions
    unsigned numParts_;                       //!< target number of parts in generated division
    std::vector<int> elem2part_;              //!< subdomain number for each element

    unsigned overlapLayers_;                  //!< number of layers around each non-overlapping subdomain
                                              //!< this is used for faster PETSc computations

private:
    std::vector<int> node2part_;         //!< subdomain number for each node

    std::vector< std::set<unsigned> > elems4nodes_; //!< transposed connectivity
                                                    //!<  - list of indices of elements for each node

    // mapping of nodes
    std::vector<int> nodesOld2new_;      //!< permutation of nodes of original mesh for new optimal numbering - old number 
                                              //!< is at new position
    std::vector<int> nodesNew2old_;      //!< permutation of nodes of original mesh for new optimal numbering - new number 
                                              //!< is at old position
    std::vector<unsigned> sub2nodeStart_;     //!< where subdomain nodes start if subdomains 0, 1, 2 has eg 4, 3, 6 nodes,
                                              //!< this array looks as [ 0 4 7 13 ] (accumulated nodes per subdomain)
    // mapping of elements
    std::vector<int> elementsOld2new_;   //!< permutation of elements of original mesh for new optimal numbering - old number 
                                              //!< is at new position
    std::vector<int> elementsNew2old_;   //!< permutation of elements of original mesh for new optimal numbering - new number 
                                              //!< is at old position
    std::vector<unsigned> sub2elementStart_;  //!< where subdomain elements of nonoverlapping division start 

    bool buildCluster_;                       //!< flag if an overlapping partition is to be created

    // status
    bool isPartitioned_;                      //!< flag if class contains already partition and renumbered mappings
};

//------------------------------------------------------------------------------
#include "PartitionerBasic.ipp"

#endif
