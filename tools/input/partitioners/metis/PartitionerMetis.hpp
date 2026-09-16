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


#ifndef tools_input_partitioners_metis_partitionermetis_h
#define tools_input_partitioners_metis_partitionermetis_h

//------------------------------------------------------------------------------
#include <mpi.h>
extern "C" {
    #include <metis.h>
}
#include <vector>
#include <string>

#include <boost/timer.hpp>

#include <tools/input/partitioners/utils/PartitionerBasic.hpp>

//------------------------------------------------------------------------------
namespace tools{
    namespace input{
        namespace partitioners{
            namespace metis{

                class PartitionerMetis;

            }
        }
    }
}



//------------------------------------------------------------------------------
/** \brief Partitioning of a mesh using METIS 
*  \details This object provides a wrapper for calling METIS division into subdomains.
*           The code was developed and tested with METIS v4.0
*           The interface for calling METIS v5.0 is different. 
*/
class tools::input::partitioners::metis::PartitionerMetis
    : public tools::input::partitioners::utils::PartitionerBasic
{

public:

    //! Initialize a partition
    PartitionerMetis( const unsigned numNodes, const unsigned numElements,
                      const corlib::shape eleShape, const unsigned numNodesPerElement );

    //! Create a partitioning of the passed mesh 
    void createPartition( std::vector< std::vector<unsigned> > & connectivity,
                          const int numParts, const unsigned overlapLayers );

private: 

#if (METIS_VER_MAJOR >= 5)
    typedef idx_t   METIS_Index_;
#else
    typedef idxtype METIS_Index_;
#endif

    // basic partitioner
    typedef tools::input::partitioners::utils::PartitionerBasic PartitionerBasic_;

    // mesh constants
    METIS_Index_ numNodesPerLinElementM_;              //!< number of nodes in corresponding LINEAR elements - constant
    METIS_Index_ elemTypeM_;                           //!< type of elememt as for METIS 
                                                      //!< ( 1 - triangle, 2 - tetra, 3 - hexa, 4 - quad )
    METIS_Index_ nCommonNodes_;                        //!< number of nodes two elements need to share to define an edge in a dual mesh graph
};

#include "PartitionerMetis.ipp"

#endif
