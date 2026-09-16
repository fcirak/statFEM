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

#ifndef tools_input_rcb_partitionerrcb_h
#define tools_input_rcb_partitionerrcb_h

#include <vector>
#include <string>

#include <tools/input/partitioners/utils/misc.hpp>
#include <tools/input/partitioners/utils/PartitionerBasic.hpp>
#include <tools/input/partitioners/rcb/algo/BasicRCB.hpp>

//------------------------------------------------------------------------------
namespace tools{
    namespace input{
        namespace partitioners{

            /// Recursive Coordinate Bisection namespace
            namespace rcb{
                namespace mesh{

                    namespace ublas = boost::numeric::ublas;
                    
                    class PartitionerRCB;
                    
                }
            }
        }
    }
}



//------------------------------------------------------------------------------
/** \brief Partitioning of a mesh using Recursive Coordinate Bisection (RCB) 
 *  \details This object provides a wrapper for calling RCB division into 
 *  subdomains.
 */
class tools::input::partitioners::rcb::mesh::PartitionerRCB
    : public tools::input::partitioners::utils::PartitionerBasic
{
//------------------------------------------------------------------------------
// Public interface
public:

    //! Initialize a partition
    PartitionerRCB( const unsigned numNodes, const unsigned numElements,
                    const unsigned numNodesPerElement );

    //! Create a partitioning of the passed mesh 
    void createPartition( std::vector< std::vector<unsigned> > & connectivity,
                          const std::vector< ublas::bounded_vector<double,3> > & coords,
                          const int numParts, const int overlapLayers );


// Private data members
private:
    //
    typedef tools::input::partitioners::utils::PartitionerBasic PartitionerBasic_;

};

//------------------------------------------------------------------------------
#include "PartitionerRCB.ipp"

#endif
