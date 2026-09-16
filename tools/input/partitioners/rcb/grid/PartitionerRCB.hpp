//
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//
//                 Computational Structural Mechanics Lab
//                         University of Cambridge
//
// This software is copyrighted and all rights are retained by the CSMLab. All
// use, disclosure, and/or reproduction of any part not expressly authorized by
// F Cirak is prohibited. (C) 2012.
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//

//! @file   rcb/grid/PartitionerRCB.hpp
//! @author Jakub Sistek, Matija Kecman, Fehmi Cirak
//! @date   3/2012

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
                namespace grid{

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
class tools::input::partitioners::rcb::grid::PartitionerRCB
    : public tools::input::partitioners::utils::PartitionerBasic
{
//------------------------------------------------------------------------------
// Public interface
public:

    //! Initialize a partition
    PartitionerRCB( const unsigned numNodes, const unsigned numElements,
                    const unsigned numNodesPerElement );

    //! Create a partitioning of the passed mesh 
    void createPartition( const std::vector< std::vector<unsigned> > & connectivity,
                          const std::vector< ublas::bounded_vector<double,3> > & coords,
                          const int numParts );

    //! Returns the elem2part_ vector
    std::vector<int> giveElement2PartMapping() const { return elem2part_; };

// Private data members
private:
    //
    typedef tools::input::partitioners::utils::PartitionerBasic PartitionerBasic_;

};

//------------------------------------------------------------------------------
#include "PartitionerRCB.ipp"

#endif
