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

//! @file   BasicRCB.hpp
//! @author Matija Kecman, Fehmi Cirak
//! @date   2011

#ifndef tools_input_partitioners_rcb_algo_basicrcb_h
#define tools_input_partitioners_rcb_algo_basicrcb_h

#include <fstream>
#include <vector>
#include <limits>
#include <bitset>
#include <set>
#include <cassert>
#include <algorithm>

#include <corlib/fuzzyEqual.hpp>

//------------------------------------------------------------------------------
namespace tools{
    namespace input{
        namespace partitioners{
            namespace rcb{
                namespace algo{

                    //----------------------------------------------------------
                    // classes/functors
                    class RCBPoint;
            
                    class RCBNode;
            
                    class RCBPointComparison;

                    //----------------------------------------------------------
                    // for convenience
                    namespace ublas = boost::numeric::ublas;
                    typedef std::vector<RCBPoint*>            RCBPointContainer;
                    typedef std::vector<RCBNode*>             RCBNodeContainer;
                    typedef typename RCBNodeContainer::const_iterator RCBNodeConstIterator;
                    typedef typename RCBPointContainer::const_iterator RCBPointConstIterator;

                    //----------------------------------------------------------
                    // This function performs the RCB algorithm
                    void performRCB( const std::vector< ublas::bounded_vector<double,3> > &,
                                     const int &,
                                     std::vector<int> & );
            
                    //----------------------------------------------------------
                    // functions
                    void findLeaves( RCBNode *, 
                                     RCBNodeContainer & );
                    
                    bool leaf( const RCBNode * );
                    
                    void getMaxBoxAndMaxDim( const RCBNode *,
                                             double &,
                                             unsigned & );
                    
                    void partitionNode( RCBNode *,
                                        double &,
                                        unsigned & );
                    
                    void extractPartitioning( const RCBNodeContainer &,
                                              std::vector<int> & );

                }
            }
        }
    }
}

//------------------------------------------------------------------------------
/** \brief An element `wrapper' for use with the RCB partitioner
 */
class tools::input::partitioners::rcb::algo::RCBPoint
{
private:
    typedef ublas::bounded_vector<double, 3> VecDim;

public:
    //! Constructor
    RCBPoint( const VecDim centerIn,
              const int pointIdIn )
        : center( centerIn ), 
          pointId( pointIdIn ) 
    { }

public:
    VecDim center;  //!< coordinates of element center
    int    pointId; //!< global ID of element
};

//------------------------------------------------------------------------------
/** \brief A node of the RCB tree. Make data a pointer?
 */
class tools::input::partitioners::rcb::algo::RCBNode
{
private:
    typedef ublas::bounded_vector<double, 3> VecDim;

public:
    //! Constructor
    RCBNode( const RCBPointContainer dataIn,
             const int levelIn )
        : data( dataIn ), level( levelIn )
    {
        childL = NULL;
        childR = NULL;
    }

    //! Destructor
    ~RCBNode(){
        if (childL) delete childL;
        if (childR) delete childR;
    }

public:
    //! Data stored at the node
    const RCBPointContainer data;
    const int level;

    //! Pointers to child nodes
    RCBNode* childL;
    RCBNode* childR;
};

//------------------------------------------------------------------------------
/** \brief RCBPoint Comparison functor.
 */
class tools::input::partitioners::rcb::algo::RCBPointComparison 
    : public std::binary_function<RCBPoint*, RCBPoint*, bool> 
{
public:
    RCBPointComparison( const int& dir ): dir_( dir ){ }
    
    bool operator()( RCBPoint* first, RCBPoint* second ) const 
    {
        return ( first->center[dir_] > second->center[dir_] );
    }
	
private:
    const int dir_; //!< Coordinate direction to be tested

};

//------------------------------------------------------------------------------
#include "BasicRCB.ipp"

#endif
