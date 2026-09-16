// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//
//                 Computational Structural Mechanics Lab
//                         University of Cambridge
//
// This software is copyrighted and all rights are retained by the University of
// Cambridge. All use, disclosure, and/or reproduction of any part not expressly
// authorised by the University of Cambridge is prohibited. (C) 2010
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

//! @file DistributedNode.hpp

#ifndef corlib_distributednode_h
#define corlib_distributednode_h

#include <iostream>
#include <limits>

//------------------------------------------------------------------------------

namespace corlib{
    template<typename NODE> class DistributedNode;
}

//------------------------------------------------------------------------------
/** \brief  Node wrapper for distributed computation 
 *  \details This object equips the NODE it wraps around with a global 
 *  identifier and an activity flag.
 *  \tparam NODE Type of node to wrap around
 */
template<typename NODE>
class corlib::DistributedNode : public NODE
{
public:
    //! Empty constructor 
    DistributedNode( )
        : ddActive_( false )
        { }

    //! Read from stream
    std::istream & readSelf( std::istream & inp ) {
        //! - let base class read its coordinates
        this -> NODE::readSelf( inp );
        //! - read global ID and activity
        unsigned globalID;
        inp >> globalID >> ddActive_;
        //! - set ID to global ID
        this -> NODE::setId( globalID );
        return inp;
    }

    //! return activity
    bool isDDActive() const { return ddActive_; }

    //! set active
    void setDDActive() { ddActive_ = true; }

private:
    bool     ddActive_;  //!< True if active (on current process)
};

#endif
