// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//
//                 Computational Structural Mechanics Lab
//                         University of Cambridge
//
// This software is copyrighted and all rights are retained by the University of
// Cambridge. All use, disclosure, and/or reproduction of any part not expressly
// authorised by the University of Cambridge is prohibited. (C) 2010
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

//! @file DofNumbererFromStarts.hpp

#ifndef corlib_dofnumbererfromstarts_h
#define corlib_dofnumbererfromstarts_h
//------------------------------------------------------------------------------
#include <functional>
#include <utility>
#include <vector>
#include <map>

#include <boost/function.hpp>
#include <boost/bind.hpp>

//------------------------------------------------------------------------------
namespace corlib{
    template<typename DOFNUM> class DofNumbererFromStarts;

    //! Convenience function given a function pointer
    template<typename NODE>
    corlib::DofNumbererFromStarts< boost::function< void(NODE *, unsigned &) > >
    dofNumberFromStartsFun( void( NODE::*f)( unsigned & ), const std::map<int,int> * starts )
    {
        typedef boost::function< void(NODE*, unsigned &) > DofNumberFun;
        return 
            corlib::DofNumbererFromStarts<DofNumberFun> ( boost::bind( f, _1, _2 ), starts );
    }

}

//------------------------------------------------------------------------------
/** \brief Functor to number degrees of freedom based on starts.
 *  \details Takes a functor providing the node's (or whatever) dof-numbering
 *  member function. In addition, a map with starting dof for global IDs, (ID,startingDof),
 *  is given.
 *  This object asks the node for ID, gets from the map from which dof the node should assign its dofnumbers,
 *  and with this start launches the dof numbering method of the node.
 *  It is assumed that the node functor assigns consequtive numbers to the 
 *  degrees of freedom beginning with the start value.
 *  \tparam DOFNUM Type of function object around 'void(NODE::*)( unsigned &)'
 */
template<typename DOFNUM >
class corlib::DofNumbererFromStarts
{
public:
    typedef typename DOFNUM::first_argument_type NodePtr;

    //! Cstor given the member function functor and array of starts
    DofNumbererFromStarts( DOFNUM dofNum, const std::map<int,int> * starts )
        : dofNum_( dofNum ), starts_( starts )
    { }

    //! Apply functor to each node (or whatever)
    void operator()( NodePtr np  )
    {
        // set start of degrees of freedom for that node according to its global ID
        unsigned iD =  np -> giveId( );

        // look for ID in the map of IDs to global dof starts
        std::map<int,int>::const_iterator pos = starts_ -> find( iD );
        FTL_VERIFY_DESCRIPTIVE( pos != starts_ -> end(),
                                "Cannot remap ID %d to global dof start. \n ", iD );
        unsigned start = pos -> second;

        // let the node number its dofs beginning from that start
        dofNum_( np, start );
        return;
    }

private:
    DOFNUM dofNum_;                    //!< Functor for dof-numbering
    const std::map<int,int> * starts_; //!< Pointer to array with start of degrees of freedom for each global node
};

#endif
