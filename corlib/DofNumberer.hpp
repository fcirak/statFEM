// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//
//                 Computational Structural Mechanics Lab
//                         University of Cambridge
//
// This software is copyrighted and all rights are retained by the University of
// Cambridge. All use, disclosure, and/or reproduction of any part not expressly
// authorised by the University of Cambridge is prohibited. (C) 2010
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

//! @file DofNumberer.hpp

#ifndef corlib_dofnumberer_h
#define corlib_dofnumberer_h
//------------------------------------------------------------------------------
#include <functional>
#include <utility>
#include <vector>

#include <boost/function.hpp>
#include <boost/bind.hpp>

//------------------------------------------------------------------------------
namespace corlib{
    template<typename DOFNUM> class DofNumberer;

    //! Convenience function given a function pointer
    template<typename NODE>
    corlib::DofNumberer<boost::function< void(NODE *, unsigned &) > >
    dofNumberFun( void( NODE::*f)( unsigned & ), const unsigned init = 0)
    {
        typedef boost::function< void(NODE*, unsigned &) > DofNumberFun;
        return 
            corlib::DofNumberer<DofNumberFun> ( boost::bind( f, _1, _2 ), init );
    }

}

//------------------------------------------------------------------------------
/** \brief Functor to consecutively number degrees of freedom
 *  \details Given a functor providing the node's (or whatever) dof-numbering
 *  member function, this object passes a reference to a counter to it.
 *  It is assumed that this functor assigns the current counter value to the 
 *  degrees of freedom and increments its value.
 *  \tparam DOFNUM Type of function object around 'void(NODE::*)( unsigned &)'
 */
template<typename DOFNUM>
class corlib::DofNumberer
{
public:
    typedef typename DOFNUM::first_argument_type NodePtr;

    //! Cstor given the member function functor and an initial value
    DofNumberer( DOFNUM dofNum, const unsigned init = 0 )
        : dofNum_( dofNum ), counter_( init )
    { }

    //! Apply functor to each node (or whatever)
    void operator()( NodePtr np  )
    {
    	dofNum_( np, counter_ );
    	return;
    }

    //! Return the valuer of the counter
    unsigned getTotalNumber() const { return counter_; }

    //! Nicer: use implicit conversion operator
    operator unsigned() { return counter_; }

private:
    DOFNUM     dofNum_;  //!< Functor for dof-numbering
    unsigned   counter_; //!< Counter
};

#endif
