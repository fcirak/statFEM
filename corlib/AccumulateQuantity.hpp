// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//
//                 Computational Structural Mechanics Lab
//                         University of Cambridge
//
// This software is copyrighted and all rights are retained by the University of
// Cambridge. All use, disclosure, and/or reproduction of any part not expressly
// authorised by the University of Cambridge is prohibited. (C) 2010
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

//! @file AccumulateQuantity.hpp

#ifndef corlib_accumulatequantity_h
#define corlib_accumulatequantity_h

#include <functional>

#include <boost/function.hpp>


//==============================================================================
// declarations
namespace corlib {
    template< typename ITEM, typename RVAL, template <typename> class OP >
    class AccumulateQuantity;    
}



//==============================================================================
// implementations

//------------------------------------------------------------------------------
/** \brief Accumulate data across a container of ITEM's 
 *  \details The data is determined by the access operator #getFun_
 *
 *  \tparam ITEM    item whose member data is accumulated (i.e. Element or Node )
 *  \tparam RVAL    value type of data
 *  \tparam OP      binary operation (i.e. plus(default), multiplies ) 
 */
template< typename ITEM, typename RVAL, template <typename> class OP = std::plus >
class corlib::AccumulateQuantity 
{
private:
    typedef boost::function< RVAL (ITEM *) > GetFun_;

public:
    //! Constructor 
    AccumulateQuantity( const GetFun_ getFun, RVAL init = RVAL() ) :
        getFun_( getFun ), value_( init )
        { }

    //! function operator which adds ITEM's quantity
    //! to global #value_
    void operator()( ITEM * ip )
    {        
        value_ = OP<RVAL> ()(value_, getFun_( ip ));
        
        return;
    }

    //! Return accumulated quantity  (with implicit conversion operator)
    //!
    //! \return quantity accumulated across ITEMs
    operator RVAL ( ) const { return value_; }

private:
    //! operator to determine ITEM's contribution
    const GetFun_  getFun_;
    //! sum of all ITEM's contribution (external variable)
    RVAL           value_;
};


//------------------------------------------------------------------------------
#endif
