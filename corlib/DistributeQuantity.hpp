// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//
//                 Computational Structural Mechanics Lab
//                         University of Cambridge
//
// This software is copyrighted and all rights are retained by the University of
// Cambridge. All use, disclosure, and/or reproduction of any part not expressly
// authorised by the University of Cambridge is prohibited. (C) 2010
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

//! @file DistributeQuantity.hpp

#ifndef corlib_distributequantity_h
#define corlib_distributequantity_h

#include <boost/function.hpp>

//==============================================================================
namespace corlib{

    template<typename ITER, typename ITEM>
    class DistributeQuantity;

}


//==============================================================================
/// Distribute a container holding #Quantity's to #Quantity's accessed 
/// by an function of an ITEM
///
/// Example: You want to distribute a std::vector with co-ordinates (=quantity)
///          to the nodes (=item) of a mesh.
///
/// \tparam ITEM    Item type having items accessible through functions
/// \tparam VAL     Quantity type
template<typename ITEM, typename VAL>
class corlib::DistributeQuantity
{
public:
    typedef VAL                                               Quantity;
    typedef boost::function<void( ITEM *, const Quantity & )> SetItemQuantity;
    typedef void                                              result_type;

public:
    /// Constructor
    ///
    /// \param[in]  setQuantity  Function to set quantity of an item
    DistributeQuantity( SetItemQuantity setQuantity )
        : setQuantity_( setQuantity )
    { }

    /// Copy quantity from container and set it in item
    ///
    /// Increment iterator, extract an item, use it to set item's quantity with
    /// the stored method #setQuantity_
    /// \param[in]  item   Pointer to mutable item
    /// \param[in]  iter   Forward iterator to container with quantities
    template<typename ITERATOR>
    void operator()( ITEM * item, ITERATOR & iter )
    {
        setQuantity_( item, *iter );
        ++iter;
    } 

private:
    /// Operation to indirectly set items with function
    SetItemQuantity   setQuantity_;

};

#endif
