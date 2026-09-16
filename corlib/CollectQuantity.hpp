// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//
//                 Computational Structural Mechanics Lab
//                         University of Cambridge
//
// This software is copyrighted and all rights are retained by the University of
// Cambridge. All use, disclosure, and/or reproduction of any part not expressly
// authorised by the University of Cambridge is prohibited. (C) 2010
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

//! @file CollectQuantity.hpp

#ifndef corlib_collectquantity_h
#define corlib_collectquantity_h

//! boost includes
#include <boost/function.hpp>
#include <boost/iterator/iterator_concepts.hpp>

//==============================================================================
// declarations
namespace corlib {
    template< typename ITEM, typename RVAL >
    class CollectQuantity;    
}



//==============================================================================
// implementations

//------------------------------------------------------------------------------
/** \brief Collect data across a container of ITEM's 
 *  \details Given an accessor of type RVAL(ITEM *), this object retrieves the
 *  datum of type RVAL from an item of type ITEM and 'dereference assigns' it
 *  to a given iterator which has to be writable and incrementable. Typically, 
 *  the iterator is either of type vector<RVAL>::iterator or 
 *  back_insert_iterator< vector<RVAL> >. 
 *
 *  \tparam ITEM    item whose member data is collected (i.e. Element or Node )
 *  \tparam RVAL    value type of data
 */
template<typename ITEM, typename RVAL>
class corlib::CollectQuantity 
{
private:
    typedef boost::function< RVAL (ITEM *) > GetFun_;

public:
    typedef void result_type;

    //! Constructor 
    CollectQuantity( const GetFun_ getFun )
        : getFun_( getFun )
    { }

    //--------------------------------------------------------------------------
    /** Function operator which passes ITEM's quantity to given iterator
     *  \tparam ITER Type of iterator (e.g, back-inserter or normal iterator)
     *  \param[in] ip  Pointer to an object which provides the datum
     *  \param[in] out Iterator which is dereferenced, assigned to and incremented
     */
    template<typename ITER>
    void operator()( ITEM * ip, ITER & out )
    {   
        // Following operation requires the concept of 'Dereference Assignment'
        // Note back_insert_iterator has value_type=void, hence the second
        // argument in the concept check is needed
        BOOST_CONCEPT_ASSERT((boost_concepts::WritableIterator<ITER,RVAL> )); 
        *out = getFun_( ip );

        // Postincrement concept is required
        BOOST_CONCEPT_ASSERT((boost_concepts::IncrementableIterator<ITER> )); 
        ++out;
        
        return;
    }

private:
    //! operator to determine ITEM's contribution of type RVAL
    const GetFun_  getFun_;
};


//------------------------------------------------------------------------------
#endif
