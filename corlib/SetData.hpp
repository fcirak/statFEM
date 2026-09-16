// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//
//                 Computational Structural Mechanics Lab
//                         University of Cambridge
//
// This software is copyrighted and all rights are retained by the University of
// Cambridge. All use, disclosure, and/or reproduction of any part not expressly
// authorised by the University of Cambridge is prohibited. (C) 2010
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

//! @file SetData.hpp

#ifndef corlib_setdata_h
#define corlib_setdata_h
//------------------------------------------------------------------------------
#include <corlib/verify.hpp>

#include <functional>
#include <map>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>

//------------------------------------------------------------------------------
namespace corlib {
    template<typename GET, typename SET> class SetData;
    template<typename SET> class SetDataFromFile;
    namespace ublas = boost::numeric::ublas;
}

//------------------------------------------------------------------------------
/** \brief Small helper to set an object datum given a function of coordinates
 *  \details Providing a function VecDof( VecDim ), this object assigns the
 *  result to a specified object according to its coordinate.
 *  \tparam GET  The function to get the datum
 *  \tparam SET  The objects mutator to set the datum
 */
template<typename GET, typename SET>
class corlib::SetData
{
    typedef typename GET::result_type         VecDof;
    typedef typename GET::argument_type       VecDim;
    typedef typename SET::first_argument_type ObjPtr;

public:
    SetData( GET get, SET set )
        : get_( get ), set_( set )
    { }

    //! Apply datum to object
    void operator()( ObjPtr objp )
    {
        VecDim X     = objp -> giveCoordinates();
        VecDof datum = get_(X);
        set_( objp, datum );

        return;
    }

private:
    GET  get_; //!< Functor to get the datum from
    SET  set_; //!< Functor to set the datum with
};

//------------------------------------------------------------------------------
/** \brief Small helper to set an object datum given an input stream
 *  \details Providing a stream of data, this object assigns the
 *  result to a specified object according to its ID.
 *  \tparam SET  The objects mutator to set the datum
 */
template<typename SET>
class corlib::SetDataFromFile
{
    typedef typename SET::first_argument_type  ObjPtr;
    typedef typename SET::second_argument_type VecDof;
private:
    typedef std::map<unsigned,VecDof>    DataMapType_;

public:
    SetDataFromFile( std::istream & inp, SET set )
        : set_( set )
    { 
        // number of data to expect
        unsigned numData;
        inp >> numData;
        inp.ignore( std::numeric_limits<std::streamsize>::max(), '\n' );
        
        // read each datum
        for ( unsigned i = 0; i < numData; i ++ ) {
            unsigned objID;
            inp >> objID;

            VecDof datum;
            for ( unsigned d = 0; d < VecDof::max_size; d ++ ) {
                inp >> datum[d];
            }
            inp.ignore( std::numeric_limits<std::streamsize>::max(), '\n' );
            dataMap_.insert( std::make_pair( objID, datum ) );
        }
    }

    //! Apply datum to object
    void operator()( ObjPtr objp )
    {
        unsigned objID = objp -> giveId();
        typename DataMapType_::iterator pos = dataMap_.find( objID );
        FTL_VERIFY_DESCRIPTIVE( pos != dataMap_.end(), "index for data not found: %d \n",  objID );

        VecDof datum = pos -> second;
        set_( objp, datum );

        return;
    }

    //! Remove storage
    void clear() { dataMap_.clear(); }

private:
    SET  set_;              //!< Functor to set the datum with
    DataMapType_ dataMap_;  //!< Storage of ID-data pairs
};

#endif
