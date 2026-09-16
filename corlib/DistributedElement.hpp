// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//
//                 Computational Structural Mechanics Lab
//                         University of Cambridge
//
// This software is copyrighted and all rights are retained by the University of
// Cambridge. All use, disclosure, and/or reproduction of any part not expressly
// authorised by the University of Cambridge is prohibited. (C) 2010
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

//! @file DistributedElement.hpp

#ifndef corlib_distributedelement_h
#define corlib_distributedelement_h
//------------------------------------------------------------------------------

namespace corlib{
    template<typename ELEMENT> class DistributedElement;
}

//------------------------------------------------------------------------------
/** \brief  Element wrapper for distributed computation 
 *  \details This object equips the ELEMENT it wraps around with an activity 
 *  flag.
 *  \tparam ELEMENT Type of element to wrap around
 */
template<typename ELEMENT>
class corlib::DistributedElement : public ELEMENT
{
public:
    //! Read from stream
    template<typename NODECONT>
    std::istream & readSelf( std::istream & inp, const NODECONT & nc ) {
        //! - let base class read its connectivity
        ELEMENT::readSelf( inp, nc );
        //! - read activity
        int ddActive;
        inp >> ddActive;
        if ( ddActive == 1 ) {
            ddState_ = DDACTIVE;
        }
        else {
            ddState_ = DDOVERLAP;
        }

        return inp;
    }

    //! switch on activity
    void setDDActive() { ddState_ = DDACTIVE; }

    //! query whether element is active
    bool isDDActive() const { return ddState_ == DDACTIVE; }

    //! switch on activity in assembly
    void setDDOverlap() { ddState_ = DDOVERLAP; }

    //! query whether element is in overlap
    bool isDDOverlap() const { return ddState_ == DDOVERLAP; }

    //! query whether element should be assembled (active or overlap)
    bool isDDAssembled() const { return ddState_ == DDACTIVE or ddState_ == DDOVERLAP; }

    //! dissable element
    void setDDInactive() { ddState_ = DDINACTIVE; }

    //! return whether element is deactivated
    bool isDDInactive() const { return ddState_ == DDINACTIVE; }

private:
    enum DDState_ { 
        DDACTIVE,      //!< Active element - belongs to the subdomain in the nonoverlapping division
        DDOVERLAP,     //!< Element belong to another subdomain, but its contribution is necessary for
                       //!<  assembling of local rows for degrees of freedom
        DDINACTIVE     //!< Element exists on the subdomain but does not contribute much
    };
    DDState_ ddState_; //!< state of the element 
};

#endif
