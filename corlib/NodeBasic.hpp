// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//
//                 Computational Structural Mechanics Lab
//                         University of Cambridge
//
// This software is copyrighted and all rights are retained by the University of
// Cambridge. All use, disclosure, and/or reproduction of any part not expressly
// authorised by the University of Cambridge is prohibited. (C) 2010
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

//! @file NodeBasic.hpp

#ifndef corlib_nodebasic_h
#define corlib_nodebasic_h
//------------------------------------------------------------------------------

#include <corlib/eigenX.hpp>

namespace corlib{

    template<unsigned DIM> class NodeBasic;

}

//------------------------------------------------------------------------------
/** \brief Base class for nodes
 *  \details Stores an ID-number and the coordinates
 *  \tparam DIM the dimension
 */
template<unsigned DIM>
class corlib::NodeBasic
    : boost::noncopyable
{
public:
    //! The spatial dimension of this node
    static const unsigned dim = DIM;   

    //! Coordinate array
    typedef eigenX::VectorSd<dim> VecDim;

    //--------------------------------------------------------------------------
    //! @name Construction
    //@{

    /** \brief Empty constructor 
     *  \details Invalidates the stored data
     */
    NodeBasic()
        : id_( std::numeric_limits<unsigned>::max() ) 
    { 
        coord_.fill( std::numeric_limits<double>::max() );
    }

    //! Read coordinates from a stream  \param[in] inp Input stream
    std::istream & readSelf( std::istream & inp )
    {
        // read dim coordinates
        for ( unsigned i = 0; i < dim; i++ ) { 
            inp >> coord_[ i ];
            if( inp.peek( ) == ',' ) inp.ignore( );
        }
        // swallow 3-dim coordinates
        for ( unsigned i = dim; i < 3; i++ ) {
            double dummy; 
            inp >> dummy;  
        } 
        return inp;
    }
    //@}
    
    //--------------------------------------------------------------------------
    //! @name Coordinates
    //@{

    /** Set coordinates \details
     *  \param[in] c  Vector containing the coordinates  */
    void setCoordinates( const VecDim & c ) { coord_ = c; }

    /** Give the coordinates \details
        \return The coordinates of this node */
    VecDim giveCoordinates( ) const { return coord_; };
    //@}

    //--------------------------------------------------------------------------
    //! @name ID
    //@{
    
    /** Give node's ID \details 
     * \return       The ID of this node */
    unsigned giveId() const { return id_; }

    /** Set the ID of this node \details
     *  \param[in] i  New ID of this node */
    void setId( const unsigned & i ) { id_ = i; return; }
    //@}
    

protected:
    //--------------------------------------------------------------------------
    unsigned        id_;           //!< identifiying number
    VecDim       coord_;           //!< coordinates
};

#endif
