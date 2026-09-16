// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//
//                 Computational Structural Mechanics Lab
//                         University of Cambridge
//
// This software is copyrighted and all rights are retained by the University of
// Cambridge. All use, disclosure, and/or reproduction of any part not expressly
// authorised by the University of Cambridge is prohibited. (C) 2010
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

//! @file SmfHead.hpp

#ifndef corlib_smfhead_h
#define corlib_smfhead_h

//------------------------------------------------------------------------------
// headers
#include <corlib/PropertiesParser.hpp>
#include <corlib/Shape.hpp>
#include <corlib/verify.hpp>
#include <iostream>
#include <limits>

//------------------------------------------------------------------------------
namespace corlib {
    class SmfHead;
}


//------------------------------------------------------------------------------
/** \brief Skims head of SMF input file and extract contained variables
 *  \details This object allows to read the head of a SMF geometry file.
 *  The head of a SMF geometry file consists of lines starting with
 *  <b>!</b> (an exclamation mark).
 *
 *  Two variables are searched for in the header:
 *     - <i>elementShape</i>: the geometric element shape, i.e. <i>line</i>,
 *                            <i>triangle</i>, etc
 *     - <i>elementNumPoints</i>: number of points (nodes) per element
 *
 *  Example 1 (plain header):
 *\verbatim
 *     ! elementShape line
 *     ! elementNumPoints 2
 *     8 8
 *     ... ...
 *\endverbatim
 *
 *
 *  Example 2 (header and comments):
 *\verbatim
 *     # My triangle mesh of a rectangle with 4x6 nodes
 *     # etc pp ...
 *     ! elementShape triangle
 *     ! elementNumPoints 3
 *     # A comment, which continues
 *     #            on the next line
 *     24 30
 *     ... ...
 *\endverbatim
 */
class corlib::SmfHead
{
public:

    //! Constructor
    SmfHead( );

    //! Destructor
    virtual ~SmfHead( ) { }

    //! Try to read the header
    //!
    //! \param[in,out] is  input stream
    //! \return            input stream without header
    std::istream & read( std::istream & is );

    //! Read and validate quantities
    //!
    //! \param[in,out]  is   input stream
    //! \param[in]      se   shape of element to validate
    //! \param[in]      np   number of element point to validate
    //! \return              input stream without header
    std::istream & readValidated( std::istream & is,
                                  const enum shape se,
                                  const unsigned np );

    //! Found a header section
    //!
    //! \return   true if at least on line starts is indicated as header line
    bool foundHeader( ) const { return foundHeader_; }

    //! Check if registered variables have been found
    //!
    //! \return   true if all registered variables have been found
    bool foundVariables( ) const;
    
    //! element shape type
    enum shape giveElementShape( ) const { return convertShapeStringToEnum( elementShape_ ); }

    //! set element shape
    //!
    //! \param[in]   es   element shape to set
    void setElementShape( const enum shape es ) { elementShape_ = convertShapeEnumToString( es ); return; }

    //! number of element nodes
    unsigned giveElementNumPoints( ) const { return elementNumPoints_; }

    //! set number of element nodes
    void setElementNumPoints( const unsigned np ) { elementNumPoints_ = np; return; }

    //! write header to stream
    //!
    //! \param[in,out]    os    output stream with added SMF header
    void write( std::ostream & os ) const;

    //! write header to stream
    //!
    //! \param[in]        es    shape of element
    //! \param[in]        np    number of element point
    //! \param[in,out]    os    output stream with added SMF header
    void write( const std::string es,
                const unsigned np,
                std::ostream & os ) const;

    //! write header to stream
    //!
    //! \param[in]        es    shape of element
    //! \param[in]        np    number of element point
    //! \param[in,out]    os    output stream with added SMF header
    void write( const enum shape es,
                const unsigned np,
                std::ostream & os ) const;

private:
    //! Character indicating a header line
    const char           headChar_;

    //! Character indicating a comment line
    const char           commentChar_;

    //! Found a header in SMF file
    bool                 foundHeader_;

    //! Parser to read in header data
    PropertiesParser     parser_;

    //! @name Header data
    //@{
    //! Input identifier for element shape
    const std::string    elementShapeName_;
    //! Shape of element
    std::string          elementShape_;

    //! Input identifier for number of nodes per element
    const std::string    elementNumPointsName_;
    //! Number of nodes per element
    unsigned             elementNumPoints_;
    //@}
    
};

//------------------------------------------------------------------------------
// implementations
#include "SmfHead.ipp"

#endif
