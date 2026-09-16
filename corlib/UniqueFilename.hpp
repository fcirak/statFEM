// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//
//                 Computational Structural Mechanics Lab
//                         University of Cambridge
//
// This software is copyrighted and all rights are retained by the University of
// Cambridge. All use, disclosure, and/or reproduction of any part not expressly
// authorised by the University of Cambridge is prohibited. (C) 2010
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

//! @file UniqueFilename.hpp
//! @todo Decide how many operator()s are really needed

#ifndef corlib_uniquefilename_h
#define corlib_uniquefilename_h
//------------------------------------------------------------------------------
#include <string>
#include <sstream>
#include <iomanip>
//------------------------------------------------------------------------------
namespace corlib{
    class UniqueFilename;
}    

//------------------------------------------------------------------------------
/** \brief Generate filename for parallel and/or dynamic VTU files
 *  \details Tool to generate filenames of the type 
 *  \a basename.XXX-YYYY.suffix or \a basename.XXX.suffix
 *  with \a XXX refering to the first (typically processor) number and 
 *  \a YYYY to the second (typically time step) number. Standard width of 
 *  the first number is 3 and of the second is 4, but both can be adjusted.
 */
class corlib::UniqueFilename
{
public:
    //! Constructor given the widths of the number strings 
    UniqueFilename( const unsigned width1 = 3,
                    const unsigned width2 = 4,
                    const unsigned width3 = 4 ) 
        : width1_( width1 ), width2_( width2 ), width3_( width3 )  { }

    //! Function call operator which generates the filename string (one number)
    std::string operator()( const std::string & basename,
                            const unsigned      number,
                            const std::string & suffix )
    {
        std::string numberString = this -> convertNumberToString( number, width1_ );
        std::string filename = basename + "." + numberString + "." + suffix;
        return filename;
    }

    //! Function call operator which generates the filename string (two numbers)
    std::string operator()( const std::string & basename,
                            const unsigned      number1,
                            const unsigned      number2,
                            const std::string & suffix )
    {
        std::string numberString1 = this -> convertNumberToString( number1, width1_ );
        std::string numberString2 = this -> convertNumberToString( number2, width2_ );
        std::string filename = basename + "." + numberString1 + "-" + numberString2
            + "." + suffix;
        return filename;
    }

    //! Function call operator which generates the filename string (three numbers)
    std::string operator()( const std::string & basename,
                            const unsigned      number1,
                            const unsigned      number2,
                            const unsigned      number3,
                            const std::string & suffix )
    {
        std::string numberString1 = this -> convertNumberToString( number1, width1_ );
        std::string numberString2 = this -> convertNumberToString( number2, width2_ );
        std::string numberString3 = this -> convertNumberToString( number3, width3_ );
        std::string filename = basename + "." + numberString1 + "." + numberString2 + "-" 
                                        + numberString3 + "." + suffix;
        return filename;
    }

public:
    //! Conversion from a number to a string of prescribed width with leading zeros
    static std::string convertNumberToString( const unsigned number,
                                              const unsigned width ) {
        std::stringstream pipe;
        pipe << std::setfill( '0' ) << std::setw( width ) << number << std::endl;
        std::string numberString;
        pipe >> numberString;
        return numberString;
    }

private:
    const unsigned width1_; //!< Field width of the first number string
    const unsigned width2_; //!< Field width of the second number string
    const unsigned width3_; //!< Field width of the third number string
};

//------------------------------------------------------------------------------
namespace corlib{

    //! convenience function for vtf-users
    std::string makeUniqueName( const int step, const int rank,
                                const char * first, const char * second )
    {
        std::string basename( first );
        std::string suffix( second );
        corlib::UniqueFilename uniqueFilename( 4, 6 );
        return uniqueFilename( basename, static_cast<unsigned>( rank ),
                               static_cast<unsigned>( step ), suffix );
    }
}

#endif
