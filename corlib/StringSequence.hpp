// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//
//                 Computational Structural Mechanics Lab
//                         University of Cambridge
//
// This software is copyrighted and all rights are retained by the University of
// Cambridge. All use, disclosure, and/or reproduction of any part not expressly
// authorised by the University of Cambridge is prohibited. (C) 2010
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

//! @file StringSequence.hpp

#ifndef corlib_stringsequence_h
#define corlib_stringsequence_h

#include <string>
#include <vector>
#include <sstream>
#include <corlib/PropertiesParser.hpp>
#include <corlib/verify.hpp>

namespace corlib {

    //==========================================================================
    /// Object for reading in a sequence of strings interspersed with comments
    ///
    /// The following example string sequences
    ///\verbatim
    ///    { hello world }
    ///    { somewhere
    ///      under the
    ///      rainbow
    ///    }
    ///    { # a collection of flags
    ///      -h world       # you need to set this one
    ///      -e someFlag
    ///      -32            # it's 32 bit
    ///      # end of flags }
    ///\endverbatim
    /// lead to 
    ///\verbatim
    ///    "hello world"
    ///    "somewhere under the rainbow"
    ///    "-h world -e someFlag -32"
    ///\endverbatim
    class StringSequence
    {
    public:
        /// Insertion operator to read string sequence from stream
        friend
        std::istream & operator>>( std::istream & is, StringSequence & strSeq )
        {
            // check first character
            char c;
            is >> c;
            FTL_VERIFY( c == StringSequence::sequenceBegin_ );

            // copy all text until reaching #sequenceEnd_ into #inputSeq
            std::stringstream inputSeq;
            is.get(c);
            while ( c != StringSequence::sequenceEnd_ ) {
                inputSeq << c;
                is.get(c);
            }

            // loop #inputSeq and skim off comments
            // read words are placed in #sequence_ separated by spaces
            while ( not inputSeq.eof() ) {
                corlib::skip_comment( inputSeq );
                std::string str;
                inputSeq >> str;
                strSeq.sequence_.push_back( str );
            }

            return is;
        }

        /// Extraction operator to write string sequence to output stream
        friend
        std::ostream & operator<<( std::ostream & os, const StringSequence & strSeq )
        {
            os << strSeq();
            return os;
        }

        /// Return string containing concatenated sequence of strings
        std::string operator()() const
        {
            std::ostringstream concat;
            std::copy( sequence_.begin(), --sequence_.end(),
                       std::ostream_iterator<std::string>( concat, " " ) );
            concat << *(--sequence_.end());
            return concat.str();
        }

    private:
        /// Character opening a sequence of strings
        static const char sequenceBegin_ = '{';
        /// Character closing a sequence of strings
        static const char sequenceEnd_   = '}';
        /// The sequence of strings
        std::vector<std::string> sequence_;
    };

}

#endif
