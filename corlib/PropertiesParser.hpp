// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//
//                 Computational Structural Mechanics Lab
//                         University of Cambridge
//
// This software is copyrighted and all rights are retained by the University of
// Cambridge. All use, disclosure, and/or reproduction of any part not expressly
// authorised by the University of Cambridge is prohibited. (C) 2010
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

//! @file PropertiesParser.hpp

#ifndef corlib_propertiesparser_h
#define corlib_propertiesparser_h
#include <map>
#include <list>
#include <string>
#include <iostream>
#include <algorithm>
#include <limits>
#include <iterator>
#include <boost/bind.hpp>
#include <boost/algorithm/string.hpp>
#include <corlib/misc.hpp>

//------------------------------------------------------------------------------
// declarations
namespace corlib {
    class PropertiesParser;

    namespace detail_{

        class SkipCommentsIstream;
        class Mutator;
        template <typename T> class TypedMutator;
    }

    std::istream & skip_comment( std::istream & in ); 
}

//! \cond  SKIPDOX

//------------------------------------------------------------------------------
// auxiliary class for skipping comments
namespace corlib{
    namespace detail_{
        
        //! Stream reads comment lines from the stream
        class SkipCommentsIstream 
        {
        public:
            //! Constructor with comment chars
            SkipCommentsIstream( char comment_begin = '#',
                                 char comment_end   = '\n') 
                : comment_begin_( comment_begin ), 
                  comment_end_( comment_end ) { }
	
            //! Main function which reads the comment line
            std::istream & operator()( std::istream & in ) const
            {
                while( true ) {
                    in >> std::ws;
                    if (in.peek() != comment_begin_) break;
                    in.get();
                    in.ignore( std::numeric_limits<std::streamsize>::max(), 
                               comment_end_ );
                }
                return in;
            } 
            
        private:
            const char comment_begin_; //!< Begin of a comment line
            const char comment_end_;   //!< End of a comment line
        };        
    }
}

//------------------------------------------------------------------------------
namespace corlib{
    // Convenience function for use of #SkipCommentsIstream
    std::istream & skip_comment( std::istream & in ) 
    {
        if( in ) {
            detail_::SkipCommentsIstream()( in );
        }
        return in;
    }
}

//------------------------------------------------------------------------------
namespace corlib {
    namespace detail_{
     
        //! Case insensitive string comparison functor
        struct ciStringCompare 
        { 
            bool operator()( const std::string& s1, const std::string& s2 ) const 
            {
                return boost::ilexicographical_compare( s1, s2 );
            }
        };
    }
}


//------------------------------------------------------------------------------
// Encapsulation of the stored references
namespace corlib{
    namespace detail_{

        /** Abstract base class for the typed mutator which encapsulates 
         *  the datum to be read         */
        class Mutator 
        {
        public:
            virtual void read ( std::istream& in )  = 0;
            virtual void print( std::ostream& out, 
                                const std::string& name ) const = 0;
            virtual bool isRead( ) const = 0;
            virtual ~Mutator( ) { }
        };

        //! Specific encapsulation of the datum with overloaded read and print 
        template<typename T>
        class TypedMutator : public detail_::Mutator 
        {
        public:
            TypedMutator( T& v ) : v_( v ), isRead_( false ) { }
            ~TypedMutator( ) { }
            void read( std::istream & in ) { in >> v_; isRead_ = true; }
            void print( std::ostream & out, const std::string & prefix ) const 
            {
                out << prefix << v_;
            }
            
            bool isRead( ) const { return isRead_; }
            
        private:
            T        &v_;      //!< Reference to datum read from stream
            bool     isRead_;  //!< True if has been read from the stream
        };
    }
} 

//!\endcond

//------------------------------------------------------------------------------
/** \brief Stream parser which searches for registered variable names and values
 *  \details Given an input stream, this object skips comments and tries to find
 *  lines with the name of a variable and its value in the stream. The variable
 *  names and references to their values are registered beforehand by the caller.
 */
class corlib::PropertiesParser 
{
public:
    //! Empty constructor
    PropertiesParser( ) { }
    
    //! Destructor deallocates dynamic memory
    ~PropertiesParser( )
    {
        // delete values of table
        std::for_each( table_.begin(), table_.end(), 
                       boost::bind( corlib::deleteFunctor(), 
                                    boost::bind( &TableType_::value_type::second, 
                                                 _1 ) ) );

        table_.clear( );
        unrecognized_.clear( );
    }

    /** Register a variable to be read from the stream
     *  \tparam T  Type of datum to be read
     *  \param[in] name   Name of the variable
     *  \param[in] t      Reference to the variable
     */
    template<class T>
    void registerPropertiesVar(const std::string& name, T& t)
    { 
        // If name is already in the table, create a warning message
        if ( table_.find(name) != table_.end( ) ) {
            std::cerr << "Properties parser called multiple times for " 
                      << name << std::endl;
        }
        else {
            // Create new encapsulation for this variable
            detail_::TypedMutator<T>* p = new detail_::TypedMutator<T>( t );
            // store in table
            table_[ name ] = p;
        }
    }

    /** Read the values of the stored variables from the input stream.
     *  Comments are first skipped, then the internal read function is called.
     *  \param[in,out] is  Input stream
     */
    void readValues( std::istream& is )
    {
        while ( is ) { 
            corlib::skip_comment( is );
            this -> readVariable( is ); 
        }
    }

    //! Read a (non-comment) line from the stream
    void readVariable( std::istream& is )
    {
        is >> std::ws;
        if( is ) {
            // get name of variable
            std::string s; 
            is >> s;
	
            if( table_.find(s) != table_.end() ) {
                // if recognised, read this variable
                table_[s] -> read( is );
            }
            else if ( !s.empty() && s != "\n") {
                // unrecognised is stored in extra list
                unrecognized_.push_back( s );
            }

        }
    }

    /** Print all the values from the table
     *  \param[in,out] out  Output stream
     *  \param[in]     pre  Prefix to variable name
     *  \param[in]     sep  Separator between variables
     */
    void printValues( std::ostream& out, const std::string& pre, 
                      const std::string& sep ) const
    {    
        TableType_::const_iterator item( table_.begin() );
        // go through all items in the table
        for( ; item != table_.end(); ++ item ) {
            (*item).second -> print( out, pre + (*item).first + sep );
            out << "\n";
        }
    }
    
    //! Predicate if unrecognised variables have been found
    bool hasUnrecognized( ) const
    {
        return ( unrecognized_.size() != 0 );
    }
    
    //! Print all unrecognised variables \param[in,out] out Output stream
    void printUnrecognized( std::ostream & out ) const
    {
        std::copy( unrecognized_.begin(), unrecognized_.end(), 
                   std::ostream_iterator<std::string>( out, "\n" ) );
    }

    //! Check if all variables have been read
    bool isEverythingRead( ) const {
        // check if number of read variables equals table size
        return ( static_cast<int>(table_.size()) ==
                 std::count_if( table_.begin(), table_.end(), 
                                boost::bind( &detail_::Mutator::isRead, 
                                             boost::bind( &TableType_::value_type::second, _1 ) ) ) );
    }
    
    //! Write out the unread variables \param[in,out] out Output stream
    void writeUnread( std::ostream & out ) const {
        // go through table and check items' states
        TableType_::const_iterator item = table_.begin();
        for ( ; item != table_.end(); ++ item ) {
            if ( not (item -> second -> isRead() ) ) {
                out << item -> first << std::endl;
            }
        }
    }
    
        
private:
    typedef std::map<std::string, detail_::Mutator*, 
                     detail_::ciStringCompare>       TableType_;
    typedef std::list<std::string>                   UnrecognizedList_;
    
    TableType_        table_;         //!< Storage of variable <name,reference> pairs
    UnrecognizedList_ unrecognized_;  //!< Storage of variables not recognised
};

#endif
