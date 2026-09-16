// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//
//                 Computational Structural Mechanics Lab
//                         University of Cambridge
//
// This software is copyrighted and all rights are retained by the University of
// Cambridge. All use, disclosure, and/or reproduction of any part not expressly
// authorised by the University of Cambridge is prohibited. (C) 2010
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

//! @file TimingProfiler.hpp

//! @todo   - add a constructor argument to set the number of decimal places
//!           to print
//!         - print a warning when the timing length too sort for the chosen
//!           accuracy to make sense
//!         - need a table that can hold anything, not just doubles
//!         - Beautify the output, search for the longest number and then insert
//!           the right number of spaces
//!         - Make sure you write out values in order they're added! It's a bit
//!           confusing when that all changes!

#ifndef aux_timingprofiler_h
#define aux_timingprofiler_h

//------------------------------------------------------------------------------
// Includes
#include <map>
#include <string>
#include <iostream>
#include <boost/timer.hpp>
#include <boost/algorithm/string/predicate.hpp>
#include <corlib/verify.hpp>

//------------------------------------------------------------------------------
// Declarations
namespace aux {

    class TimingProfiler;

    namespace detail_{
        //! Case insensitive string comparison functor
        struct ciStringCompare
        {
            bool operator()(const std::string& s1, const std::string& s2) const
            {
                return boost::algorithm::ilexicographical_compare( s1, s2 );
            }
        };
    }
}

//------------------------------------------------------------------------------
/** \brief Timer class stores the variables being timed and times each of them
 *  between start and stop calls. Writes a nicely formatted output file. Usage:
 *
 *  aux::TimingProfiler prof;
 *  ...
 *  prof.record( "degree", degree );      // record a variable
 *  ...
 *  prof.startTimer( "Solve" );           // start timng a procedure
 *  ... <solve system> ...
 *  prof.stopTimer( "Solve" );            // end timing a procedure
 *  ...
 *  std::ofstream timing( "timing.dat" );
 *  FTL_VERIFY( timing.is_open() );
 *  prof.writeOutput( timing );           // write results to a file
 *  timing.close( );
 *
 */

class aux::TimingProfiler
{
public:
    typedef std::map< std::string, double, detail_::ciStringCompare > Table;
    typedef Table::const_iterator TableIter;

public:

    // Constructor
    TimingProfiler( ) : stopwatch_() { }

    // Destructor
    ~TimingProfiler( ) {
        FTL_VERIFY_DESCRIPTIVE( starts_.size() == times_.size(), "Did not stop all timers " );
    }

    // Record a parameter/variable
    void record( const std::string name, const double value  ){
        FTL_VERIFY_DESCRIPTIVE( recorded_.count( name ) == 0, "Tried to duplicate variable " );
        recorded_[ name ] = value;
        return;
    }

    // Start Timer method
    void startTimer( const std::string name ){
        FTL_VERIFY_DESCRIPTIVE( starts_.count( name ) == 0, "Tried to duplicate timer " );
        starts_[ name ] = stopwatch_.elapsed();
        return;
     }

    // Stop Timer method
    void stopTimer( const std::string name ){
        FTL_VERIFY_DESCRIPTIVE( starts_.count( name ) != 0, "Tried to stop non-existant timer " );
        times_[ name ] = stopwatch_.elapsed() - starts_[ name ];
        return;
    }

    // Write output
    void writeOutput( std::ostream & out, const bool header = true ) {
        TableIter iter;

        // write header
        if( header ) {
            writeHeader_( out, recorded_ );
            writeHeader_( out, times_ );
            out << std::endl;
        }

        // write data
        writeData_( out, recorded_ );
        writeData_( out, times_ );
        out << std::endl;

        return;
    }

private:
    // Write table header
    void writeHeader_( std::ostream & out, Table table ){
        TableIter iter;

        iter = table.begin();
        for ( ; iter != table.end(); ++ iter )
            out << "# " << (*iter).first << " ";

        return;
    }

    // write table data
    void writeData_( std::ostream & out, Table table ){
        TableIter iter;

        iter = table.begin();
        for ( ; iter != table.end(); ++ iter )
            out << (*iter).second << " ";

        return;
    }


private:
    boost::timer stopwatch_;
    Table recorded_, starts_, times_;
};

#endif // aux_timingprofiler_h
