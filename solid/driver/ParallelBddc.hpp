// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//
//                 Computational Structural Mechanics Lab
//                         University of Cambridge
//
// This software is copyrighted and all rights are retained by the University of
// Cambridge. All use, disclosure, and/or reproduction of any part not expressly
// authorised by the University of Cambridge is prohibited. (C) 2010
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

//! @file ParallelBddc.hpp

#ifndef solid_driver_parallelbddc_h
#define solid_driver_parallelbddc_h

// includes for parallel
#include <mpi.h>

#include <solid/driver/Parallel.hpp>
#include <corlib/SystemSolveBddc.hpp>


//==============================================================================
namespace solid {
    namespace driver {

        namespace ublas = boost::numeric::ublas;

        namespace detail_ {

            /// Hide call to constructor with a single unsigned
            template< >
            void allocateSolver<corlib::SystemSolveBddc>( corlib::SystemSolveBddc * & solver,
                                                          const unsigned numDofs )
            {
                FTL_VERIFY_DESCRIPTIVE( false, "This should never be called\n" );
            }

        }

        /// Specialisation of parallel driver for BDDCML
        template< typename DRIVERSEQ >
        class Parallel<DRIVERSEQ,corlib::SystemSolveBddc>;

    }
}

//==============================================================================
template< typename DRIVERSEQ >
class solid::driver::Parallel<DRIVERSEQ,corlib::SystemSolveBddc>
    : public solid::driver::ParallelBasic<DRIVERSEQ>
{
private:
    typedef ParallelBasic<DRIVERSEQ>                  DriverParallelBasic_;

public:
    typedef typename DriverParallelBasic_::Node       Node;
    typedef typename DriverParallelBasic_::Solver     Solver;  // == corlib::SystemSolveBddc
    BOOST_STATIC_ASSERT(( boost::is_same<Solver,corlib::SystemSolveBddc>::value ));

public:
    /// C'tor for partitioned mesh
    Parallel( std::istream & smf ) : DriverParallelBasic_( smf )
    { this -> initialise_(); }

    /// C'tor for complete mesh to be partitioned
    Parallel( std::istream & smf, const unsigned overlapLayers ) 
        : DriverParallelBasic_( smf, overlapLayers )
    { this -> initialise_(); }

    /// Empty c'tor to read in mesh later
    Parallel() : DriverParallelBasic_()
    { this -> initialise_(); }

protected:
    /// Common initialisation of member variables
    void initialise_()
    {
        numNodesSub_ = 0;
        numDofsSub_  = 0;
    }

public:
    /// Numbering of DOFs
    virtual unsigned numberDofs()
    {
        numDofs_ = this -> generateDofs();

        // count number of nodes on subdomain
        numNodesSub_ = this -> countLocalNodes();

        // number of subdomain dofs
        numDofsSub_  = numNodesSub_ * Node::dof;

        //
        return numDofs_;
    }

    virtual void allocateSolver()
    {
        if ( numDofs_ == 0 )
            this -> numberDofs();

        // specify type of matrix ( GENERAL, SPD, or SYMMETRICGENERAL )
        const typename Solver::MatrixType matrixType = Solver::GENERAL;

        FTL_VERIFY( numDofs_ > 0 );
        FTL_VERIFY( numDofsSub_ > 0 );
        FTL_VERIFY( not solver_ );
        solver_ = new Solver( numDofs_, numDofsSub_, matrixType );

        // load subdomain mesh into BDDCML solver
        solver_ -> loadMesh( mesh_ );
    }
    
    virtual void allocateSolver( const unsigned size )
    {
        FTL_VERIFY( false );
    }

    /// Finalise assembly of global system
    virtual void finishAssembly() { ; }  // do nothing

private:
    using DriverParallelBasic_::mesh_;
    using DriverParallelBasic_::solver_;
    using DriverParallelBasic_::numDofs_;
    
    int       numNodesSub_;
    unsigned  numDofsSub_;
    
};

#endif
