// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//
//                 Computational Structural Mechanics Lab
//                         University of Cambridge
//
// This software is copyrighted and all rights are retained by the University of
// Cambridge. All use, disclosure, and/or reproduction of any part not expressly
// authorised by the University of Cambridge is prohibited. (C) 2010
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

//! @file ParallelPetsc.hpp

#ifndef solid_driver_parallelpetsc_h
#define solid_driver_parallelpetsc_h

#include <boost/static_assert.hpp>
#include <boost/type_traits.hpp>

#include <corlib/SystemSolvePetsc.hpp>
#include <solid/driver/Parallel.hpp>


//==============================================================================
namespace solid {
    namespace driver {

        namespace ublas = boost::numeric::ublas;

        /// Specialisation of parallel driver for PETSc 
        template< typename DRIVERSEQ >
        class Parallel<DRIVERSEQ,corlib::SystemSolvePetsc>;

    }
}

//==============================================================================
template< typename DRIVERSEQ >
class solid::driver::Parallel<DRIVERSEQ,corlib::SystemSolvePetsc>
    : public ParallelBasic<DRIVERSEQ>
{
private:
    typedef ParallelBasic<DRIVERSEQ>                  DriverParallelBasic_;

public:
    typedef typename DriverParallelBasic_::Node       Node;
    typedef typename DriverParallelBasic_::Mesh       Mesh;
    typedef typename DriverParallelBasic_::Solver     Solver;  // == corlib::SystemSolvePetsc
    BOOST_STATIC_ASSERT(( boost::is_same<Solver,corlib::SystemSolvePetsc>::value ));

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
        numLocalNodes_ = 0;
        numLocalDofs_  = 0;
    }

public:
    /// Numbering of DOFs
    virtual unsigned numberDofs()
    {
        numDofs_ = this -> DriverParallelBasic_::generateDofs();

        // count number of local nodes
        numLocalNodes_ = this -> DriverParallelBasic_::countLocalNodes( boost::bind( &Node::isDDActive, _1 ) );

        // number of local dofs
        numLocalDofs_ = numLocalNodes_ * Node::dof;

        //
        return numDofs_;
    }

    /// Allocate the linear solver WRT to #numLocalDofs_ of the processor
    virtual void allocateSolver()
    {
        if ( numDofs_ == 0 )
            this -> numberDofs();

        FTL_VERIFY( not solver_ );
        FTL_VERIFY( numDofs_ > 0 );
        FTL_VERIFY( numLocalDofs_ > 0 );
        solver_ = new Solver( numDofs_, numLocalDofs_ );
    }

    /// Create linear algebraic solver
    virtual void allocateSolver( const unsigned size )
    {
        FTL_VERIFY( false );
    }

protected:
    using DriverParallelBasic_::numDofs_;
    using DriverParallelBasic_::solver_;

    /// Number of local nodes
    int       numLocalNodes_;
    /// Number of local DOFs 
    unsigned  numLocalDofs_;

};

#endif
