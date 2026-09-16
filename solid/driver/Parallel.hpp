// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//
//                 Computational Structural Mechanics Lab
//                         University of Cambridge
//
// This software is copyrighted and all rights are retained by the University of
// Cambridge. All use, disclosure, and/or reproduction of any part not expressly
// authorised by the University of Cambridge is prohibited. (C) 2010
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

//! @file Parallel.hpp

#ifndef solid_driver_parallel_h
#define solid_driver_parallel_h

#include <sstream>

// includes for parallel
#include <mpi.h>

// includes from corlib
#include <corlib/GenerateDofStarts.hpp>
#include <corlib/DofNumbererFromStarts.hpp>
#include <corlib/misc.hpp>
#include <corlib/AccumulateQuantity.hpp>
#include <corlib/VTUwriter.hpp>

#include <solid/driver/Static.hpp>

#include <tools/input/smf2vtu/readwrite.hpp>
#include <tools/input/partitioners/rcb/mesh/PartitionerRCB.hpp>

//==============================================================================
namespace solid {
    namespace driver {

        namespace ublas = boost::numeric::ublas;

        template< typename DRIVERSEQ > class ParallelBasic;

        /// Parallel driver which is specialised based on SOLVER type
        ///
        /// \tparam  DRIVERSEQ    Sequential driver type to wrap around
        /// \tparam  SOLVER       Repeated solver type to specialise
        template< typename DRIVERSEQ, typename SOLVER > class Parallel;
        
        /// Partition global mesh with RCB --- front end for tools/input/RCB
        ///
        /// \tparam ISTREAM   Input stream type
        /// \tparam OSTREAM   Output stream type
        ///
        /// \param[in]   myShape        Basic shape of element domain
        /// \param[in]   smf            SMF stream containing nodes and element connectivity
        /// \param[in]   lists          Optional node lists
        /// \param[in]   overlapLayers  Depth of overlapping (ghosting)
        /// \param[in]   numSub         Number of partitions
        /// \param[in]   iPart          Index of my partition (processor's rank)
        /// \param[out]  partSmf        My partitioned SMF file (ie PSMF)
        /// \param[out]  partLists      My optionally partitioned node lists
        template< typename ISTREAM, typename OSTREAM >
        static
        void partitionWithRCB( const enum corlib::shape myShape,
                               std::istream & smf,
                               std::vector<ISTREAM *> & lists,
                               const unsigned overlapLayers,
                               const int numPart,
                               const int iPart,
                               OSTREAM & partSmf,
                               std::vector<OSTREAM *> & partLists );


    }
}

//==============================================================================
/// Basic driver is for parallel solid computations
///
/// This object contains methods for driving in parallel shared by the
/// specfic solvers. 
///
/// \tparam DRIVERSEQ  Sequential driver
template< typename DRIVERSEQ >
class solid::driver::ParallelBasic : public DRIVERSEQ
{
private:
    typedef DRIVERSEQ                                         DriverSeq_;
    typedef solid::driver::ParallelBasic<DRIVERSEQ>           DriverParallelBasic_;

public:
    typedef typename DriverSeq_::Element                      Element;
    typedef typename DriverSeq_::Node                         Node;

    typedef typename DriverSeq_::Mesh                         Mesh;

    typedef typename DriverSeq_::Solver                       Solver;

    typedef typename DriverSeq_::VTUwriter                    VTUwriter;
    typedef corlib::UniqueFilename                            UniqueFilename;
    typedef corlib::VTUparaAnim                               VTUparaAnim;
    

public:
    //--------------------------------------------------------------------------
    /// @name Coarse level methods
    //@{
    
    /// Constructor with partitioned mesh
    ParallelBasic( std::istream & smf )
        : DriverSeq_( smf ),
          commAll_( MPI_COMM_WORLD ),
          animationPara_( NULL )
    { }

    /// Constructor with global mesh which is partitioned with RCB
    ParallelBasic( std::istream & smf, const unsigned overlapLayers )
        : DriverSeq_(),  // empty mesh
          commAll_( MPI_COMM_WORLD ),
          animationPara_( NULL )
    {
        std::stringstream subSmf;
        std::vector<std::stringstream *> lists;  // dummy
        std::vector<std::stringstream *> subLists;  // dummy
        this -> accessOutStream_()
            << "Partitioning mesh" << std::endl;
        solid::driver::partitionWithRCB( Element::myShape,
                                         smf, lists, overlapLayers,
                                         this -> getNumProcessors(),
                                         this -> getRank(),
                                         subSmf, subLists );
        mesh_.read( subSmf );  // fixed interface
    }

    /// Empty constructor read mesh later (with sequential driver interface)
    ParallelBasic()
        : DriverSeq_(),
          commAll_( MPI_COMM_WORLD ),
          animationPara_( NULL )
    { }

    /// D'tor
    virtual ~ParallelBasic()
    {
        if ( animationPara_ ) { this -> deAllocateAnimation(); }
    }

    //@}

    //--------------------------------------------------------------------------
    /// @name Fine level methods
    //@{


    /// Get number of processors in realm (MPI communicator)
    int getNumProcessors() const
    { 
        int nProc = -1;
        MPI_Comm_size( commAll_, &nProc );
        return nProc;
    }

    /// Get ID of processor
    int getRank() const
    {
        int rank = -1;
        MPI_Comm_rank( commAll_, &rank );
        return rank;
    }
        
    /// Wrapper to determine total number of DOFs
    ///
    /// \return Total number of DOFs
    unsigned generateDofs();

    ///
    std::pair<unsigned,unsigned> generateDofStarts( MPI_Comm commAll, std::map<int,int> & numNodeDofsGlobal )
    {
        return corlib::generateDofStarts( mesh_, commAll, numNodeDofsGlobal );
    }

    ///
    void dofNumberFromStartsFun( const std::map<int,int> & numNodeDofsGlobal )
    {
        mesh_.iterateOverNodes( corlib::dofNumberFromStartsFun( &Node::numberDOFs, &numNodeDofsGlobal ) );
    }

    /// Count number of local nodes
    int countLocalNodes( boost::function<bool( Node * )> condition = corlib::Positive() )
    {
        typedef corlib::AccumulateQuantity<Node,unsigned> Accumulator;
        Accumulator numNodesLocCounter( boost::bind( corlib::getOne ), 0 );
        const int numLocalNodes = mesh_.iterateOverNodesWithPredicate( numNodesLocCounter,
                                                                       condition );
        return numLocalNodes;
    }

    /// Clear global system of equations: matrix and RHS
    virtual void clearSystem()
    {
        this -> DriverSeq_::clearSystem();
        this -> clearSolution();
    }

    /// Create global residual displacement vector
    void clearSolution()
    {
        solver_ -> clearSol();
    }

    /// Inform user why/how linear system solve did (may be iterative solver)
    virtual void informAboutSolve() const
    {
        const unsigned numberOfIterations = solver_ -> giveNumIterations();
        const int convergedReason         = solver_ -> giveConvergedReason();
        this -> accessOutStream_() 
            << "    Solver converged in   " << numberOfIterations << " iterations"
            << " for reason   " << convergedReason << std::endl;
    }

    /// Enhance writing of element data with MPI ranks
    virtual void writeElementData( VTUwriter & vtuwriter ) const
    {
        const int rank = this -> getRank();
        corlib::AssignConstantNumber<const Element> procId( rank, "procID" );
        vtuwriter.writeElementQuantity( procId );

        //typedef boost::function< bool( const Element * )> BoolFun;
        //BoolFun giveDDActive = boost::bind( &Element::isDDActive, _1 );
        //corlib::Accessor< BoolFun > accessDDActive( giveDDActive, "IsDDActive" );
        //vtuwriter.writeElementQuantity( accessDDActive );

        this -> DriverSeq_::writeElementData( vtuwriter );
    }

    /// Create animation writer
    virtual void allocateAnimation( const std::string basename )
    {
        const int nProc = this -> getNumProcessors();
        animationPara_ = new VTUparaAnim( basename, "vtu", nProc );
    }

    /// Destroy animation writer
    virtual void deAllocateAnimation() 
    {
        delete animationPara_; animationPara_ = NULL;
    }

    /// Write step (ie VTU file)
    ///
    /// \param[in]  time  Time or characteristic load
    /// \param[in]  step  Time step or load step index
    virtual void writeAnimationStep( const double time, const unsigned step )
    {
        if ( animationPara_ ) {
            if ( this->getRank() == 0 ) animationPara_ -> storeSnapshots( time, step );
            const std::string vtuFile = UniqueFilename(3)( animationPara_->getBasename(),
                                                           this->getRank(), step, "vtu" );
            this -> writeMeshData( vtuFile );
        }
    }

    ///
    virtual void writeAnimationFile( const std::string animationFileName ) const
    {
        if ( animationPara_ ) {
            const std::string animationFile( animationFileName.c_str() );
            std::ofstream anim( animationFile.c_str( ) );
            animationPara_ -> writeAnimationFile( anim );
            anim.close();
        }
    }
    
    //@}

protected:
    /// Overload access to STDOUT enabling it only on master process
    virtual std::ostream & accessOutStream_() const
    {
        DriverParallelBasic_::coutDummy_.str( "" );  // empty dummy std out
        return ( this->getRank() == 0 ) ? std::cout : 
            DriverParallelBasic_::coutDummy_;
    }

protected:
    /// Solid mesh
    using DriverSeq_::mesh_;
    /// Linear algebraic solver
    using DriverSeq_::solver_;
    /// Used as total number of Dofs across all processors
    using DriverSeq_::numDofs_;

    /// MPI communicator
    MPI_Comm           commAll_;

    /// Animation writer for parallel environment
    VTUparaAnim *      animationPara_;

private:
    /// Dummy replacement for STDOUT on non-master processor
    static std::ostringstream coutDummy_;
};

template< typename DRIVERSEQ >
std::ostringstream solid::driver::ParallelBasic<DRIVERSEQ>::coutDummy_;

//------------------------------------------------------------------------------
#include "Parallel.ipp"

#endif
