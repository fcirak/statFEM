// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//
//                 Computational Structural Mechanics Lab
//                         University of Cambridge
//
// This software is copyrighted and all rights are retained by the University of
// Cambridge. All use, disclosure, and/or reproduction of any part not expressly
// authorised by the University of Cambridge is prohibited. (C) 2010
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

//! @file VectorAssembler.hpp

#ifndef corlib_vectorassembler_h
#define corlib_vectorassembler_h

//------------------------------------------------------------------------------
//! system includes
#include <vector>
//! boost includes
#include <boost/function.hpp>
//! Eigen includes
#include <Eigen/Core>
//! corlib includes
#include <corlib/verify.hpp>
#include <corlib/eigenX.hpp>

//------------------------------------------------------------------------------
namespace corlib{

    template<typename GIVVEC, typename GIVDOF, typename SOLVER>
    class VectorAssembler;

    //--------------------------------------------------------------------------
    //! Convenience function to generate a vector assembler
    template<typename VEC, typename NODE, typename SOLVER>
    corlib::VectorAssembler< boost::function<VEC(  const NODE *)>, 
                             boost::function<void( const NODE *, 
                             std::vector<unsigned> & )>,
                             SOLVER>
    vectorAssemblerFun( VEC ( NODE::*f )() const,
                        void( NODE::*g )( std::vector<unsigned> & ) const,
                        SOLVER * solver, 
                        const double factor = 1. )
    {
        //----------------------------------------------------------------------
        typedef boost::function<VEC(  const NODE *)>                    GiveVec;
        typedef boost::function<void( const NODE *, 
                                      std::vector<unsigned> & )>        GiveDof;

        //----------------------------------------------------------------------
        GiveVec gv = f;
        GiveDof gd = g;
        
        //----------------------------------------------------------------------
        // return vector assembler object
        return corlib::VectorAssembler<GiveVec,GiveDof,
                                       SOLVER>( gv, gd, solver, factor );
    }
}

//------------------------------------------------------------------------------
/** \brief Assembly of a local vector to the system's vector
 *  \details By means of two functors which provide the local vector and the
 *  corresponding degree of freedom indices, the object is assembled into the
 *  system vector using an insertion functor.
 *  \tparam GIVVEC Type of vector giving functor
 *  \tparam GIVDOF Type of dof index giving functor
 *  \tparam SOLVER Type of solver receiving the vector
 */
template<typename GIVVEC, typename GIVDOF, typename SOLVER> 
class corlib::VectorAssembler
{
public:
    //! @name Template parameters
    //@{
    typedef      GIVVEC             GiveVector;
    typedef      GIVDOF             GiveDof;
    typedef      SOLVER             Solver;
    //@}

    //! Extract node pointer type from functor
    typedef typename GiveVector::argument_type NodePtr;

    //! Cstor with initialisers for all private data
    VectorAssembler( GiveVector   giveVector,
                     GiveDof      giveDof,
                     Solver * solver, 
                     const double factor = 1. )
        : giveVector_( giveVector ),
          giveDof_(    giveDof    ),
          solver_(     solver     ),
          factor_(     factor     )
    { 
        FTL_VERIFY( solver_ not_eq NULL );
    }
    
    //----------------------------------------------------------------------
    /** Overloaded function call to perform the local->global assembly
     *  \param[in] np  Pointer to node which delivers the vector
     */
    void operator()( const NodePtr np ) const
    {
        //! - get corresponding dof indices
        std::vector<unsigned> dofs;
        giveDof_( np, dofs  );

        //! - retrieve local vector from object
        Eigen::VectorXd forceVec = giveVector_( np );
        forceVec *= factor_;

        //! - pass on to system solver 
        solver_ -> insertToRhs( forceVec, dofs );
        return;
    }

private:
    //--------------------------------------------------------------------------
    GiveVector     giveVector_;   //!< Functor giving the local vector
    GiveDof        giveDof_;      //!< Functor giving the dof indices
    Solver       * solver_;       //!< Pointer to system solver
    const double   factor_;       //!< Scalar multiplier
};


#endif
