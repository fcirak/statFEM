// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//
//                 Computational Structural Mechanics Lab
//                         University of Cambridge
//
// This software is copyrighted and all rights are retained by the University of
// Cambridge. All use, disclosure, and/or reproduction of any part not expressly
// authorised by the University of Cambridge is prohibited. (C) 2010
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

//! @file Distributor.hpp

#ifndef distributor_h
#define distributor_h

//------------------------------------------------------------------------------
#include <vector>
#include <type_traits>

#include <boost/function.hpp>
#include <boost/bind.hpp>
#include <Eigen/Core>

#include <corlib/verify.hpp>
#include <corlib/eigenX.hpp>

//------------------------------------------------------------------------------
namespace corlib{

    template<typename GIVSOL, typename GIVDOF, typename STOSOL> 
    class Distributor;

    //--------------------------------------------------------------------------
    //!Convenience function to provide a distributor using the function pointers
    template<typename VEC, typename SOLVER, typename NODE>
    corlib::Distributor<boost::function<void( const std::vector<unsigned> &,
                                              Eigen::VectorXd & )>,
                        boost::function<void( const NODE *,
                                              std::vector<unsigned> & )>,
                        boost::function<void( NODE *, const VEC &)> >
    distributorFun( const SOLVER * solver, 
                    void (NODE::*f)( std::vector<unsigned> & ) const,
                    void (NODE::*g)( const VEC & ),
                    const double factor = 1.0 )
    {
        FTL_VERIFY( solver not_eq NULL );

        //----------------------------------------------------------------------
        typedef boost::function<void( const std::vector<unsigned> &,
                                      Eigen::VectorXd & )>         GiveSolution;
        typedef boost::function<void( const NODE*,
                                      std::vector<unsigned> &)>         GiveDof;
        typedef boost::function<void( NODE *, const VEC & )>      StoreSolution;

        //----------------------------------------------------------------------
        GiveSolution  gs = boost::bind( &SOLVER::giveSolution, solver, _1, _2 );
        GiveDof       gd = f;
        StoreSolution ss = g;

        //----------------------------------------------------------------------
        return 
            corlib::Distributor<GiveSolution,GiveDof,StoreSolution>(gs, gd, ss, factor);
    }
}

//------------------------------------------------------------------------------
/** \brief Solution distribution object
 *  \details This object queries the nodal solution from the system based
 *  on a specific functor for the nodal dof indices (also given by a functor)
 *  and stores (using a functor) the solution at the node.
 *  \tparam GIVSOL  Solution giving functor
 *  \tparam GIVDOF  Dof index giving functor
 *  \tparam STOSOL  Solution storing functor
 */
template<typename GIVSOL, typename GIVDOF, typename STOSOL> 
class corlib::Distributor
{
public:
    typedef GIVSOL GiveSolution;
    typedef GIVDOF GiveDof;
    typedef STOSOL StoreSolution;

    // Object pointer (normal node) as first argument of the solution storage
    typedef typename StoreSolution::first_argument_type  NodePtr;
    // Type of nodal solution vector as second argument
    typedef typename StoreSolution::second_argument_type VecTypeRef;
    // Get rid of the reference qualifier
    typedef typename std::remove_reference<VecTypeRef>::type  Vec;

    //! Cstor given the functors and a scalar multiplier
    Distributor( GiveSolution  giveSolution,
                 GiveDof       giveDof, 
                 StoreSolution storeSolution,
                 const double  factor = 1. )
        : giveSolution_(  giveSolution  ), 
          giveDof_(       giveDof       ),
          storeSolution_( storeSolution ),
          factor_(        factor        ) {  }

    //! Overloaded function call acting on a node pointer
    void operator()( NodePtr np ) const
    {
        // get dof indices
        std::vector<unsigned> dofVec( Vec::RowsAtCompileTime );
        giveDof_( np, dofVec );

        // get solution from system
        static constexpr int size = Vec::RowsAtCompileTime;
        Eigen::VectorXd solVec( size );
        giveSolution_( dofVec, solVec );

        // convert and pass to node
        Vec sol( factor_ * solVec );
        storeSolution_( np, sol );

        return;
    }

private:
    GiveSolution  giveSolution_; //!< Solution giving functor
    GiveDof       giveDof_;      //!< Dof index giving functor
    StoreSolution storeSolution_;//!< Solution storing functor
    const double  factor_;       //!< Scalar multiplier
};

//------------------------------------------------------------------------------
#endif
