// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//
//                 Computational Structural Mechanics Lab
//                         University of Cambridge
//
// This software is copyrighted and all rights are retained by the University of
// Cambridge. All use, disclosure, and/or reproduction of any part not expressly
// authorised by the University of Cambridge is prohibited. (C) 2010
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

//! @file   misc.hpp

#ifndef corlib_misc_h
#define corlib_misc_h
//------------------------------------------------------------------------------
#include <cstdlib>
#include <cassert>
#include <utility>
#include <functional>
#include <iomanip>
#include <sstream>
#include <boost/function.hpp>
#include <corlib/verify.hpp>
//------------------------------------------------------------------------------


namespace corlib {
    template<typename INP, typename FUN, typename PRED>
    FUN for_each_if( INP first, INP last, FUN f, PRED p );

    //! Indicate the side of the basis transform
    enum Side{ LEFT, RIGHT, NONE };
}

//------------------------------------------------------------------------------

namespace corlib {
    //--------------------------------------------------------------------------
    //! General destructor functor
    struct deleteFunctor
    {
        typedef void result_type;

        template <typename T>
        void operator( )( const T* p ) const
        {
            assert( p != NULL );
            delete p;
        }
    };

    //--------------------------------------------------------------------------
    //! Helper function returning argument
    //!
    //! Useful in conjunction with corlib::CollectQuantity to collect element
    //! pointers (or similar) itself (rather than a depending quantity).
    template< typename T >
    struct ReturnArg
    {
        typedef T result_type;
        T & operator()( T & t ) { return t; }
    };

    //--------------------------------------------------------------------------
    //! \brief Let compiler compute \f$ 2^N \f$ \tparam N The exponent
    template<unsigned N> 
    struct twoToTheN{ static const unsigned value = 2 * twoToTheN<N-1>::value; };
    //! \cond SKIPDOX
    template<> struct twoToTheN<0>{ static const unsigned value = 1; };
    //! \endcond

    //--------------------------------------------------------------------------
    //!\brief Let the compiler compute \f$M^N\f$ \tparam M base, \tparam exponent
    template<unsigned M, unsigned N>
    struct mToTheN{ static const unsigned value = M * mToTheN<M,N-1>::value; };
    //! \cond SKIPDOX
    template<unsigned M> struct mToTheN<M,0>{ static const unsigned value = 1; };
    //! \endcond

    //--------------------------------------------------------------------------
    //! \brief Let the compiler compute \f$ \sum_{i=1}^N i \f$ \tparam N
    template<unsigned N>
    struct sumUpTo{ static const unsigned value = N + sumUpTo<N-1>::value; };
    //! \cond SKIPDOX
    template<> struct sumUpTo<1>{ static const unsigned value = 1; };
    //! \endcond

    //--------------------------------------------------------------------------
    //! \brief Let the compiler compute \f$ N! \f$ \tparam N
    template<unsigned N>
    struct factorial{
        FTL_STATIC_ASSERT_MSG( (N < 13), "Risk of integer overflow" );   
        static const unsigned value = N * factorial<N-1>::value; };
    //! \cond SKIPDOX
    template<> struct factorial<0>{ static const unsigned value = 1; };
    //! \endcond

    //--------------------------------------------------------------------------
    //! \brief Let the compiler compute a number (double) to the power 'P'
    template<unsigned P>
    struct IntegralPower 
    {
        double operator()( const double x ) const {
            return IntegralPower<P-1>()(x) * x;
        }
    };
    //! \cond SKIPDOX
    template<>
    struct IntegralPower<0> { 
        double operator()( const double x ) const { return 1.; }
    };
    //! \endcond

    //--------------------------------------------------------------------------
    //! \brief Let the compiler compute the N-th root of a RADICANT
    template<unsigned RADICANT, unsigned N, int AUX=-1, bool less = false> 
    struct NthRoot 
    {
        // call with condition (AUX^N <= RADICANT)
        static const unsigned value= NthRoot<RADICANT,N,AUX-1,
                                             (mToTheN<AUX,N>::value<=RADICANT) >::value; 
    };
    
    //! \cond SKIPDOX
    template<unsigned RADICANT, unsigned N, int AUX> 
    struct NthRoot<RADICANT, N, AUX, true>
    {
        // here AUX^N <= RADICANT is true
        static const unsigned value = AUX+1; 
    };

    template<unsigned RADICANT, unsigned N> 
    struct NthRoot<RADICANT, N, -1, false> 
    {
        static const unsigned value = NthRoot<RADICANT, N, RADICANT, false>::value;
    }; 
    //! \endcond



    //--------------------------------------------------------------------------
    //! for_each with a predicate
    template<typename INP, typename FUN, typename PRED>
    FUN for_each_if( INP first, INP last, FUN f, PRED p )
    {
        for ( ; first != last; ++ first ) {
            if ( p( *first ) ) f( *first );
        }
        return f;
    }

    //--------------------------------------------------------------------------
    //! transform with predicate
    //!
    //! \param[in,out] first1    Input iterator to the initial position of the first sequence.
    //! \param[in]     last1     Input iterator to the final position of the first sequence.
    //! \param[in,out] result    Output iterator to the initial position of the range
    //!                          where function results are stored.
    //! \param[in]     trans     Unary function taking one element as argument,
    //!                          and returning some result value.
    //! \param[in]     pred      Unary predicate taking an element in the range as argument,
    //!                          and returning a value indicating the falsehood or truth
    //!                          of some condition applied to it.
    template<typename INP, typename OUT, typename FUN, typename PRED>
    OUT transform_if( INP first1, INP last1, OUT result, FUN trans, PRED pred )
    {
        for ( ; first1 != last1; ++first1 ) {
            if ( pred( *first1 ) )
                *result++ = trans( *first1 );
        }
        return result;
    }

    //--------------------------------------------------------------------------
    //! Return if a sequence [first,last) contains a provided value.
    //! 
    //! \param[in]     first     Input iterators to the initial position in a sequence.
    //! \param[in]     last      Input iterators to the final position in a sequence.
    //!                          The range used is [first,last), which contains all
    //!                          the elements betweenfirst and last, including the
    //!                          element pointed by first but not the element pointed by last.
    //! \param[in]     value     Value to be compared to each of the elements.
    //! \return                  True, if sequence [first,last) constains value,
    //!                          otherwise false.
    template< typename INP, typename VAL >
    bool contains( INP first, INP last, const VAL & value )
    {
        return ( std::find( first, last, value ) != last );
    }

    //--------------------------------------------------------------------------
    //! Return the position of a value with respect to lower
    //! limit of a sequence [first,last)
    //! 
    //! \param[in]     first     Input iterators to the initial position in a sequence.
    //! \param[in]     last      Input iterators to the final position in a sequence.
    //!                          The range used is [first,last), which contains all
    //!                          the elements betweenfirst and last, including the
    //!                          element pointed by first but not the element pointed by last.
    //! \param[in]     value     Value to be compared to each of the elements.
    //! \return                  Position of value in sequence [first,last). If it
    //!                          is not contained, the return value is the length of
    //!                          the sequence.
    template< typename INP, typename VAL >
    unsigned position( INP first, INP last, const VAL & value )
    {
        return std::distance( first, std::find( first, last, value ) );
    }

    //--------------------------------------------------------------------------
    //! Dummy predicate which always returns true
    struct Positive
    {
        template<typename T>
        bool operator()( const T ) { return true; }
    };

    //--------------------------------------------------------------------------
    //! Dummy predicate which always returns false
    struct Negative
    {
        template<typename T>
        bool operator()( const T ) { return false; }
    };

    //--------------------------------------------------------------------------
    //! Check if operand is not NULL
    struct ValidPointer
    {
        template<typename T>
        bool operator()( const T t)
        {
            return not(t == static_cast<T>(NULL) );
        }
    };

    //--------------------------------------------------------------------------
    //! Class template wrapper for std::min function template
    //!
    //! Useful for corlib::AccumulateQuantity
    template< typename T >
    struct Min : boost::function<const T&(const T&, const T&)>
    {
        const T & operator()( const T & x, const T & y ) const
        {
            return std::min( x, y );
        }
    };

    //--------------------------------------------------------------------------
    //! Class template wrapper for std::max function template
    //!
    //! Useful for corlib::AccumulateQuantity
    template< typename T >
    struct Max : boost::function<const T&(const T&, const T&)>
    {
        const T & operator()( const T & x, const T & y ) const
        {
            return std::max( x, y );
        }
    };

    //--------------------------------------------------------------------------
    //! Object for implicit singleton destruction
    template<typename SINGLETON>
    class SingletonDestroyer
    {
    public:
        SingletonDestroyer() : singleton_( NULL ) { }
        ~SingletonDestroyer() { delete singleton_; }

        void setSingleton( SINGLETON * s ) { singleton_ = s; }
        void releaseSingleton() { singleton_ = NULL; }
    private:
        SINGLETON * singleton_;
    };

    //--------------------------------------------------------------------------
    //! Convert seconds to a string of type hour:minutes:seconds
    std::string secondsToHMMSS( const double & seconds )
    {
        unsigned secondsInt = static_cast<unsigned>(seconds);
        std::stringstream pipe;
        const unsigned hours   = secondsInt/3600; secondsInt = secondsInt % 3600;
        const unsigned minutes = secondsInt/60;   secondsInt = secondsInt % 60;
        const double secDouble = seconds - 3600.*hours - 60.*minutes;
        pipe << hours << ":" 
             << std::setw(2) << std::setfill('0') << minutes << ":" 
             << ( secondsInt < 10 ? '0' : '\0' )  << secDouble << std::endl;
        return pipe.str();
    }

    //--------------------------------------------------------------------------
    //! Convert bool to unsigned to be used for accumulating based on bool functions
    unsigned bool2unsigned( bool state )
    {
        return state ? static_cast<unsigned>(1) : static_cast<unsigned>(0);
    }

    //--------------------------------------------------------------------------
    //! Return unsigned one
    unsigned getOne( )
    {
        return 1;
    }

    //--------------------------------------------------------------------------
    /** Connect given data points (t_1, f_1) and (t_2, f_2) by a fifth-order
     *  polynomial which has zero first- and second-order derivatives at the
     *  interval boundaries. This function evaluates the polynomial and its
     *  first derivative at any given point t.
     *
     * <pre>
     *     +++--------+--------+--------+---------+--------+--------+--------+++
     *     |             max(dot{f})  ---  #####                   f(t) ****** |
     *     +                             ###   ###            dot{f}(t) ###### +
     *     |                           ###       ###                           |
     *     +                          #             #                          +
     *     |                         ##             ##                         |
     *     +                        #                 #                        +
     *     |                       ##                 ##                       |
     *     +                       #                   #                       +
     *     |                      #                     #                      |
     * f_2 +---------------       #      ---------      #   ------  ************
     *     |                     #                      ************           |
     *     +                    #                    **** #                    +
     *     |                   #                  ***      #                   |
     *     +                  #                ***          #                  +
     *     |                  #              ***            #                  |
     *     +                ##            ***                ##                +
     *     |               #           ***                     #               |
     *     +              ##        ***                        ##              +
     *     |             #      *****                            #             |
     * f_1 *********************--------+---------+--------+------##############
     *               t_1                                           t_2     --> t   
     * </pre>
     * Note that \f$ \max_{t} \dot{f} = 30/16 = 1.875 \f$.
     *
     * \param[in]  t1   lower interval boundary
     * \param[in]  t2   upper interval boundary
     * \param[in]  f1   first function value
     * \param[in]  f2   second function value
     * \param[in]  t    evaluation coordinate
     * \param[out] f    function value of polynomial
     * \param[out] fDot value of first derivative 
     */
    void smoothConnect( const double t1, const double t2, 
                        const double f1, const double f2,
                        const double t,
                        double & f, double & fDot )
    {
        if ( t <= t1 ) {  // below lower boundary
            f    = f1;
            fDot = 0.;
        }
        else if ( t >= t2 ) { // above upper boundary
            f    = f2;
            fDot = 0.;
        }
        else { // evaluation in the interval
            const double df = f2 - f1;
            const double dt = t2 - t1;
            const double tr = (t - t1) / dt;
            f    = f1 + df * ( 6. * (tr*tr*tr*tr*tr) - 15. * (tr*tr*tr*tr) + 10. * (tr*tr*tr) );
            fDot = df / dt * (30. * (tr*tr*tr*tr)    - 60. * (tr*tr*tr)    + 30. * (tr*tr)    );
        }
    }
    

}// end namespace corlib
#endif
