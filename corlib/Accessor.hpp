// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//
//                 Computational Structural Mechanics Lab
//                         University of Cambridge
//
// This software is copyrighted and all rights are retained by the University of
// Cambridge. All use, disclosure, and/or reproduction of any part not expressly
// authorised by the University of Cambridge is prohibited. (C) 2010
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

//! @file Accessor.hpp

#ifndef corlib_accessor_h
#define corlib_accessor_h

#include <boost/function.hpp>
#include <corlib/eigenX.hpp>

//------------------------------------------------------------------------------
namespace corlib {

    template<typename FUN> class Accessor;

    namespace detail_{
        template<typename RET> class AccessorTraits;
    }

    //--------------------------------------------------------------------------
    /** Helper function to directly get the Accessor object by means 
     *  of a function pointer --- sufficient for most purposes    
     *  \tparam RES  result type of the function object   (type of datum)
     *  \tparam TP   argument type of the function object (type of donator)
     */
    template<typename RES, typename TP>
    corlib::Accessor< boost::function< RES( const TP* )> >
    accessorFun( RES (TP::*f)() const, std::string name = "" ) {
        return corlib::Accessor<boost::function< RES( const TP * )> >( f, name );
    }

}

//------------------------------------------------------------------------------
/** \brief Generic accessor object which retrieves a datum from an object
 *  \details Given a function object 'f' with signature 
 *                   result_type f( const argument_type * )
 *   this object queries the datum of type 'result_type' from an object of type
 *   'argument_type'. Note that argument_type and result_type have to be 
 *   provided by the function object. This will be the case for std- (and boost-)
 *   function objects. Moreover, in most applications it will be sufficient to 
 *   use the above defined helper function 'accessorFun'.
 *  \tparam FUN  Type of function object 
 */
template<typename FUN>
class corlib::Accessor : public FUN
{
public:
    typedef typename FUN::result_type   ReturnType;
    typedef typename FUN::argument_type DonatorType;
    typedef typename detail_::AccessorTraits<ReturnType>::ValueType ValueType;

    //! Provide information of the number of components
    static const unsigned numComponents = detail_::AccessorTraits<ReturnType>::size;
    //! Provide a zero element for accumulation
    static const ReturnType zero() { return detail_::AccessorTraits<ReturnType>::zero(); }

    //! Constructor with object member function and name of the quantity
    Accessor( FUN func, std::string name = "" )
        : func_( func ), name_( name )
    { }

    //! The function call by accessing the objets member function
    ReturnType operator()( DonatorType objPtr ) const
    {
        return func_( objPtr );
    }

    //! Return the name of the quantity
    std::string name() const { return name_; }

private:
    FUN func_;                              //!< Access function object
    const std::string name_;                //!< Quantity's name
};

//------------------------------------------------------------------------------
namespace corlib {

    namespace detail_{

        //--------------------------------------------------------------------------
        /** \brief Traits class for the size of the ReturnValue and its zero element
         *  \details This class provides some necessary data to the Accessor object.
         *  It provides the information of the number of components the return_type
         *  has and a zero element for accumulation functions.
         *  Note that, by default, the type RET is assumed to be a 
         *  eigenX::VectorSd<RowsAtCompileTime> which provides the max_size
         *  const-expression and the clear() function.
         *  For other data types as return values, specializations of this traits
         *  class have to be provided.
         *  \tparam RET  The return value of the Accessor object
         */
        template<typename RET>
        class AccessorTraits
        {
        public:
            typedef double ValueType;
            static const unsigned size = RET::RowsAtCompileTime;
            static RET zero() { RET aux; aux.setZero(); return aux; }
        };

        // \cond SKIPDOX
        // specialization: RET=eigenX::MatrixSd<3>
        template<> class AccessorTraits<eigenX::MatrixSd<3> >
        {
        public:
            typedef eigenX::MatrixSd<3> ReturnType;
            typedef double ValueType;
            static const unsigned size = 3*3;
            static ReturnType zero() { ReturnType aux; aux.setZero(); return aux; }
        };

        // specialization: RET=eigenX::MatrixSd<3,2>
        template<> class AccessorTraits<eigenX::MatrixSd<3,2> >
        {
        public:
            typedef eigenX::MatrixSd<3,2> ReturnType;
            typedef double ValueType;
            static const unsigned size = 3*2;
            static ReturnType zero() { ReturnType aux; aux.setZero(); return aux; }
        };

        // specialization: RET=double
        template<> class AccessorTraits<double>
        {
        public:
            typedef double ValueType;
            static const unsigned size = 1;
            static double zero() { return 0.; }
        };

        // specialization: RET=unsigned
        template<> class AccessorTraits<unsigned>
        {
        public:
            typedef unsigned ValueType;
            static const unsigned size = 1;
            static unsigned zero() { return 0; }
        };

        // specialization: RET=int
        template<> class AccessorTraits<int>
        {
        public:
            typedef int ValueType;
            static const unsigned size = 1;
            static  int zero() { return 0; }
        };

        // spezialization: RET=bool
        template<> class AccessorTraits<bool>
        {
        public:
            typedef bool ValueType;
            static const unsigned size = 1;
            static bool zero() { return 0; }
        };
        // \endcond

    }
}

#endif
