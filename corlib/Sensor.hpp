// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//
//                 Computational Structural Mechanics Lab
//                         University of Cambridge
//
// This software is copyrighted and all rights are retained by the University of
// Cambridge. All use, disclosure, and/or reproduction of any part not expressly
// authorised by the University of Cambridge is prohibited. (C) 2010
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

//! @file Sensor.hpp

#ifndef corlib_sensor_h
#define corlib_sensor_h
//------------------------------------------------------------------------------
#include <utility>
#include <vector>
#include <algorithm>
#include <ostream>

#include <Eigen/Core>

//------------------------------------------------------------------------------
namespace corlib {
    template< typename OP >  class Sensor;

}

//------------------------------------------------------------------------------
/** \brief %Sensor equipped with a pointer to something and an operator for it
 *  \details This little object has access to an object which stores some 
 *  computational data, e.g., a node or an element. Using a given operator
 *  it can ask for this data and store it in a container. Such an operator
 *  could be the nodal displacements.
 *  \tparam OP   Type of operator acting on PTR
 */
template< typename OP >
class corlib::Sensor
{
public:
    typedef std::pair< double, double >    Datum_;
    typedef typename OP::DonatorType       ObjectPtr_;

    //! Constructor with a pointer to the object and the corr. functor
    Sensor( ObjectPtr_ object, OP op ) : objPtr_ ( object ), operator_(op) { }

    //! Empty destructor
    ~Sensor( ) { }

    //--------------------------------------------------------------------------
    /** Listen for a value corresponding to given abscissa and direction.
     *  Note that this works only for vector type quantities!!
     *  \param[in] abscissa The abscissa to the value it stores
     *  \param[in] comp     The component direction                  */
    void listenComponent( const double abscissa, const unsigned comp )
    {
        // expected return type from OP is Eigen vector
        const double value = (operator_( objPtr_ ))( comp );
        data_.push_back( Datum_( abscissa, value ) );
    }

    void listenValue( const double abscissa )
    {
        // expected return type from OP is double
        const double value = (operator_( objPtr_ ));
        data_.push_back( Datum_( abscissa, value ) );
    }

    //--------------------------------------------------------------------------
    /** Listen for a vectorial component corresponding to the given abscissa and
     *  store the 2-norm of it.
     *  \param[in] abscissa The abscissa corresponding to the 2-norm it stores */
    void listenNorm( const double abscissa )
    {
        // expected return type from OP is Eigen vector
        const double norm = (operator_( objPtr_) ).norm( );
        data_.push_back( Datum_( abscissa, norm ) );
    }

    //--------------------------------------------------------------------------
    /** Write the stored pairs of type <double,double> to a stream
     *  \param[in,out] out Output stream
     *  \retval        out Output stream     */
    std::ostream & write( std::ostream & out  ) const
    {
        // since there is no way to out<<std::pair, the classical way
        for ( unsigned i = 0; i < data_.size( ); i ++ )
            out << data_[i].first << "  " << data_[i].second << std :: endl;
        return out;
    }

private:
    ObjectPtr_               objPtr_;  //!< Pointer to some object (node,element..)
    OP                     operator_;  //!< Operator on that object 
    std::vector< Datum_ >      data_;  //!< Container with the data to store
};



#endif
