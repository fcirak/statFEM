// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//
//                 Computational Structural Mechanics Lab
//                         University of Cambridge
//
// This software is copyrighted and all rights are retained by the University of
// Cambridge. All use, disclosure, and/or reproduction of any part not expressly
// authorised by the University of Cambridge is prohibited. (C) 2010
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

//! @file Quadrature.hpp

#ifndef corlib_quadrature_h
#define corlib_quadrature_h
//------------------------------------------------------------------------------
#include <array>
#include <utility>

#include <corlib/Shape.hpp>
#include <corlib/TensorQuadrature.hpp>
#include <corlib/GaussLegendre.hpp>
#include <corlib/SimplexQuadrature.hpp>
#include <corlib/eigenX.hpp>

//------------------------------------------------------------------------------
namespace corlib {
    template< corlib::shape SHAPE, unsigned QPS > class Quadrature;
    
    namespace detail_{
        //! Helper to generate a quadrature object according to element type
        template<bool isSimplex, unsigned DIM, unsigned QPS> 
        struct QuadratureConstructor;
    }
}

//------------------------------------------------------------------------------
/** \brief Container for numerical integration rules
 *  \details Storage of an array of pairs of weights and points
 *  \tparam SHAPE the shape of the element (LINE, TRIANGLE, etc.)
 *  \tparam QPS   the number of quadrature points
 */
//------------------------------------------------------------------------------
template< corlib::shape SHAPE, unsigned QPS>
class corlib::Quadrature 
{
public:

    static const corlib::shape theShape   = SHAPE;     
    static const unsigned      dim        = corlib::ShapeTraits<SHAPE>::dim; 
    static const bool          isASimplex = (corlib::Simplex<dim>::theShape == theShape);
    static const unsigned      numPoints  = detail_::QuadratureConstructor<isASimplex,dim,QPS>::numPoints;
    static unsigned size(  ) { return numPoints; }     
    
    typedef eigenX::VectorSd< dim >                  VecDim;
    typedef std::pair< double, VecDim >              WeightedPoint;
    typedef std::array< WeightedPoint, numPoints >   QuadratureArray;
    typedef typename QuadratureArray::const_iterator QuadIter;

public:
    //! Constructor using the constructor helper
    Quadrature( )
    {
        detail_::QuadratureConstructor<isASimplex, dim, QPS>()( weightsAndPoints_ );
    }

    //! Empty destructor 
    ~Quadrature( ) { }

    //! Return iterator pointing to the begin of the weighted points
    QuadIter begin( ) const { return weightsAndPoints_.begin( ); }
    //! Return iterator pointing to the end of the weighted points
    QuadIter end(   ) const { return weightsAndPoints_.end( ); }

protected:
    QuadratureArray weightsAndPoints_;  //! container with the weights and points
    
};

//! \cond SKIPDOX
namespace corlib{
    namespace detail_{
        
        // Integration over a hyper cube (n>=1)
        template<unsigned DIM,unsigned QPS>
        struct QuadratureConstructor<false,DIM,QPS>
        {
            // compute the integral part of the DIM-th root of QPS
            static const unsigned numPoints1D = corlib::NthRoot<        QPS,DIM>::value;
            // total number of points is the previous result to the DIM-th power
            static const unsigned numPoints   = corlib::mToTheN<numPoints1D,DIM>::value;

            void operator()( std::array< std::pair<double, eigenX::VectorSd<DIM> >, numPoints > &
                             weightsAndPoints )
            {
                corlib::detail_::TensorQuadratureConstructor< corlib::GaussLegendre<numPoints1D>, DIM >()( weightsAndPoints );
            }
            
        };

     
        // Integration over a n-Simplex (n>1)
        template<unsigned DIM, unsigned QPS>
        struct QuadratureConstructor<true,DIM,QPS>
        {
            static const unsigned numPoints = QPS;

            void operator()( std::array< std::pair<double, eigenX::VectorSd<DIM> >, QPS > &
                             weightsAndPoints )
            {
                corlib::detail_::SimplexQuadratureConstructor<DIM,QPS>()( weightsAndPoints );
            }
            
        };

    }
}
//! \endcond

//------------------------------------------------------------------------------
#endif
