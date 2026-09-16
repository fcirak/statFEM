// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//
//                 Computational Structural Mechanics Lab
//                         University of Cambridge
//
// This software is copyrighted and all rights are retained by the University of
// Cambridge. All use, disclosure, and/or reproduction of any part not expressly
// authorised by the University of Cambridge is prohibited. (C) 2010
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

//! @file TensorQuadrature.hpp

#ifndef corlib_tensorquadrature_h
#define corlib_tensorquadrature_h
//------------------------------------------------------------------------------
//! system includes
#include <utility>
#include <array>
//! corlib includes
#include <corlib/misc.hpp>
#include <corlib/Shape.hpp>
#include <corlib/eigenX.hpp>

//------------------------------------------------------------------------------
namespace corlib{

    template<typename QUAD, unsigned DIM> class TensorQuadrature;

    namespace detail_{
        //! Overcome the fact, that member functions cannot be partially spec.
        template<typename QUAD,unsigned DIM> struct TensorQuadratureConstructor;
    }
}

//------------------------------------------------------------------------------
/** \brief   %Quadrature for tensor-product (0,1)^d elements (hypercubes)
 *  \details  Wrapper for any one-dimensional (0,1)-quadrature rule which 
 *   generates the points and weights for the corresponding d-fold tensor-product
 *  \param QUAD Underlying one-dimensional (0,1)-quadrature rule
 *  \param DIM  Dimension d of the d-fold tensor-product
 */
template<typename QUAD, unsigned DIM>
class corlib::TensorQuadrature
{
public:
    static const unsigned dim         = DIM;
    static const unsigned numPoints1D = QUAD::numPoints;
    static const unsigned numPoints   = corlib::mToTheN<numPoints1D,dim>::value;
    static const corlib::shape theShape = corlib::HyperCube<dim>::theShape;

    typedef eigenX::VectorSd<dim>                                                  VecDim;
    typedef std::pair<double, VecDim>                                    PairDoubleVecDim;
    typedef typename std::array<PairDoubleVecDim, numPoints>::const_iterator     QuadIter;

    //! Begin of array iterator
    QuadIter begin() const  { return weightsAndPoints_.begin(); }
    //! End of array iterator
    QuadIter end()   const  { return weightsAndPoints_.end();   }

    //! Constructor which generates the array of (weight,point)-pairs
    TensorQuadrature()
    {
        // Use constructor helper object defined below
        corlib::detail_::TensorQuadratureConstructor<QUAD,DIM> tqc;
        tqc( weightsAndPoints_ );
    }

private:
    //! Array of weight and point pairs for the tensor product rule
    std::array< PairDoubleVecDim, numPoints > weightsAndPoints_;
};

//! \cond SKIPDOX
//------------------------------------------------------------------------------
namespace corlib{
    namespace detail_{

        //------------------------------------------------------------------
        //! Constructor for one-dimensional 'tensor'-product rule 
        template<typename QUAD>
        struct TensorQuadratureConstructor<QUAD,1>
        {
            void operator()( std::array< std::pair<double, eigenX::VectorSd<1> >,
                                           QUAD::numPoints> & weightsAndPoints  )
            {
                // underlying one-dimensional quadrature rule
                QUAD quad1D;

                // fill array with weights and points from QUAD
                for ( unsigned i = 0; i < QUAD::numPoints; i ++ ) {
                    weightsAndPoints[i].first        = quad1D[i].first;
                    (weightsAndPoints[i].second)[0]  = quad1D[i].second;
                }

                return;
            }
        };
            
        //------------------------------------------------------------------            
        //! Constructor for general dimensions
        template<typename QUAD, unsigned DIM>
        struct TensorQuadratureConstructor
        {
            void operator()( std::array< std::pair<double, eigenX::VectorSd<DIM> >,
                             corlib::mToTheN<QUAD::numPoints,DIM>::value > & weightsAndPoints )
            {
                // lower-dimensional quadrature
                TensorQuadrature<QUAD,DIM-1> lowerQuad;
                const unsigned lowerNumPoints = TensorQuadrature<QUAD,DIM-1>::numPoints;
                    
                // copy lower-dimensional quadrature into arrays 
                for ( unsigned outer = 0; outer < QUAD::numPoints; outer ++ ) {
                    unsigned index = outer * lowerNumPoints;
                    for ( typename TensorQuadrature<QUAD,DIM-1>::QuadIter lower = lowerQuad.begin();
                          lower != lowerQuad.end(); ++ lower ) {
                            
                        // copy weights
                        weightsAndPoints[index].first = lower->first;
                        // copy points into [0,DIM-1) subrange of the DIM-points
                        ( weightsAndPoints[index].second ).head(DIM-1)
                            = ( lower->second ).head(DIM-1);
                        index ++;
                    }
                }

                // underlying one-dimensional quadrature rule
                QUAD quad1D;

                // make outer product for the new rule
                for ( unsigned outer = 0; outer < QUAD::numPoints; outer ++ ) {
                    for ( unsigned inner = 0; inner < lowerNumPoints; inner ++ ) {
                        const unsigned index = outer * lowerNumPoints + inner;
                        weightsAndPoints[index].first           *= quad1D[outer].first;
                        (weightsAndPoints[index].second)[DIM-1]  = quad1D[outer].second;

                    }
                }
            }
        };

    } 
}
//! \endcond

#endif
