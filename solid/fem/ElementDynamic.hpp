// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//
//                 Computational Structural Mechanics Lab
//                         University of Cambridge
//
// This software is copyrighted and all rights are retained by the University of
// Cambridge. All use, disclosure, and/or reproduction of any part not expressly
// authorised by the University of Cambridge is prohibited. (C) 2010
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

//! @file ElementDynamic.hpp

#ifndef solid_fem_elementdynamic_h
#define solid_fem_elementdynamic_h
//------------------------------------------------------------------------------
#include <solid/fem/ElementStatic.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
//------------------------------------------------------------------------------
namespace solid{
    namespace fem{
        template<typename BELEMENT, typename MAT> class ElementDynamic;

        namespace ublas = boost::numeric::ublas;
    }
}

//------------------------------------------------------------------------------
/** \brief   Extension of ElementStatic for dynamic problems
 * \details In addition to the element for solids, this element comutes and 
 * stores the mass matrix. 
 * \tparam BELEMENT Type of basis element
 * \tparam MAT      Type of material the element occupies
 */
//------------------------------------------------------------------------------
template<typename BELEMENT, typename MAT>
class solid::fem::ElementDynamic 
    : public solid::fem::ElementStatic<BELEMENT,MAT>
{
public:
    //! @name Convenience typedefs
    //@{
    typedef          BELEMENT                                    BasisElement;
    typedef typename BasisElement::Node                          Node;
    typedef typename solid::fem::ElementStatic<BasisElement,MAT> ElementStatic;
    typedef typename ElementStatic::VecLDim                      VecLDim;
    typedef typename ElementStatic::MatDofNN                     MatDofNN;
    //@}

    //! @name Attributes used here
    //@{
    static const unsigned     dof          = ElementStatic::dof;
    static const unsigned     numNodes     = ElementStatic::numNodes;
    //@}

public:

    //! Compute element mass matrix
    void massIntegrand( const VecLDim & xi, const double & weight,
                        ublas::matrix<double> & ) const;

    /** @name Accessors for nodal quantity */
    //@{
    //! Return nodal velocities
    void nodalVelocities( MatDofNN & V ) const;

    //! Return nodal accelerations
    void nodalAccelerations( MatDofNN & A ) const;
    //@}

    //! Compute nodal forces due to body force
    //! this is a temporary fix; else icc compiler 
    //! doesn't compile BodyForce
    template< typename FUNC > 
    void bodyForce( const VecLDim & xi, const double & weight, FUNC func,
                    const double & factor ){
        ElementStatic::bodyForce(xi, weight, func, factor);
    }

};
//------------------------------------------------------------------------------
#include "ElementDynamic.ipp"

//------------------------------------------------------------------------------
#endif
