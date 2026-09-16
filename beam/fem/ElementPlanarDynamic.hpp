// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//
//                 Computational Structural Mechanics Lab
//                         University of Cambridge
//
// This software is copyrighted and all rights are retained by the University of
// Cambridge. All use, disclosure, and/or reproduction of any part not expressly
// authorised by the University of Cambridge is prohibited. (C) 2010
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

//! @file ElementPlanarDynamic.hpp

#ifndef beam_fem_elementplanardynamic_h
#define beam_fem_elementplanardynamic_h
//------------------------------------------------------------------------------
#include <Eigen/Core>
#include <corlib/eigenX.hpp>
#include <beam/fem/ElementPlanarStatic.hpp>

//------------------------------------------------------------------------------
namespace beam{
    namespace fem{

        template< typename NODE, typename SFUN, typename MAT >
        class ElementPlanarDynamic;

        namespace eigenX = corlib::eigenX;
    }
}

//------------------------------------------------------------------------------
/**\ingroup elements
 * \brief   Extension of ElementPlanarStatic for dynamic problems
 *
 * \details In addition to the element for beams, this element computes and 
 * stores the mass matrix. 
 *
 * \tparam NODE  Type of node
 * \tparam SFUN  Shape function for the displacements
 * \tparam MAT   Type of material the element occupies
 */
//------------------------------------------------------------------------------
template< typename NODE, typename SFUN, typename MAT >
class beam::fem::ElementPlanarDynamic : 
    public beam::fem::ElementPlanarStatic< NODE, SFUN, MAT >
{
private:
    //! Type-access to static element
    typedef beam::fem::ElementPlanarStatic<NODE, SFUN, MAT> ElementPlanarStatic_;

public:
    //! @name Used inherited types 
    //@{
    typedef typename ElementPlanarStatic_::VecNF            VecNF;
    typedef typename ElementPlanarStatic_::VecLDim          VecLDim;
    typedef typename ElementPlanarStatic_::MatDofNF         MatDofNF; 
    //@}

    //! @name Used attributes
    //@{
    static const unsigned     dof               = NODE::dof;
    static const unsigned     numNodesSN        = ElementPlanarStatic_::numNodesSN;
    static const unsigned     matSize           = ElementPlanarStatic_::matSize;
    //@}

public:
    //--------------------------------------------------------------------------
    /** @name Set-up of element */
    //@{
    //! Constructor
    ElementPlanarDynamic() 
        : elemMass_( Eigen::MatrixXd::Zero( matSize, matSize ) )
    {
        return;
    }

    //--------------------------------------------------------------------------
    /** @name Mass matrix related things */
    //@{
    //! Clear the element stiffness matrix
    void clearMass()
    { 
        elemMass_.setZero();
        return;
    }

    //! Compute element mass matrix
    void massIntegrand( const VecLDim& xi, const double& weight );

    //! Add element mass matrix
    void addMassMatrix( Eigen::MatrixXd & matrix ) const
    { 
        matrix += elemMass_;
    }
    //@}

public:
    //--------------------------------------------------------------------------
    /** @name Accessors for nodal quantity */
    //@{
    //! Return nodal velocities
    void nodalVelocities( MatDofNF & vel ) const;

    //! Return nodal accelerations
    void nodalAccelerations( MatDofNF & acc ) const;
    //@}

    //! Interpolate the nodal velocity at local coordinate
    typename NODE::VecDof giveVelocity( const VecLDim& xi ) const;

protected:
    using ElementPlanarStatic_::material_;
    using ElementPlanarStatic_::area_;

    //! Element mass matrix
    typename ElementPlanarStatic_::ElemMat_      elemMass_;
};
//------------------------------------------------------------------------------
#include "ElementPlanarDynamic.ipp"

//------------------------------------------------------------------------------
#endif
