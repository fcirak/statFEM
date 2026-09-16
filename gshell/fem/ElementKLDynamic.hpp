//
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//
//                 Computational Structural Mechanics Lab
//                         University of Cambridge
//
// This software is copyrighted and all rights are retained by the CSMLab. All
// use, disclosure, and/or reproduction of any part not expressly authorized by
// F Cirak is prohibited. (C) 2010.
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//

//! @author Burkhard Bornemann, Kosala Bandara
//! @date   2011

#ifndef gshell_fem_elementkldynamic_h
#define gshell_fem_elementkldynamic_h

#include <Eigen/Core>

#include <gshell/fem/ElementKLStatic.hpp>

//------------------------------------------------------------------------------
namespace gshell {
    namespace fem {
        
        template< typename BELEMENT, typename MAT >
        class ElementKLDynamic;

    }
}

//==============================================================================
template< typename BELEMENT, typename MAT >
class gshell::fem::ElementKLDynamic : 
    public gshell::fem::ElementKLStatic< BELEMENT, MAT >
{
private:
    typedef gshell::fem::ElementKLStatic< BELEMENT, MAT >          ElementKLStatic_;

public:
    typedef BELEMENT                                               BasisElement;
    typedef typename BasisElement::ShapeFun                        ShapeFun;
    typedef typename BasisElement::Node                            Node;
    typedef MAT                                                    Material;

    //! Simplex shape of element -- can be corlib::TRIANGLE or corlib::QUADRILATERAL
    static const corlib::shape myShape   = ElementKLStatic_::myShape;
    //! Number of vertices of simplex shape -- can be 3 or 4
    static const unsigned      numVertices = ElementKLStatic_::numVertices;
    //! Dimension of embedding space -- always 3
    static const unsigned      dim       = ElementKLStatic_::dim;  // always 3
    //! Total number of DOFs at nodes -- always 3
    static const unsigned      dof       = ElementKLStatic_::dof;  // always 3
    //! Dimension of shell manifold -- always 2
    static const unsigned      localDim  = ElementKLStatic_::localDim; // always 2

    typedef typename ElementKLStatic_::VecNF                       VecNF;
    typedef typename ElementKLStatic_::VecDim                      VecDim;
    typedef typename ElementKLStatic_::VecDof                      VecDof;
    typedef typename ElementKLStatic_::VecLDim                     VecLDim;
    typedef typename ElementKLStatic_::MatDimNF                    MatDimNF;
    typedef typename ElementKLStatic_::MatLDimNF                   MatLDimNF;
    typedef typename ElementKLStatic_::MatDofNF                    MatDofNF;

public:
    //! Constructor
    ElementKLDynamic() : ElementKLStatic_() { return; }

    //! @name Repeated functions to help binding
    //!
    //! <b>NOTE:</b> Intel compiler has problem to detect NodeStatic_'s function when binding
    //@{
    void getDofIndices( std::vector< unsigned > & dofIndices ) const
    {
        this->ElementKLStatic_::getDofIndices( dofIndices );
    }
    
    unsigned numFunctions( ) const { return ElementKLStatic_::numFunctions( ); }

    Node * giveSupportNodePtr( const unsigned ind ) const
    {
        return this->ElementKLStatic_::giveSupportNodePtr( ind );
    }

    Eigen::MatrixXd giveStiffnessMatrix( ) const
    { 
        return this->ElementKLStatic_::giveStiffnessMatrix( );
    }

    void nodalIncrements( MatDimNF & deltaU ) const
    {
        this->ElementKLStatic_::nodalIncrements( deltaU ); return;
    }
    //@}

    //! Compute element mass matrix
    void massIntegrand( const VecLDim & xi, const double & weight,
                        Eigen::MatrixXd & elemMass ) const;

    //! Compute lumped element mass matrix
    void lumpedMassIntegrand( const VecLDim & xi, const double & weight,
                              Eigen::MatrixXd & lumpedMass ) const;

public:
    //! @name Accessors for nodal quantities
    //@{
    //! Return nodal velocities
    //!
    //! \param[out]   vel   Matrix with the nodal mid-surface velocities as columns
    //!                     \f$[v^A{}_K]\f$
    void nodalVelocities( MatDimNF & vel ) const;

    //! Return nodal accelerations
    //!
    //! \param[out]   acc   Matrix with the nodal mid-surface accelerations as columns
    //!                     \f$[a^A{}_K]\f$
    void nodalAccelerations( MatDimNF & acc ) const;
    //@}

protected:
    using ElementKLStatic_::material_;
    using ElementKLStatic_::thickness_;

};

//------------------------------------------------------------------------------
#include "ElementKLDynamic.ipp"

#endif
