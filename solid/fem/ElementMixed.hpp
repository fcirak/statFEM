// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//
//                 Computational Structural Mechanics Lab
//                         University of Cambridge
//
// This software is copyrighted and all rights are retained by the University of
// Cambridge. All use, disclosure, and/or reproduction of any part not expressly
// authorised by the University of Cambridge is prohibited. (C) 2010
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

//! @file ElementMixed.hpp

#ifndef solid_fem_elementmixed_h
#define solid_fem_elementmixed_h
//------------------------------------------------------------------------------
#include <solid/fem/ElementStatic.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>

namespace solid{
    namespace fem{
        
        template< typename NODE, typename SFUN, typename SFUNP, typename MAT>
        class ElementMixed;

        namespace ublas = boost::numeric::ublas;
    }
}

//------------------------------------------------------------------------------
/**\ingroup elements
 * \brief Element for a mixed (u,p)-formulation
 * \details This object extends the ElementStatic by using additional pressure
 *  degrees of freedom. Currently, this element is hard-coded as a Taylor-Hood
 *  element. It overloads the stiffness kernel of ElementStatic and provides the
 *  functionality to compute the mixed form b(u,p) and some stabilization form.
 *  \tparam NODE   Type of nodes (should have a pressure degree of freedom)
 *  \tparam SFUNU  Type of shape function for the displacements
 *  \tparam SFUNP  Type of shape function for the pressure variables
 *  \tparam MAT    Type of material the element occupies
 */
//------------------------------------------------------------------------------
template< typename NODE, typename SFUN, typename SFUNP, typename MAT>
class solid::fem::ElementMixed : 
    public solid::fem::ElementStatic<NODE, SFUN, MAT>
{
private:
    typedef solid::fem::ElementStatic<NODE,SFUN,MAT> ElementStatic_;

public:

    static const unsigned    dim           = ElementStatic_::dim;
    static const unsigned    dof           = ElementStatic_::dof;
    static const unsigned    localDim      = ElementStatic_::localDim;
    static const unsigned    numNodes      = SFUN::numFunctions;
    static const unsigned    numPressNodes = SFUNP::numFunctions;
    static const unsigned    matSize       = numNodes * dof;
    static const unsigned    matSizeP      = numPressNodes;
    
    typedef typename ElementStatic_::MatDimNN                        MatDimNN;
    typedef typename ElementStatic_::MatLDimNN                       MatLDimNN;
    typedef typename ElementStatic_::VecDim                          VecDim;
    typedef typename ElementStatic_::VecLDim                         VecLDim;
    typedef typename ElementStatic_::Mat3x3                          Mat3x3;

private:
    typedef ublas::bounded_vector<double, numPressNodes>            VecNPN_;
    typedef ublas::bounded_matrix<double, matSize,  matSizeP>       MixedMat_;
    typedef ublas::bounded_matrix<double, matSizeP, matSizeP>       StabilMat_;
    typedef corlib::ElementBasic<NODE,numNodes,SFUN::myShape>       ElementBasic_;
    
public:
    //! Constructor sets pressure flag to true for the first nodes
    void markPressureNodes( ) 
    {
        for ( unsigned v = 0; v < numPressNodes; v ++ ) 
            ElementBasic_::nodes_[ v ] -> setPressureFlag( );
    }

protected:
    //! Evaluate the pressure shape function
    void sfunPress_( const VecLDim & xi, VecNPN_ & psi ) const
    {
        sfunp_.evaluate( xi, psi );
    }

public:
    //! Give pressure dof indices
    void getDofIndicesP( std::vector<unsigned> & dofIndices ) const
    {
        std::vector<unsigned>::iterator iter = dofIndices.begin();

        for ( unsigned v = 0; v < numPressNodes; v ++ ) { 
            std::vector<unsigned> aux( dof );
            ElementBasic_::nodes_[v] -> copyDofArrayP( aux );
            iter = std::copy( aux.begin(), aux.end(), iter );
        }
        return;
    }

    //! Clear the entries of the mixed matrix
    void clearMixedMat( ){  mixedMat_.clear(); return; }
    
    //! Clear the entries of the stabilization matrix
    void clearStabilMat(){ stabilMat_.clear(); return; }

    //! Interpolate the pressure values to non-pressure nodes
    void interpolatePressure();

    //! Kernel function of a(u,v)
    void stiffnessIntegrand( const VecLDim & xi, const double & weight );

    //! Kernel function of b(u,p)
    void mixedFormIntegrand( const VecLDim & xi, const double & weight );

    //! So far, \f$ \int p q dX\f$  which normalizes the pressure
    void meanPressIntegrand( const VecLDim & xi, const double & weight );

    //! Return the Mixed form b(u,p)
    ublas::matrix<double> giveGradientMatrix()   const { return mixedMat_; }
    ublas::matrix<double> giveDivergenceMatrix() const { return trans(mixedMat_); }

    //! Return the stabilization form c(p,q)
    ublas::matrix<double> giveStabilMatrix(   ) const { return stabilMat_; }


protected:
    using ElementStatic_::material_;
    using ElementStatic_::elemStiff_;

    SFUNP         sfunp_;      //!< Shape function for the pressure
    MixedMat_     mixedMat_;   //!< Element matrix of b(u,p)
    StabilMat_    stabilMat_;  //!< Element matrix of c(p,q)
};

//------------------------------------------------------------------------------
#include "ElementMixed.ipp"
//------------------------------------------------------------------------------
#endif
