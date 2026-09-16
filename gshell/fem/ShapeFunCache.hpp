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
//! @date   2010

#ifndef gshell_fem_shapefuncache_h
#define gshell_fem_shapefuncache_h

#include <array>
#include <tuple>
#include <unordered_map>

#include <corlib/fuzzyEqual.hpp>
#include <corlib/Shape.hpp>

#include <Eigen/Core>
#include <corlib/eigenX.hpp>

//------------------------------------------------------------------------------
namespace gshell {
    namespace fem {

        class VecLDimEqual;

        class VecLDimHash;

        template< corlib::shape SHAPE >
        class ShapeFunCache;

        namespace eigenX = corlib::eigenX;

    }
}

//------------------------------------------------------------------------------
/// Equality comparator of two surf vectors
///
/// Comparator to determine if two "local" surf vectors
/// are (approximately) equal
class gshell::fem::VecLDimEqual : 
    std::binary_function< eigenX::VectorSd<2>,
                          eigenX::VectorSd<2>,
                          bool >
{
public:
    typedef eigenX::VectorSd<2> VecLDim;

    bool operator()( const VecLDim & a, const VecLDim & b ) const
    {
        return ( corlib::fuzzyEqual( a[ 0 ], b[ 0 ] ) and 
                 corlib::fuzzyEqual( a[ 1 ], b[ 1 ] ) );
    }
};

//------------------------------------------------------------------------------
/// Hasher for eigenX::VectorSd<2> vectors
class gshell::fem::VecLDimHash :
    std::unary_function< eigenX::VectorSd<2>, std::size_t >
{
public:
    typedef eigenX::VectorSd<2> VecLDim;

    std::size_t operator()( VecLDim v ) const
    {
        return sizeof( VecLDim );
    }
};

//------------------------------------------------------------------------------
/// Interface to subdivision shape functions and their parametric derivatives
///
/// This class provides an interface to subdivision shape functions. The
/// shape functions and their derivatives can be evaluated at parametric
/// co-ordinate \$\xi\f$.
///
/// It stores (if requested) the shape functions to prevent repeated computation
/// using (potentially intensive) subdivison.
template< corlib::shape SHAPE >
class gshell::fem::ShapeFunCache
{
public:
    /// Simplex shape of element -- can be corlib::TRIANGLE or corlib::QUADRILATERAL
    static const corlib::shape myShape   = SHAPE;
    /// Number of vertices of simplex shape -- can be 3 or 4
    static const unsigned      numVertices = corlib::ShapeTraits<myShape>::numVertices;
    /// Dimension of shell manifold -- always 2
    static const unsigned      localDim  = corlib::ShapeTraits<myShape>::dim;
    /// Number of effective second derivatives subtracting duplicities due to symmetry -- always 3
    static const unsigned      sDim      = (localDim+1)*localDim/2;
    /// Map to get Voigt indices \f$\mathcal{A}\f$ from pair of local indices \f$(\alpha,\beta)
    static const unsigned      voigtForward[localDim][localDim];// = { { 0, 1 }, { 1, 2 } };

    typedef eigenX::VectorSd<localDim>              VecLDim;

    /// Vector for holding shape functions at evaluation point
    typedef Eigen::VectorXd                         VecNF;
    /// Matrix for holding first derivatives shape functions at evaluation point
    typedef Eigen::MatrixXd                         MatLDimNF;
    /// Matrix for holding second derivatives shape functions at evaluation point.
    /// Second derivatives are stored in Voigt-pattern
    /// \f$[ \phi_{,1,1}^K \phi_{,1,2}^K \phi_{,2,2}^K ]^T\f$
    typedef Eigen::MatrixXd                         MatSDimNF;

protected:
    /// Container for storing evaluated shape functions
    /// and their first and second derivatives
    typedef std::tuple< VecNF, MatLDimNF, MatSDimNF >              ShapeFunGradHess_;
    /// Map holding lcoal/parametric point and shape functions and their
    /// derivatives evaluated there
    typedef std::unordered_map< VecLDim, ShapeFunGradHess_,
                                  VecLDimHash, VecLDimEqual >      MapLCoordToShapeFun_;
    /// Const iterator on #MapLCoordToShapeFun_
    typedef typename MapLCoordToShapeFun_::const_iterator          MapLCoordToShapeFunConstIter_;

public:
    /// Empty constructor
    ShapeFunCache() { }

    /// Number of shape functions
    unsigned numFunctions() const
    {
        unsigned numFunctions = 0;
        if ( not evaluatedFun_.empty() ) {
            numFunctions = std::get< 0 >( evaluatedFun_.begin()->second ).size();
        }
        return numFunctions;
    }

    //! Populate evaluated shape functions and their derivatives
    template< typename SFUN, typename QUADRATURE >
    void populateEvaluatedShapeFunctions( SFUN * shapeFun,
                                          const QUADRATURE & quadrature );

    //! Populate limit surface coefficients
    template< typename SFUN >
    void populateLimitCoefficientsAtVertices( SFUN * shapeFun );

    //! Populate tangent coefficients at vertices
    template< typename SFUN >
    void populateTangentCoefficientsAtVertices( SFUN * shapeFun );

public:
    /// @name Evaluation of shape function and gradient
    //@{

    /// Return values of shape function at xi
    /// \param[in]  xi   Local cooordinate at which the functions is evaluated
    /// \param[out] phi  Result of all shape functions at xi 
    void evaluate( const VecLDim & xi, VecNF & phi ) const;

    /// Return gradient of shape function at xi
    //!
    //! \param[in]  xi       Local coordinate at which the gradient is evaluated
    //!                      \f$\xi^\alpha\f$
    //! \param[out] dPhiDXi  First derivative of shape functions
    //!                      \f$ [\varphi^K{}_{,\alpha}]^T = [\partial\varphi^K / \partial\xi^\alpha]^T \f$
    void evaluateGradient( const VecLDim & xi, MatLDimNF & dPhiDXi ) const;

    /// Return 2nd derivative of shape function at xi
    //!
    //! \param[in]  xi         Local coordinate at which the gradient is evaluated
    //!                        \f$\xi^\alpha\f$
    //! \param[out] ddPhiDDXi  2nd derivative of shape functions
    //!                        \f$\varphi^K{}_{,\alpha\beta}\f$
    //!                        \f$=\partial^2\varphi^K / \partial\xi^\alpha\partial\xi^\beta\f$.
    //!                        The result is stored
    //!                        \f$[\varphi^K{}_{,11}, \varphi^K{}_{,12}, \varphi^K{}_{,22}]^T\f$.
    void evaluateHessian( const VecLDim & xi, MatSDimNF & ddPhiDDXi ) const;

    /// Return value, 1st and 2nd derivatives of shape functions at xi
    //!
    //! \param[in]  xi         Local coordinate at which the gradient is evaluated
    //!                        \f$\xi^\alpha\f$
    //! \param[out] phi        Result of all shape functions \f$[\phi^K]\f$ at xi 
    //! \param[out] dPhiDXi    First derivative of shape functions
    //!                        \f$[\partial\varphi^K / \partial\xi^\alpha]^T\f$
    //! \param[out] ddPhiDDXi  2nd derivative of shape functions
    //!                        \f$\varphi^K{}_{,\alpha\beta}\f$
    //!                        \f$=\partial^2\varphi^K / \partial\xi^\alpha\partial\xi^\beta\f$.
    //!                        The result is stored
    //!                        \f$[\varphi^K{}_{,11}, \varphi^K{}_{,12}, \varphi^K{}_{,22}]^T\f$.
    void evaluateGradHess( const VecLDim & xi, VecNF & phi,
                           MatLDimNF & dPhiDXi, MatSDimNF & ddPhiDDXi ) const;

    /// Return limit position
    ///
    /// \param[in]     vIndex    Local facet-wise index of node
    /// \param[out]    coeff     Coefficients of limit position w.r.t. vertices
    void evaluate( const unsigned vIndex, VecNF & coeff ) const
    {
        coeff = limitCoeffs_[ vIndex ];
        return;
    }

    /// Return tangents
    ///
    /// \param[in]     vIndex    Local facet-wise index of node
    /// \param[out]    coeff     Coefficients of tangents w.r.t. vertices
    void evaluateGradient( const unsigned vIndex, MatLDimNF & coeff ) const
    {
        coeff = tangCoeffs_[ vIndex ];
        return;
    }

    //@}

private:
    /// Map from local co-ordinates (e.g. Gauss points) to
    /// evaluated shape functions and their derivatives
    MapLCoordToShapeFun_     evaluatedFun_;

    /// Limit coefficients at vertices
    std::array< VecNF, numVertices >      limitCoeffs_;

    /// Tangent coefficients at vertices
    std::array< MatLDimNF, numVertices >  tangCoeffs_;
};

//------------------------------------------------------------------------------
#include "ShapeFunCache.ipp"

#endif
