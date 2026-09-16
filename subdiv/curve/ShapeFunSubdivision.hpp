//
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//
//                 Computational Structural Mechanics Lab
//                         University of Cambridge
//
// This software is copyrighted and all rights are retained by the CSMLab. All
// use, disclosure, and/or reproduction of any part not expressly authorized by
// F Cirak is prohibited. (C) 2011.
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//

//! @author Kosala Bandara, Burkhard Bornemann
//! @date   2011

#ifndef subdiv_curve_shapefunsubdivision_h
#define subdiv_curve_shapefunsubdivision_h

#include <iostream>

#include <corlib/misc.hpp>
#include <corlib/Shape.hpp>

#include <subdiv/curve/collect.hpp>
#include <subdiv/curve/ShapeDim.hpp>
#include <subdiv/curve/Mesh.hpp>

#include <Eigen/Core>
#include <corlib/eigenX.hpp>

namespace subdiv{
    namespace curve{

        template< typename EDGE, typename SUBDIV >
        class ShapeFunSubdivision;

    }
}

//==============================================================================
/// Shape functions using subdivision
///
/// \tparam     EDGE    Edge type
/// \tparam    SUBDIV   Subdivision scheme type
template< typename EDGE, typename SUBDIV >
class subdiv::curve::ShapeFunSubdivision
{
public:
    typedef EDGE                                              Edge;
    typedef SUBDIV                                            Subdiv;

    static const corlib::shape myShape    = Edge::myShape;
    static const unsigned numVertices     = Edge::numVertices;
    static const unsigned numNeighbors    = Edge::numNeighbors;

    /// number of shape functions over regular patch
    static const unsigned numFunctions    = Subdiv::numFunctions;
    /// manifold dimension
    static const unsigned localDim        = corlib::ShapeTraits< myShape >::dim;
    /// Number of effective second derivatives subtracting duplicities due to symmetry
    static const unsigned sDim            = localDim*(localDim+1)/2;

    /// Real vector
    typedef eigenX::VectorSd< localDim >                   VecLDim;
    /// Real square matrix
    typedef eigenX::MatrixSd< localDim, localDim >         MatLDimLDim;
    /// number of shape functions of regular patch
    typedef eigenX::VectorSd< numFunctions >               VecNF;
    ///
    typedef eigenX::MatrixSd< localDim, numFunctions >     MatLDimNF;
    ///
    typedef eigenX::MatrixSd< sDim, numFunctions >         MatSDimNF;
    ///
    typedef Eigen::Matrix< VecNF, localDim, localDim >     MatVecNFLDimLDim;
    /// General real vector
    typedef Eigen::VectorXd                                Vec;
    /// General real matrix
    typedef Eigen::MatrixXd                                Mat;

    /// 
    typedef typename Edge::Vertex::VecDim                  VecDim;

    typedef std::vector< Edge * >                          VecFPtr;
    typedef typename VecFPtr::iterator                     VecFPtrIter;

    typedef typename Edge::Vertex                          Vertex;
    typedef std::set< Vertex * >                           SetVPtr;
    typedef typename SetVPtr::iterator                     SetVPtrIter;

private:
    /// Type of spline over regular patch
    typedef typename Subdiv::Spline                        Spline_;

    /// Edge in patch mesh
    typedef Edge< Vertex, Edge >                           Edge_;
    /// Edge tree in patch mesh
    typedef EdgeTree< myShape, Vertex >                    EdgeTree_;
    /// Mesh type of patch mesh
    typedef Mesh< Vertex, Edge_, EdgeTree_ >               Mesh_;

    typedef typename Vertex::VecP                          VecP_;
    typedef std::array< Vertex *, numFunctions >           VecVPtrNF_;

public:
    /// Constructor
    ///
    /// \param[in]   f        Edge on which subdivision shape functions are computed
    /// \param[in]   subdiv   Subdivision scheme
    ShapeFunSubdivision( Edge * f, Subdiv * subdiv );

    /// Return number of patch vertices
    unsigned numPatchVertices( ) const { return patchVertices_.size( ); }

    /// Return iterator to begin of patch vertices container
    SetVPtrIter patchVerticesBegin( ) const { return patchVertices_.begin( ); }
    /// Return iterator to end of patch vertices container
    SetVPtrIter patchVerticesEnd( ) const { return patchVertices_.end( ); }

    //@name Evaluate shape functions and derivatives at given point
    //@{

    /// Evaluate shape functions \f$\phi^k\f$ at point \f$\xi\f$
    void evaluate( const VecLDim & xi,
                   Vec & phi );
    /// Evaluate 1st derivative of shape functions \f$\phi_{,\xi}^k\f$ at point \f$\xi\f$
    void evaluateGradient( const VecLDim & xi,
                           Mat & dPhiDXi );
    /// Evaluate 2nd derivative of shape functions \f$\phi_{,\xi,\xi}^k\f$ at point \f$\xi\f$
    void evaluateHessian( const VecLDim & xi,
                          Mat & ddPhiDDXi );

    /// At point \f$\xi\f$ evaluate shape functions \f$\phi^k\f$
    /// and their derivatives \f$\phi_{,\xi}^k\f$, \f$\phi_{,\xi,\xi}^k\f$
    ///
    /// \param[in]  xi         Local coordinate at which evaluation is done
    ///                        \f$\xi^\alpha\f$
    /// \param[out] phi        Result of all shape functions
    ///                        \f$[\phi^K]\f$
    /// \param[out] dPhiDXi    First derivative of shape functions
    ///                        \f$[\partial\varphi^K / \partial\xi^\alpha]^T\f$
    /// \param[out] ddPhiDDXi  2nd derivative of shape functions
    ///                        \f$\varphi^K{}_{,\alpha\beta}\f$
    ///                        \f$=\partial^2\varphi^K / \partial\xi^\alpha\partial\xi^\beta\f$.
    ///                        The result is stored
    ///                        \f$[\varphi^K{}_{,11}, \varphi^K{}_{,12}, \varphi^K{}_{,22}]^T\f$.
    void evaluateGradHess( const VecLDim & xi,
                           Vec & phi,
                           Mat & dPhiDXi,
                           Mat & ddPhiDDXi );

    //@}


    /// Set positions
    ///
    /// \param[in]   coords   Positions to be assumed
    void setPositions( const Eigen::MatrixXd & coords );

    /// Tangents
    ///
    /// \param[in]     vIndex    Local element-wise index of node
    /// \param[out]    coeff0    Coefficients of first tangent w.r.t. vertices
    /// \param[out]    coeff1    Coefficients of second tangent w.r.t. vertices
    void computeTangents( const unsigned vertexId,
                          Vec & coeff0, Vec & coeff1 ){}

    //@name Debug methods
    //@{

    /// Print patch's topology
    void write( std::ostream & os );

    /// A temporary job for debugging purposes
    void tempJob( );

    //@}

private:
    /// Create SMF and TG file stream containing the patch.
    ///
    /// The SMF stream contains
    ///    1. #patchVertices_  : IDs = [ 0, patchVertices_.size() )
    ///    2. #extraVertices_  : IDs = [ patchVertices_.size(),
    ///                                  patchVertices_.size() + extraVertices_.size() )
    ///    3. #patchEdges_    : first is #edge_
    ///    4. #extraEdges_
    ///
    /// \param[out]  smf   SMF stream containing vertex co-ordinates
    ///                    and edge connectivity
    /// \param[out]  tg    TG stream containing vertex tags
    void writePatch_( std::ostream & smf,
                      std::ostream & tg );

    /// Compute number of needed subdivision to have parametric evaluation
    /// point inside of a regular patch
    ///
    /// \param[in]   xi   Local evaluation co-ordinate
    unsigned neededSubdivisions_( const VecLDim & xi );

    /// Create mesh for local patch
    ///
    /// \param[in,out] patch    The Patch mesh to be filled
    void formPatchMesh_( Mesh_ & patch );

    /// Compute subdivision matrix mapping vertices/DOFs in
    /// patch to regular sub-region
    ///
    /// \param[in]     xi       Evaluation point with respect to #edge_
    /// \param[out]    subMat   Subdivision matrix
    /// \param[out]    subXi    Transformed evaluation point
    ///                         with respect to regular sub-region   
    /// \param[out]    subJac   Transformed Jacobian at evaluation point
    ///                         with respect to regular sub-region   
    void computeSubdivisionMatrix_( const VecLDim & xi,
                                    Mat & subMat,
                                    VecLDim & subXi,
                                    MatLDimLDim & subJac );

private:
    /// 
    Edge *            edge_;
    /// Subdivision scheme
    Subdiv *          subdiv_;
    /// The limit surface spline
    Spline_           spline_;

    /// Patch of edges, i.e. one-neighbourhood of #edge_ (including #edge_ itself)
    VecFPtr           patchEdges_;
    /// Vertices in patch. These make up the degree-of-freedom.
    SetVPtr           patchVertices_;
    /// Additional edges, i.e. two-neighbourhood of #edge_, but without #patchEdges_
    VecFPtr           extraEdges_;
    /// Additional vertices attached only to #extraEdges_
    SetVPtr           extraVertices_;
};

//------------------------------------------------------------------------------
#include "ShapeFunSubdivision.ipp"

#endif
