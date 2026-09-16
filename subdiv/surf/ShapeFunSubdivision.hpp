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

#ifndef subdiv_surf_shapefunsubdivision_h
#define subdiv_surf_shapefunsubdivision_h

#include <iostream>

#include <Eigen/Core>
#include <corlib/eigenX.hpp>

#include <corlib/misc.hpp>

#include <subdiv/surf/collect.hpp>
#include <subdiv/surf/Edge.hpp>
#include <subdiv/surf/FacetTree.hpp>
#include <subdiv/surf/Mesh.hpp>
#include <subdiv/surf/ShapeFunHelpers.hpp>

//------------------------------------------------------------------------------
namespace subdiv{
    namespace surf{

        template< typename FACET, typename SUBDIV >
        class ShapeFunSubdivision;

        namespace eigenX = corlib::eigenX;
    }
}

//==============================================================================
/// Shape functions using subdivision
///
/// \tparam    FACET    Facet type
/// \tparam    SUBDIV   Subdivision scheme type
template< typename FACET, typename SUBDIV >
class subdiv::surf::ShapeFunSubdivision
{
public:
    typedef FACET                                             Facet;
    typedef SUBDIV                                            Subdiv;

    static const corlib::shape myShape    = Facet::myShape;
    static const unsigned numVertices     = Facet::numVertices;
    static const unsigned numNeighbors    = Facet::numNeighbors;
    static const unsigned numEdges        = Facet::numEdges;
    /// number of shape functions over regular patch
    static const unsigned numFunctions    = Subdiv::numFunctions;
    /// manifold dimension
    static const unsigned localDim        = corlib::ShapeTraits< myShape >::dim;
    /// Number of effective second derivatives subtracting duplicities due to symmetry
    static const unsigned sDim            = localDim*(localDim+1)/2;
    //! Map to get Voigt indices from pair of local indices
    //! \f$(\alpha,\beta) \mapsto \mathcal{A}\f$
    //by \f$[localDim*\alpha+\beta]\f$ for any \f$\alpha,\beta=1,2\f$
    static const unsigned voigtForward[localDim][localDim];// = { { 0, 1 }, { 1, 2 } };

    static const unsigned numRings        = Subdiv::numRings;
    static const unsigned numExtraRings   = Subdiv::numExtraRings;

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
    typedef typename Facet::Vertex::VecDim                            VecDim;

    typedef std::vector< Facet * >                                    VecFPtr;
    typedef typename VecFPtr::iterator                                VecFPtrIter;

    typedef typename Facet::Vertex                                    Vertex;
    typedef std::set< Vertex * >                                      SetVPtr;
    typedef typename SetVPtr::iterator                                SetVPtrIter;
    typedef typename SetVPtr::const_iterator                          SetVPtrCIter;

private:
    /// Type of spline over regular patch
    typedef typename Subdiv::Spline                                   Spline_;

    /// Edge in patch mesh
    typedef subdiv::surf::Edge< Vertex, Facet >                       Edge_;
    /// Facet tree in patch mesh
    typedef subdiv::surf::FacetTree< myShape, Vertex >                FacetTree_;
    /// Mesh type of patch mesh
    typedef subdiv::surf::Mesh< Vertex, Edge_, Facet, FacetTree_ >    Mesh_;

    typedef typename Vertex::Point                                    Point_;
    typedef std::array< Vertex *, numFunctions >                      VecVPtrNF_;

public:
    /// Constructor
    ///
    /// \param[in]   f                  Facet on which subdivision shape functions 
    ///                                 are computed
    /// \param[in]   independentPatch   Make independent patch of surface
    ShapeFunSubdivision( Facet * f, const bool independentPatch = false );

    /// Destructor
    ~ShapeFunSubdivision( )
    {
        if ( indiePatch_ ) { delete indiePatch_; indiePatch_ = NULL; }
    }

    /// Return number of patch vertices
    unsigned numPatchVertices( ) const { return patchVertices_.size( ); }

    /// Return iterator to begin of patch vertices container
    SetVPtrCIter patchVerticesBegin( ) const { return patchVertices_.begin( ); }
    /// Return iterator to end of patch vertices container
    SetVPtrCIter patchVerticesEnd( ) const { return patchVertices_.end( ); }

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

    /// Limit position
    ///
    /// \param[in]     vIndex    Local facet-wise index of node
    /// \param[out]    coeff     Coefficients of limit position w.r.t. vertices
    void evaluate( const unsigned vertexId,
                   Vec & coeff );

    /// Tangents
    ///
    /// \param[in]     vIndex    Local facet-wise index of node
    /// \param[out]    coeff     Coefficients of tangents w.r.t. vertices
    void evaluateGradient( const unsigned vertexId,
                           Mat & coeff );

    //@}

    /// Set positions
    ///
    /// \param[in]   coords   Positions to be assumed
    void setPositions( const Eigen::MatrixXd & coords );

    //@name Debug methods
    //@{
    /// Print patch's topology
    void write( std::ostream & os );
    //@}

private:
    /// Make independent patch to surface mesh
    void liberatePatch_( );

    /// Create SMF and TG file stream containing the patch.
    ///
    /// The SMF stream contains
    ///    1. #patchVertices_  : IDs = [ 0, patchVertices_.size() )
    ///    2. #extraVertices_  : IDs = [ patchVertices_.size(),
    ///                                  patchVertices_.size() + extraVertices_.size() )
    ///    3. #patchFacets_    : first is #facet_
    ///    4. #extraFacets_
    ///
    /// \param[out]  smf   SMF stream containing vertex co-ordinates
    ///                    and facet connectivity
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
    /// \param[in]     xi       Evaluation point with respect to #facet_
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
    Facet *           facet_;
    /// Subdivision scheme
    Subdiv            subdiv_;
    /// The limit surface spline
    Spline_           spline_;

    /// Patch of facets, i.e. one-neighbourhood of #facet_ (including #facet_ itself)
    VecFPtr           patchFacets_;
    /// Vertices in patch. These make up the degree-of-freedom.
    SetVPtr           patchVertices_;
    /// Additional facets, i.e. two-neighbourhood of #facet_, but without #patchFacets_
    VecFPtr           extraFacets_;
    /// Additional vertices attached only to #extraFacets_
    SetVPtr           extraVertices_;
    
    /// Local patch to surface. Only exists if patch is copied at construction.
    Mesh_ *           indiePatch_;
};

//------------------------------------------------------------------------------
#include "ShapeFunSubdivision.ipp"

#endif
