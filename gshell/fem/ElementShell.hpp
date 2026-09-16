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
//! @date   2012

#ifndef gshell_fem_elementshell_h
#define gshell_fem_elementshell_h

#include <iostream>
#include <cassert>
#include <array>

#include <Eigen/Core>
#include <corlib/eigenX.hpp>

#include <corlib/Shape.hpp>

//------------------------------------------------------------------------------
namespace gshell{
    namespace fem{

        template< typename NODE, typename SFUN > class ElementShell;

        namespace eigenX = corlib::eigenX;
    }
}

//==============================================================================
/** \brief Basic element template
 *
 *  \details  Has pointers to #nodes_ and #supportNodes_ which can be set, returned 
 *  and accessed. The #nodes_ are the vertices (or corners) of the element's shape;
 *  The #supportNodes_ are the nodes used for spanning the (polynomial) space over the 
 *  element's shape together with the assocatied shape functions #shapeFun_;
 *
 *  The shape function object is used for the geometry description
 *  \f[
 *       x(\xi) = \sum_k x^k \phi^k(\xi)
 *  \f]
 *  where the \f$ x^k \f$ are the coordinates of the stored #supportNodes_ and
 *  \f$ \phi^k \f$ the shape functions.
 *
 *  Note, the number of #nodes_ is fixed and shaped by #myShape, but the number of
 *  #supportNodes_ is variable depending on the connectivity of the mesh.
 *
 *  \tparam NODE   Type of node 
 *  \tparam SFUN   Type of shape function
 */
template< typename NODE, typename SFUN >
class gshell::fem::ElementShell
{
public:
    //--------------------------------------------------------------------------
    //! @name Basic typedefs
    //@{
    typedef NODE Node;
    typedef SFUN ShapeFun;
    //@}

    //--------------------------------------------------------------------------
    //! @name Basic attributes
    //@{
    static const unsigned      dim      = Node::dim;
    static const corlib::shape myShape  = ShapeFun::myShape;
    static const unsigned      localDim = ShapeFun::localDim;
    //! Number of effective second derivatives subtracting duplicities due to symmetry -- always 3
    static const unsigned      sDim     = ShapeFun::sDim;
    //! Map to get Voigt indices from pair of local indices
    //! \f$(\alpha,\beta) \mapsto \mathcal{A}\f$
    static const unsigned      voigtForward[localDim][localDim];
    //! Number of vertices of simplex shape -- can be 3 or 4
    static const unsigned      numVertices = ShapeFun::numVertices;   // this clashes
    static const unsigned      numNodes = numVertices;                // this clashes
    //@}

    //--------------------------------------------------------------------------
    //! @name Convenience typedefs, publicly available
    //@{
    typedef eigenX::VectorSd<dim>                   VecDim;
    typedef eigenX::VectorSd<localDim>              VecLDim;
    typedef eigenX::MatrixSd<dim, numNodes>         MatDimNN;
    typedef Eigen::VectorXd                         VecNF;
    typedef Eigen::MatrixXd                         MatDimNF;
    typedef Eigen::MatrixXd                         MatLDimNF;
    typedef Eigen::MatrixXd                         MatSDimNF;
    //@}

    //--------------------------------------------------------------------------
    //! @name Node storage
    //@{
    typedef std::array< Node *, numNodes >                       NodeArray;
    typedef typename NodeArray::iterator                         NodeIterator;
    typedef typename NodeArray::const_iterator                   NodeConstIterator;
    typedef std::vector< Node * >                                NodeVector;
    //@}

public:
    //--------------------------------------------------------------------------
    //! @name Construction methods
    //@{

    //! Default constructor
    ElementShell() : shapeFun_(),
                     surfaceArea_( 0. )
    { 
        this -> clearNodes_();
    }

    //! Set shape function cache
    //!
    //! Note: This may be called several times.
    //!
    //! \tparam     SFUNPARAMETRIC  A parametric shape function
    //! \tparam     QUADRATURE      Quadrature type
    //! \param[in]  shapeFun        Pointer to shape function
    //! \param[in]  quadrature      Quadrature object
    template< typename SFUNPARAMETRIC, typename QUADRATURE >
    void populateShapeFunctionCache( SFUNPARAMETRIC * shapeFun,
                                     const QUADRATURE & quadrature )
    {
        shapeFun_.populateEvaluatedShapeFunctions( shapeFun, quadrature );
        shapeFun_.populateLimitCoefficientsAtVertices( shapeFun );
        shapeFun_.populateTangentCoefficientsAtVertices( shapeFun );
    }

    //! Store limit subdivision weights at node
    //!
    //! Note: This is only carried out for each node once.
    void setLimitCoefficientsAtNodes();

    //! Store unique subdivision weights for tangents at node
    //!
    //! Note: This is only carried out for each node once.
    void setUniqueTangentCoefficientsAtNodes();

    //! Set the node pointer ind to node
    //!
    //! \param[in] index Number of the nodes position
    //! \param[in] node  Pointer to node which will be the new node
    void setNodePtr( const unsigned ind, Node * node );

    //! Read connectivity from stream inp, given access to the global node vector
    //!
    //! \param  inp     Input stream
    //! \param  nodes   Container with node pointers
    std::istream & readSelf( std::istream & inp, 
                             const std::vector< Node * > & nodes );

    //! Set support node pointers
    //!
    //! \param[in]  supportNodes  Container with support node pointers
    void setSupportNodes( const NodeVector & supportNodes );

    //@}

    //! Virtual destructor to please the compiler
    virtual ~ElementShell()
    {
        this -> clearNodes_();
    }

    //--------------------------------------------------------------------------
    //! @name Accessors to vertex nodes, their indices and coordinates
    //@{
    
    //! Give node pointer for given index \details
    //! \param[in] ind     Index of the desired node 
    //! \return            Pointer to the node with number #ind
    Node * giveNodePtr( const unsigned ind ) const { return nodes_[ ind ]; }

    //! Return iterator to begin of node array
    NodeIterator nodesBegin() { return nodes_.begin(); }
    //! Return iterator to end of node array
    NodeIterator nodesEnd() { return nodes_.end(); }

    //! Return const iterator to begin of node array
    NodeConstIterator nodesBegin() const { return nodes_.begin(); }
    //! Return const iterator to end of node array
    NodeConstIterator nodesEnd() const { return nodes_.end(); }

    //! Write the IDs of the nodes to a stream
    //!
    //! \param[in,out] out Output stream
    //! \return            Output stream
    std::ostream & writeNodeIndices( std::ostream & out ) const;

    //! Write the IDs of the support nodes to a stream
    //!
    //! \param[in,out]  os   Output stream
    //! \return              Output stream
    std::ostream & writeSupportNodeIndices( std::ostream & out ) const;

    //! Assign the node indices to a given iterator.
    //!
    //! \tparam OUT     Type of the output iterator
    //! \param[in,out]  oIter Iterator to a storage of the node indices
    //! \return               The manipulated iterator
    template< typename OUT >
    OUT giveNodeIndices( OUT iter ) const;

    //! Apply a given functor to all nodes by means of the for_each algorithm
    //! (const version)
    //!
    //! \tparam OP       Type of node functor
    //! \param[in] op    Functor applicable to NODE*
    //! \return          Given functor op
    template< typename OP >
    OP iterateOverNodes( OP op ) const;
    
    //! Apply a given functor to all nodes by means of the for_each algorithm
    //! (non-const version)
    //!
    //! \tparam OP       Type of node functor
    //! \param[in] op    Functor applicable to NODE*
    //! \return          Given functor op
    template< typename OP > 
    OP iterateOverNodes( OP op );
    
    //! Collect matrix containing the co-ordinates of the vertices
    //!
    //! \param[out]   X    Reference co-ordinates of supporting nodes
    //!                    \f$[X^A{}_K]\f$
    void nodalCoordinates( MatDimNN & X ) const;

    //! Collect matrix containing the co-ordinates of the support nodes
    //!
    //! \param[out]   X    Reference co-ordinates of supporting nodes
    //!                    \f$[X^A{}_K]\f$
    void supportNodeCoordinates( MatDimNF & X ) const;
    //@}

    //--------------------------------------------------------------------------
    //! @name Accessors to support nodes 
    //@{

    //! Return number of supporting nodes \f$NF\f$
    unsigned numFunctions() const { return supportNodes_.size(); }

    //! Give node pointer for given index \details
    //!
    //! \param[in] ind       Local index of the desired node 
    //! \retval    NodePtr   Pointer to the node with number ind */
    Node * giveSupportNodePtr( const unsigned ind ) const
    {
        return supportNodes_[ ind ];
    }

    //@}

    //--------------------------------------------------------------------------
    //! @name Evaluation of shape function and gradient 
    //@{

    //! Evaluates the shape function at the given coordinate xi
    //!
    //! \param[in]  xi   Local cooordinate at which the functions is evaluated
    //! \param[out] phi  Result of all shape functions at xi 
    void sfun( const VecLDim & xi, VecNF & phi ) const;

    //! Compute the gradient w.r.t xi of the shape functions at the given xi:
    //!
    //! \param[in]  xi       Local coordinate at which the gradient is evaluated
    //!                      \f$\xi^\alpha\f$
    //! \param[out] dPhiDXi  First derivative of shape functions
    //!                      \f$ [\varphi^K{}_{,\alpha}]^T = [\partial\varphi^K / \partial\xi^\alpha]^T \f$
    void sfunGrad( const VecLDim & xi, MatLDimNF & dPhiDXi ) const;

    //! Compute the 2nd derivative w.r.t xi of the shape functions at the given xi:
    //!
    //! \param[in]  xi         Local coordinate at which the gradient is evaluated
    //!                        \f$\xi^\alpha\f$
    //! \param[out] ddPhiDDXi  2nd derivative of shape functions
    //!                        \f$\varphi^K{}_{,\alpha\beta}\f$
    //!                        \f$=\partial^2\varphi^K / \partial\xi^\alpha\partial\xi^\beta\f$.
    //!                        The result is stored
    //!                        \f$[\varphi^K{}_{,11}, \varphi^K{}_{,12}, \varphi^K{}_{,22}]^T\f$.
    void sfunHess( const VecLDim & xi, MatSDimNF & ddPhiDDXi ) const;

    //! Return value, 1st and 2nd derivatives of shape functions at xi
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
    void sfunGradHess( const VecLDim & xi, VecNF & phi,
                       MatLDimNF & dPhiDXi, MatSDimNF & ddPhiDDXi ) const;

    /// Return limit position
    ///
    /// \param[in]     nIndex    Local element-wise index of node
    /// \param[out]    coeff     Coefficients of limit position w.r.t. vertices
    void sfun( const unsigned nIndex, VecNF & coeff ) const;

    /// Return tangents
    ///
    /// \param[in]     nIndex    Local element-wise index of node
    /// \param[out]    coeff     Coefficients of tangents w.r.t. vertices
    void sfunGrad( const unsigned nIndex, MatLDimNF & coeff ) const;

    //@}

    //--------------------------------------------------------------------------
    //! @name Kinematic quantities
    //@{

    //! Interpolate reference surface (co-ordinate) at local point
    //!
    //! \param[in]      xi     Local point
    //! \param[out]     x      Interpolated co-ordinate
    void geometry( const VecLDim & xi, VecDim & x ) const;

    //! Interpolate reference surface (co-ordinate) at local point
    //!
    //! \param[in]      xi     Local point
    //! \param[in]      conf   Configuration
    //! \return                Interpolated co-ordinate
    VecDim giveCoordinate( const VecLDim & xi )
    {
        VecDim x;
        this -> geometry( xi, x );
        return x;
    }

    //! Interpolate reference surface (co-ordinate) at local corner
    //! 
    //! \param[in]     nIndex    Local element-wise index of node    
    //! \param[out]    x         Interpolated co-ordinate
    void geometry( const unsigned nIndex, VecDim & x ) const;

    //! Give surface metric w.r.t. reference configuration
    double jacobian( const VecLDim & xi ) const;

//    //! Compute global derivatives and return jacobian
//    double globalDerivatives( const VecLDim & xi, MatDimNF & dPhiDX ) const;

    //@}
    
    //--------------------------------------------------------------------------
    //! @name Element surface area
    //@{

    //! Clear area
    void clearArea( ) { surfaceArea_ = 0.; return; }

    //! Return the element's internal strain energy
    double giveArea( ) const { return surfaceArea_; }

    //! In-mid-surface area integrand
    //!
    //! \param[in] xi      Local coordinate at which this function is evaluated
    //! \param[in] weight  Corresponding quadrature weight
    void areaIntegrand( const VecLDim & xi, const double & weight );

    //@}

private:
    //! Set all node pointers to NULL
    void clearNodes_( );
    
private:
    //--------------------------------------------------------------------------
    NodeArray     nodes_;        //!< pointers to its vertex nodes
    NodeVector    supportNodes_; //!< pointers to the nodes 
    ShapeFun      shapeFun_;     //!< shape function object

    double        surfaceArea_;         
};

#include "ElementShell.ipp"

#endif
