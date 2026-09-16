// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//
//                 Computational Structural Mechanics Lab
//                         University of Cambridge
//
// This software is copyrighted and all rights are retained by the University of
// Cambridge. All use, disclosure, and/or reproduction of any part not expressly
// authorised by the University of Cambridge is prohibited. (C) 2010
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

//! @file ElementBasic.hpp

#ifndef corlib_elementbasic_h
#define corlib_elementbasic_h
//------------------------------------------------------------------------------

#include <array>
#include <boost/utility.hpp>

#include <Eigen/Core>
#include <corlib/Shape.hpp>
#include <corlib/eigenX.hpp>

namespace corlib{
    template< typename NODE, typename SFUN> class ElementBasic;

}

//------------------------------------------------------------------------------
/** \brief Basic element template
 *  \details  Has pointers to its node which can be set, returned 
 *  and accessed. It holds a shape function object which is used for the 
 *  geometry description
 *  \f[
 *       x(\xi) = \sum_k x^k \phi^k(\xi)
 *  \f]
 *  where the \f$ x^k \f$ are the coordinates of the stored nodes and
 *  \f$ \phi^k \f$ the shape functions.
 *  \tparam NODE   Type of node 
 *  \tparam SFUN   Type of shape function
 */
template<typename NODE, typename SFUN>
class corlib::ElementBasic
    : public boost::noncopyable
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
    static const unsigned      localDim = corlib::ShapeTraits<myShape>::dim;
    static const unsigned      numNodes = ShapeFun::numFunctions;
    //@}

    //--------------------------------------------------------------------------
    //! @name Convenience typedefs, publicly available
    //@{
    typedef eigenX::VectorSd<dim>                                VecDim;
    typedef eigenX::VectorSd<localDim>                           VecLDim;
    typedef eigenX::VectorSd<numNodes>                           VecNN;
    typedef eigenX::MatrixSd<dim, numNodes>                      MatDimNN;
    typedef eigenX::MatrixSd<localDim, numNodes>                 MatLDimNN;
    typedef Eigen::Matrix<VecNN, localDim, localDim>             MatVecNNLDimLDim;
    typedef Eigen::Matrix<VecNN, dim, dim>                       MatVecNNDimDim;
    //@}

    //--------------------------------------------------------------------------
    //! @name Node storage
    //@{
    typedef          std::array<Node*, numNodes>   NodeArray;
    typedef typename NodeArray::iterator           NodeIterator;
    typedef typename NodeArray::const_iterator     NodeConstIterator;
    //@}

public:
    //--------------------------------------------------------------------------
    //! @name Construction methods
    //@{

    //! Default constructor which sets the node pointers to NULL
    ElementBasic() { this -> clearNodes_(); }

    //! Read connectivity from stream inp, given access to the global node vector
    std::istream & readSelf( std::istream & inp, 
                             const std::vector<Node*> & nodes );
    //@}

    //! Virtual destructor to please the compiler
    virtual ~ElementBasic() { this -> clearNodes_(); }

    //--------------------------------------------------------------------------
    //! @name Accessors to Nodes, their indices and coordinates
    //@{
    
    //! Return iterator to begin of node array
    NodeIterator nodesBegin() { return nodes_.begin(); }
    //! Return iterator to end of node array
    NodeIterator nodesEnd(  ) { return nodes_.end(); }

    //! Return const iterator to begin of node array
    NodeConstIterator nodesBegin() const { return nodes_.begin(); }
    //! Return const iterator to end of node array
    NodeConstIterator nodesEnd(  ) const { return nodes_.end(); }

    //! Iterate over nodes
    template<typename OP> OP iterateOverNodes(OP op ) const;

    //! Iterate over nodes
    template<typename OP> OP iterateOverNodes(OP op );
    //@}

    //--------------------------------------------------------------------------
    //! @name Evaluation of shape function and gradient 
    //@{
    //! Return values of shape function at xi
    void sfun( const VecLDim & xi, VecNN & phi ) const;
    //! Return gradients of shape function at xi
    void sfunGrad( const VecLDim & xi, MatLDimNN & dphiDxi ) const;
    //! Return Hessians of shape functions at xi
    void sfunHess( const VecLDim & xi, MatVecNNLDimLDim & ddPhiDDXi ) const;
    //@}

    //--------------------------------------------------------------------------
    //! @name Kinematic quantities
    //@{
    //! Give Jacobian determinant
    double jacobian( const VecLDim & xi ) const;
    //! Compute global derivatives and return jacobian
    double globalDerivatives( const VecLDim & xi, MatDimNN & dphiDx ) const;
    //! Compute global derivatives of second order and return jacobian
    double globalSecondDerivatives( const VecLDim & xi, 
                                    MatVecNNDimDim & ddPhiDdx ) const;
    //! Geometry function 
    void geometry( const VecLDim & xi, VecDim & x ) const;
    //@}

private:
    //! Set all node pointers to NULL
    void clearNodes_( );
    //! Get nodal coordinates
    void nodalCoordinates_( MatDimNN & X ) const;

private:
    //--------------------------------------------------------------------------
    std::array<Node*, numNodes> nodes_;    //!< pointers to its nodes
};

#include "ElementBasic.ipp"

#endif
