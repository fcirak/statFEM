// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//
//                 Computational Structural Mechanics Lab
//                         University of Cambridge
//
// This software is copyrighted and all rights are retained by the University of
// Cambridge. All use, disclosure, and/or reproduction of any part not expressly
// authorised by the University of Cambridge is prohibited. (C) 2010
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

//! @file MatrixLumper.hpp

#ifndef corlib_matrixlumper_h
#define corlib_matrixlumper_h

#include <boost/function.hpp>
#include <boost/bind.hpp>
#include <boost/type_traits.hpp>

#include <Eigen/Core>

#include <corlib/eigenX.hpp>
#include <corlib/ComputeElementMatrix.hpp>

//------------------------------------------------------------------------------
namespace corlib{

    template< typename RETMAT, typename RETNNOD, typename RETNOD,
              typename RETROWIND, typename STOVAL >
    class MatrixLumper;

    template< typename ELEMENT, typename QUAD >
    class MatrixComputeAndLumpFun;

    template< typename ELEMENT >
    class MatrixGiveAndLumpFun;

    namespace eigenX = corlib::eigenX;
    
}

//==============================================================================
/// Convenience function to make a MatrixLumper using element and node member function
///
/// This convenience functions provides an interface to construct a 
/// #MatrixLumper providing element and node member functions.
/// This is the typical interface used for lumping mass matrices, which are integrated
/// on-the-fly using quadrature QUAD, and "assembling" them on the nodes.
///
/// \tparam  ELEMENT   Element type
/// \tparam  QUAD      Quadrature type
template< typename ELEMENT, typename QUAD >
class corlib::MatrixComputeAndLumpFun
    : public boost::function< void( ELEMENT * ) >
{
public:
    // convenience type definitions
    typedef boost::function< void ( const ELEMENT *, 
                                    const typename ELEMENT::VecLDim &,
                                    const double &, 
                                    Eigen::MatrixXd & ) >                   MatrixIntegrand;
    typedef corlib::ComputeElementMatrix< QUAD, MatrixIntegrand >           MatrixComputer;
    typedef boost::function< unsigned ( const ELEMENT * ) >                 ReturnNumNodes;
    typedef boost::function< typename ELEMENT::Node * ( const ELEMENT *,
                                                        const unsigned ) >  ReturnNodePtr;
    typedef boost::function< unsigned ( const unsigned ) >                  ReturnRowIndex;
    typedef boost::function< void ( typename ELEMENT::Node *,
                                    const double & ) >                      StoreLumpedValue;
    typedef corlib::MatrixLumper< MatrixComputer, ReturnNumNodes,
                                  ReturnNodePtr, ReturnRowIndex,
                                  StoreLumpedValue >                        MatrixLumper;

    //! Constructor
    MatrixComputeAndLumpFun( MatrixIntegrand  matrixIntegrandFun,
                             ReturnNumNodes   returnNumNodesFun,
                             ReturnNodePtr    returnNodePtrFun,
                             StoreLumpedValue storeLumpedValueFun,
                             QUAD &           quadrature,
                             const unsigned   component = 0 )
        : matrixIntegrandFun_ ( matrixIntegrandFun ),
          matrixComputerFun_  ( quadrature, matrixIntegrandFun_ ),
          returnNumNodesFun_  ( returnNumNodesFun ),
          returnNodePtrFun_   ( returnNodePtrFun ),
          returnRowIndexFun_  ( boost::bind( std::plus< unsigned >(),
                                             boost::bind( std::multiplies< unsigned >(),
                                                          _1, ELEMENT::Node::dof ),
                                             component ) ),
          storeLumpedValueFun_( storeLumpedValueFun ),
          matrixLumper_       ( matrixComputerFun_, returnNumNodesFun_, returnNodePtrFun_,
                                returnRowIndexFun_, storeLumpedValueFun_ )
    { }

    //! Overloaded function call operator applied to an element pointer
    void operator()( ELEMENT * element )
    {
        matrixLumper_( element );
    }

private:
    MatrixIntegrand   matrixIntegrandFun_;     //!< Functor of the element's integrand
    MatrixComputer    matrixComputerFun_;      //!< Functor for matrix computation
    ReturnNumNodes    returnNumNodesFun_;      //!< Element method to get its number of (supporting) nodes / shape functions
    ReturnNodePtr     returnNodePtrFun_;       //!< Element method to get its node pointer by element-wise index
    ReturnRowIndex    returnRowIndexFun_;      //!< Function to get row index to lump
    StoreLumpedValue  storeLumpedValueFun_;    //!< Node method to assemble lumped value (e.g. point mass)
    MatrixLumper      matrixLumper_;           //!< The lumper
};

//==============================================================================
/// Convenience function to make a MatrixLumper using element and node member function
///
/// This convenience functions provides an interface to construct a 
/// #MatrixLumper providing element and node member functions.
/// This is the typical interface used for lumping mass matrices, which are held
/// readily on the element, and "assembling" them on the nodes.
///
/// \tparam  ELEMENT   Element type
template< typename ELEMENT >
class corlib::MatrixGiveAndLumpFun
    : public boost::function< void( ELEMENT* ) >
{
public:
    // convenience type definitions
    typedef boost::function< void ( const ELEMENT *,
                                    Eigen::MatrixXd & ) >                   ReturnMatrix;
    typedef boost::function< unsigned ( const ELEMENT * ) >                 ReturnNumNodes;
    typedef boost::function< typename ELEMENT::Node * ( const ELEMENT *,
                                                        const unsigned ) >  ReturnNodePtr;
    typedef boost::function< unsigned ( const unsigned ) >                  ReturnRowIndex;
    typedef boost::function< void ( typename ELEMENT::Node *,
                                    const double & ) >                      StoreLumpedValue;
    typedef corlib::MatrixLumper< ReturnMatrix, ReturnNumNodes,
                                  ReturnNodePtr, ReturnRowIndex,
                                  StoreLumpedValue >                        MatrixLumper;

    //! Constructor
    MatrixGiveAndLumpFun( ReturnMatrix     returnMatrixFun,
                          ReturnNumNodes   returnNumNodesFun,
                          ReturnNodePtr    returnNodePtrFun,
                          StoreLumpedValue storeLumpedValueFun,
                          const unsigned   component = 0 )
        : returnMatrixFun_    ( returnMatrixFun ),
          returnNumNodesFun_  ( returnNumNodesFun ),
          returnNodePtrFun_   ( returnNodePtrFun ),
          returnRowIndexFun_  ( boost::bind( std::plus< unsigned >(),
                                             boost::bind( std::multiplies< unsigned >(),
                                                          _1, ELEMENT::Node::dof ),
                                             component ) ),
          storeLumpedValueFun_( storeLumpedValueFun ),
          matrixLumper_       ( returnMatrixFun_, returnNumNodesFun_, returnNodePtrFun_,
                                returnRowIndexFun_, storeLumpedValueFun_ )
    { }

    //! Overloaded function call operator applied to an element pointer
    void operator()( ELEMENT * element )
    {
        matrixLumper_( element );
    }

private:
    ReturnMatrix      returnMatrixFun_;        //!< Element method to retrieve matrix (e.g. mass matrix)
    ReturnNumNodes    returnNumNodesFun_;      //!< Element method to get its number of (supporting) nodes / shape functions
    ReturnNodePtr     returnNodePtrFun_;       //!< Element method to get its node pointer by element-wise index
    ReturnRowIndex    returnRowIndexFun_;      //!< Function to get row index to lump
    StoreLumpedValue  storeLumpedValueFun_;    //!< Node method to assemble lumped value (e.g. point mass)
    MatrixLumper      matrixLumper_;           //!< The lumper
};

//==============================================================================
/// Lumps element matrices into scalars and stores it at nodes
///
/// The functor takes the first four function objects to accumulate
/// the row of an element matrix and the fifth function stores the
/// value at an element node.
///
/// Note, the element matrix must be a square matrix. Moreover, its dimensions
/// need to be #returnNumNodes_ * Element::Node::dof.
///
/// \tparam  RETMAT     Function type to return element's matrix
/// \tparam  RETNNOD    Function type to get element's number of element nodes
/// \tparam  RETNOD     Function type to get element's node pointer by element-wise index
/// \tparam  RETROWIND  Function type to get row index to lump
/// \tparam  STOVAL     Function type to store lumped value at node
template< typename RETMAT, typename RETNNOD, typename RETNOD, typename RETROWIND, typename STOVAL >
class corlib::MatrixLumper
    : public boost::function< void( typename RETMAT::first_argument_type ) >
{
public:
    typedef RETMAT                                    ReturnMatrix;
    typedef RETNNOD                                   ReturnNumNodes;
    typedef RETNOD                                    ReturnNodePtr;
    typedef RETROWIND                                 ReturnRowIndex;
    typedef STOVAL                                    StoreLumpedValue;

    typedef typename ReturnMatrix::first_argument_type           ElementPtr;
    typedef typename ReturnNodePtr::result_type                  NodePtr;
    typedef typename boost::remove_pointer< ElementPtr >::type   Element;

public:
    /// Constructor storing the different functions
    MatrixLumper( ReturnMatrix      returnMatrix,   
                  ReturnNumNodes    returnNumNodes, 
                  ReturnNodePtr     returnNodePtr,
                  ReturnRowIndex    returnRowIndex,
                  StoreLumpedValue  storeLumpedValue )
        : returnMatrix_    ( returnMatrix ),
          returnNumNodes_  ( returnNumNodes ), 
          returnNodePtr_   ( returnNodePtr ),
          returnRowIndex_  ( returnRowIndex ),
          storeLumpedValue_( storeLumpedValue )
    { }

    /// Method to lump values
    void operator()( ElementPtr elementPtr )
    {
        // get matrix
        const unsigned matSize = returnNumNodes_( elementPtr ) * Element::Node::dof;
        Eigen::MatrixXd elementMatrix( matSize, matSize, 0. );
        returnMatrix_( elementPtr, elementMatrix );

        // for each element node
        for ( unsigned n = 0; n < returnNumNodes_( elementPtr ); ++n ) {
            unsigned rowIndex = returnRowIndex_( n );
            const Eigen::VectorXd matrixRow = elementMatrix.row( rowIndex );
            const double lumpedRow =
                std::accumulate( matrixRow.begin( ), matrixRow.end( ), 0.0 );
            // assemble lumped mass on node if incremental function
            NodePtr nodePtr = returnNodePtr_( elementPtr, n );
            storeLumpedValue_( nodePtr, lumpedRow );
        }
    }

private:
    ReturnMatrix      returnMatrix_;
    ReturnNumNodes    returnNumNodes_;
    ReturnNodePtr     returnNodePtr_;
    ReturnRowIndex    returnRowIndex_;
    StoreLumpedValue  storeLumpedValue_;

};

#endif
