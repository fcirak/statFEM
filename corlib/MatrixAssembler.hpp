// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//
//                 Computational Structural Mechanics Lab
//                         University of Cambridge
//
// This software is copyrighted and all rights are retained by the University of
// Cambridge. All use, disclosure, and/or reproduction of any part not expressly
// authorised by the University of Cambridge is prohibited. (C) 2010
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

//! @file MatrixAssembler.hpp

#ifndef corlib_matrixassembler_h
#define corlib_matrixassembler_h

#include <boost/function.hpp>
#include <boost/bind.hpp>
#include <boost/type_traits.hpp>

#include <Eigen/Core>

#include <corlib/ComputeElementMatrix.hpp>
#include <corlib/verify.hpp>
#include <corlib/eigenX.hpp>

namespace corlib{

    //--------------------------------------------------------------------------
    template<typename GIVMAT, typename GIVDOF, typename INSMAT, bool TRANS>
    class MatrixAssembler;


    //--------------------------------------------------------------------------
    /** \brief    Construct and store a MatrixAssembler for quicker coding
     *  \details  Creates all necessary object for the assembly and an 
     *  assembler object. Application via function call operator.
     *  \tparam ELEMENT  Type of element 
     *  \tparam QUAD     Type of quadrature
     *  \tparam SOLVER   Type of solver
     *  \tparam TRANS    True for transposition before assembly
     */
    template<typename ELEMENT, typename QUAD, typename SOLVER, bool TRANS=false>
    struct MatrixComputeAndAssembleFun
        : public boost::function< void( ELEMENT* ) >
    {
    public:
        //! @name Convenience typedefs
        //@{
        typedef boost::function<void(ELEMENT*, 
                                     std::vector<unsigned>&)>          GiveDof;
        typedef boost::function<void(ELEMENT*,
                                     const typename ELEMENT::VecLDim&,
                                     const double&,
                                     Eigen::MatrixXd&)>              Integrand;
        typedef corlib::ComputeElementMatrix<QUAD, Integrand>   MatrixComputer;
        typedef corlib::MatrixAssembler<MatrixComputer,GiveDof,
                                        SOLVER,TRANS>          MatrixAssembler;
        //@}

        //----------------------------------------------------------------------
        /*  \brief Cstor for the setup of a MatrixAssembler object
         *  \details  Since the setup of a MatrixAssembler requires quite many
         *  lines of code in the application, which repeat among most
         *  applications, this object attempts to shorten the procedure.
         *  In the constructor all objects, to which the MatrixAssembler holds
         *  references, are constructed and then the assembler object is 
         *  created. The function call operator applies the assembler to an
         *  element pointer.
         *  \param[in] intFun      Function pointer to the integrand
         *  \param[in] rowDofFun   Function pointer to the row dof numbers
         *  \param[in] colDofFun   Function pointer to the column dof numbers
         *  \param[in] quadrature  Reference to the quadrature created in main
         *  \param[in] solver      Pointer to the solver created in main
         */
        MatrixComputeAndAssembleFun( Integrand   intFun,
                                     GiveDof     rowDofFun,
                                     GiveDof     colDofFun,
                                     QUAD      & quadrature,
                                     SOLVER    * solver, 
                                     const double factor = 1. )
            : integrand_(       intFun ),
              matrixComputer_(  quadrature, integrand_ ),
              rowDofs_(         rowDofFun ), 
              colDofs_(         colDofFun ), 
              matrixAssembler_( matrixComputer_, rowDofs_, colDofs_, solver, factor )
        {}

        //! Overloaded function call operator applied to an element pointer
        void operator()( ELEMENT * ep )
        {
            matrixAssembler_( ep );
        }

    private:
        Integrand       integrand_;       //!< Functor of the element's integrand
        MatrixComputer  matrixComputer_;  //!< Functor for matrix computation
        GiveDof         rowDofs_;         //!< Functor for row dof numbers
        GiveDof         colDofs_;         //!< Functor for col dof numbers
        MatrixAssembler matrixAssembler_; //!< Functor to assemble the matrix
    };

    //--------------------------------------------------------------------------
    /** \brief    Construct and store a MatrixAssembler for quicker coding
     *  \details  Creates all necessary object for the assembly and an 
     *  assembler object. Application via function call operator.
     *  \tparam ELEMENT  Type of element 
     *  \tparam QUAD     Type of quadrature
     *  \tparam SOLVER   Type of solver
     */
    template<typename ELEMENT, typename SOLVER, bool TRANS=false>
    struct MatrixGiveAndAssembleFun
        : public boost::function< void( ELEMENT* ) >
    {
    public:
        //! @name Convenience typedefs
        //@{
        typedef boost::function<void(ELEMENT*, 
                                     std::vector<unsigned>&)>           GiveDof;
        typedef boost::function<void( const ELEMENT*, 
                                      Eigen::MatrixXd & )>          MatrixGiver;
        typedef corlib::MatrixAssembler<MatrixGiver,GiveDof,
                                        SOLVER,TRANS>           MatrixAssembler;

        //@}

        //----------------------------------------------------------------------
        /*  \brief Cstor for the setup of a MatrixAssembler object
         */
        MatrixGiveAndAssembleFun( MatrixGiver   giveFun,
                                  GiveDof       rowDofFun, 
                                  GiveDof       colDofFun,
                                  SOLVER      * solver,
                                  const double factor = 1. )
            : giver_(           giveFun ),
              rowDofs_(         rowDofFun ), 
              colDofs_(         colDofFun ), 
              matrixAssembler_( giver_, rowDofs_, colDofs_, solver, factor )
        {}

        //! Overloaded function call operator applied to an element pointer
        void operator()( ELEMENT * ep )
        {
            matrixAssembler_( ep );
        }

    private:
        MatrixGiver     giver_;           //!< Functor to give a matrix
        GiveDof         rowDofs_;         //!< Functor for row dof numbers
        GiveDof         colDofs_;         //!< Functor for col dof numbers
        MatrixAssembler matrixAssembler_; //!< Functor to assemble the matrix
    };

}

//------------------------------------------------------------------------------
/** \brief Matrix assembler based on functors
 *  \details Using the functors which provide the element's matrix and the 
 *  corresponding row and column degrees of freedom, this object collects all
 *  this data and uses the fourth functor, the matrix acceptor, for passing
 *  the data on. An additional factor is used for multiplying the matrix with
 *  a scalar.
 *  \tparam GIVMAT  Matrix donator type: functor providing the matrix
 *  \tparam GIVDOF  Degree of freedom donator type: gives the dof indices
 *  \tparam INSMAT  Matrix acceptor type: receives the matrices and assembles
 */
template<typename GIVMAT, typename GIVDOF, typename SOLVER, bool TRANS = false>
class corlib::MatrixAssembler
    : public boost::function< void( typename GIVMAT::first_argument_type ) >
{
public:
    //! @name Template parameters
    //@{
    typedef GIVMAT  GiveMatrix;
    typedef GIVDOF  GiveDof;
    typedef SOLVER  Solver;
    static const bool isTransposed = TRANS;
    //@}

    //! Constructor with function pointer to element matrix access
    MatrixAssembler( GiveMatrix   & giveMatrix, 
                     GiveDof      & rowDof,
                     GiveDof      & colDof,
                     SOLVER       * solver, 
                     const double factor = 1. )
        : giveMatrix_(   giveMatrix ),
          rowDof_(       rowDof ),
          colDof_(       colDof ),
          solver_(       solver ), 
          factor_(       factor )        
    {
        FTL_VERIFY( solver_ not_eq NULL );
    }

    //! Convenience typedef
    typedef typename GiveMatrix::first_argument_type ElementPtr;

    //! Function call operator for assembly
    void operator()( ElementPtr ep ) 
    {
        //! - get dof indices for rows and columns
        std::vector<unsigned> rowDofIndices, colDofIndices;
        rowDof_( ep, rowDofIndices );
        colDof_( ep, colDofIndices );

        //! - deduce matrix size from dofs: in the case of isTransposed==true
        //!   the element matrix is the transposed of the target matrix
        const unsigned numRows = 
            ( isTransposed ? colDofIndices.size() : rowDofIndices.size() );
        const unsigned numCols = 
            ( isTransposed ? rowDofIndices.size() : colDofIndices.size() );

        //! - storage for the element stiffness matrix
        Eigen::MatrixXd elementMatrix( numRows, numCols );
        elementMatrix.setZero();

        //! - obtain matrix from functor
        giveMatrix_( ep, elementMatrix );

        //! - multiply with factor
        elementMatrix *= factor_;

        //! - pass matrix to the solver:
        //!   inline transpose is used for isTransposed==true
        if ( TRANS ) 
            solver_ -> insertToMatrix( elementMatrix.transpose( ),
                                       rowDofIndices, colDofIndices );
        else
            solver_ -> insertToMatrix( elementMatrix, 
                                       rowDofIndices, colDofIndices );

        return;
    }

private:
    GiveMatrix    & giveMatrix_;   //!< Matrix giving functor
    GiveDof       & rowDof_;       //!< Row dof index giving functor
    GiveDof       & colDof_;       //!< Column dof index giving functor
    Solver        * solver_;       //!< Pointer to solver object

    //! Factor to multiply the matrix with
    const double factor_;
};
//------------------------------------------------------------------------------

#endif


