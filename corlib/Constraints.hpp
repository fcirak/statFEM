// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//
//                 Computational Structural Mechanics Lab
//                         University of Cambridge
//
// This software is copyrighted and all rights are retained by the University of
// Cambridge. All use, disclosure, and/or reproduction of any part not expressly
// authorised by the University of Cambridge is prohibited. (C) 2010
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

//! @file Constraints.hpp
//! @todo   store contraints as vector<pair<Node*,Constraint> > globally

#ifndef corlib_constraints_h
#define corlib_constraints_h

//! system includes
#include <map>
#include <set>
#include <vector>
#include <utility>
#include <iterator>
#include <limits>
//! boost includes
#include <boost/function.hpp>
//! corlib includes
#include <corlib/verify.hpp>
#include <corlib/eigenX.hpp>
//------------------------------------------------------------------------------
namespace corlib{
    class NodalConstraint;

    //! @name Helper Functors to generate constraints
    //@{
    template<typename NODE> class ConstraintsFromFile;
    template<typename NODE, typename FUNC> class ConstraintsFromFunction;
    template<typename NODE, typename FUNC> class ConstraintsByID;
    template<typename NODE, typename SYSM> class ConstraintFunctor;
    //@}

}

//------------------------------------------------------------------------------
/** \brief Handles the constraints for one node
 *  \details A nodal constraint is represented by an activity bit and a value 
 *   for each degree of freedom. Each component can be set individually and
 *   the constraint can be applied to a given node.
 */
class corlib::NodalConstraint
{
public:
    //! Empty constructor 
    NodalConstraint() {}

    //! Constructor given the direction and value of constraint
    NodalConstraint( const unsigned dir, const double value )
    {
        constraints_.push_back( std::make_pair( dir, value ) );
    }

    //! Constructor given the whole set of values 
    template<unsigned DOF>
    NodalConstraint( const eigenX::VectorSd<DOF> & values )
    { 
        for ( unsigned d = 0; d < DOF; d ++ )
            constraints_.push_back( std::make_pair( d, values[d] ) );
    }
    
    //! Set a specific component of the constraint
    void setComponent( const unsigned dir, const double value )
    {
        // see if dir is already constrained
        unsigned index; // index of already constrained component
        const bool exists = this -> exists_( dir, index );

        if ( exists ) constraints_[ index ].second = value; // over-write
        else          constraints_.push_back( std::make_pair( dir, value ) );

        return;
    }

    //! Set all components 
    template<unsigned DOF>
    void setAll( const eigenX::VectorSd<DOF> & values )
    {
        constraints_.clear();
        for ( unsigned d = 0; d < DOF; d ++ )
            constraints_.push_back( std::make_pair( d, values[d] ) );
        return;
    }

    //! Fix completely
    void homogeneous( const unsigned dof )
    {
        constraints_.clear();
        for ( unsigned d = 0; d < dof; d ++ ) 
            constraints_.push_back( std::make_pair( d, 0. ) );
        return;
    }

    //! Apply the constraint to a node
    template<typename NODE>
    void applyToNode( NODE * np ) const
    {
        for ( unsigned s = 0; s < constraints_.size(); s ++ ) {
            np -> storeConstraint( constraints_[s].first, constraints_[s].second );
        }

        return;
    }

    //! Check total activity
    bool isActive() const { return constraints_.size(); } 

private:
    //! Check if direction has already been constrained
    bool exists_( const unsigned dir, unsigned & index)
    {
        bool result = false;
        for ( unsigned s = 0; s < constraints_.size(); s ++ ) {
            if ( dir == constraints_[s].first ) {
                index  = s;
                result = true;
            }
        }
        return result;
    }

private:
    //! Storage of pairs <component,value>
    std::vector< std::pair<unsigned,double> > constraints_;
};

//------------------------------------------------------------------------------
/** \brief Handle constraints that are given by an input file
 *  \details This constraint handler is based on an input stream with the
 *   <pre>
 *   format
 *   numConstraints
 *   N1  D1  V1
 *   .   .   .
 *   </pre>
 *   with the number of constraints given in the first line and each following 
 *   line contains the node number (N), the direction (D) and the
 *   value (V) to apply.
 *   \tparam NODE  The type of node, constraints are applied to
 */
template<typename NODE>
class corlib::ConstraintsFromFile
    : public boost::function<void( NODE * ) > 
{
private:
    static const unsigned dof_ = NODE::dof;
public:
    //! Constructor with given file stream inp
    ConstraintsFromFile( std::istream & inp ) 
    {
        //! get number of constraints
        unsigned numConstraints = 0;
        inp >> numConstraints;
        //! go through contraints
        for ( unsigned c = 0; c < numConstraints; c ++ ) {
            unsigned nodenum, component;
            double value;
            inp >> nodenum >> component >> value;
            if ( component < dof_ ) {
                typename MapNodeID2Constraint_::iterator pos = constraintContainer_.find( nodenum );
                if ( pos != constraintContainer_.end() ) {
                    (pos -> second).setComponent( component, value );
                }
                else {
                    constraintContainer_[ nodenum ] = corlib::NodalConstraint( component, value );
                }
            }
            inp.ignore( std::numeric_limits<std::streamsize>::max(), '\n' );
        }
    }

    //! Overloaded function for application of constraints to the nodes
    void operator()( NODE* np ) 
    {
        const unsigned id = np -> giveId();
        typename MapNodeID2Constraint_::iterator pos = constraintContainer_.find( id );
        if ( pos != constraintContainer_.end() ) {
            (pos -> second).applyToNode( np );
        }
        return;
    }

private:
    typedef std::map<unsigned, NodalConstraint >  MapNodeID2Constraint_;
    MapNodeID2Constraint_  constraintContainer_; //!< Container of constraints
};

//------------------------------------------------------------------------------
/** \brief Constraint handler using a function object
 *  \details Using a function object which returns a 'Constraint' object using
 *  a coordinate as argument.
 *  \tparam NODE  Type of node to apply the constraint to
 *  \tparam FUNC  Type of passed function object
 */
template<typename NODE, typename FUNC>
class corlib::ConstraintsFromFunction
    : public boost::function<void( NODE * ) >
{
public:
    //! Constructor given a suitable function object
    ConstraintsFromFunction( FUNC function ) : function_( function ) { }

    //! Overloaded function call in order to apply the constraint
    void operator()( NODE * np ) 
    {
        typename NODE::VecDim X = np -> giveCoordinates();
        corlib::NodalConstraint cstr = function_( X );
        if ( cstr.isActive() ) {
            cstr.applyToNode( np );
        }
        return;
    }

private:
    FUNC function_; //!< Function object of type 'Constraint = f(X)'
};

//------------------------------------------------------------------------------
/** \brief Select constraint nodes by ID-list and apply a function to them
 *  \details Passing a stream (containing a list of node numbers) to the
 *   constructor, this object memorizes the node IDs. Its function call operator
 *   applies to each of the nodes a function of type Constraint<dof>(VecDim X) 
 *   and applies the outcome as constraint.
 *  \tparam NODE Type of node to apply constraint to
 *  \tparam FUNC Type of function providing the constraint
 */
template<typename NODE, typename FUNC>
class corlib::ConstraintsByID 
    : public boost::function<void( NODE * ) >
{
public:
    //! Constructor given a suitable function object
    ConstraintsByID( std::istream & inp, 
                     FUNC function )
        : function_( function )
    {
        // fill set of constraint node ids with numbers from input stream
        std::istream_iterator<unsigned> inpIter( inp );
        std::istream_iterator<unsigned> eofIter;
 
        // commented out as it not passed for parallel applications when subdomain has no Dirichlet BC
        //FTL_VERIFY(inpIter != eofIter);

        constraintNodes_.insert(inpIter, eofIter);
    }

    //! Overloaded function call in order to apply the constraint
    void operator()( NODE * np ) 
    {
        const unsigned id = np -> giveId();
        if ( constraintNodes_.find( id ) != constraintNodes_.end() ) {
            typename NODE::VecDim X = np -> giveCoordinates();
            corlib::NodalConstraint constraint = function_( X );
            if ( constraint.isActive() ) {
                constraint.applyToNode( np );
            }
        }
        return;
    }

    unsigned numNodes() const { return constraintNodes_.size(); }

private:
    std::set<unsigned>            constraintNodes_;
    FUNC                          function_;
};

//------------------------------------------------------------------------------
//! Collects the constraints from the nodes and passes them to the system matrix
template<typename NODE, typename SYSM>
class corlib::ConstraintFunctor 
    : public boost::function<void( NODE * ) >
{
private:
    typedef SYSM *                                      SysMatPtr_;
    typedef std::vector< std::pair<unsigned, double> >  ConstraintVec_;

public:
    ConstraintFunctor( SysMatPtr_ sm, 
                       const double f = 1., 
                       const double s = 1. ) 
        : sysmat_( sm ), factor_( f ), scalar_( s )  { }

    ~ConstraintFunctor( ) 
    {
        sysmat_ = NULL;
        constraints_.clear( );
    }


    //! Collect the constraints in a temporary container 
    void operator( ) ( NODE * np ) 
    { 
        // get constraints from the node
        np -> giveConstraints( constraints_ ); 
        return;
    }

    //! Apply constraints to system matrix
    void applyConstraints( )
    {
        if ( not constraints_.empty( ) )
            sysmat_ -> applyConstraints( constraints_, factor_, scalar_ ); 
        return;
    }

private:
    SysMatPtr_      sysmat_;       //!< Access to the system matrix handler
    ConstraintVec_  constraints_;  //!< Temporary storage of constraints
    double          factor_;       //!< multiplier for scaling the constraints
    double          scalar_;       //!< weighting factor in the stiffness matrix 
};


#endif
