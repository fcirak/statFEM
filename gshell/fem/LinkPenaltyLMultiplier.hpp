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

#ifndef gshell_fem_linkpenaltylmultiplier_h
#define gshell_fem_linkpenaltylmultiplier_h

//------------------------------------------------------------------------------
#include <iostream>
#include <array>

#include <Eigen/Core>
#include <corlib/eigenX.hpp>

namespace gshell {
    namespace fem {

        template < typename NODE, unsigned DOF >
        class LinkPenaltyLMultiplier;

        namespace eigenX = corlib::eigenX;
    }
}

//==============================================================================
/// Basic link object which can be used to add a constraint of nodal DOFs
/// with:
///    - penalty method             (then #penalty_ > 0.)
///    - Lagrange multiplier        (then #dofMult_ must be set)
///    - augmented Lagrange method  (then both of above)
///
/// Residual and its gradient matrix:
///
/// \tparam  NODE   Node type whose DOFs are constraint
/// \tparam  DOF    Number of Lagrange multipliers
template< typename NODE, unsigned DOF >
class gshell::fem::LinkPenaltyLMultiplier
{
public:
    /// Node type
    typedef NODE                                              Node;
    /// Number of primary DOFs :  Lagrange multiplier(s)
    static const unsigned dof = DOF;

    /// Vector of primary DOFs: the one Lagrange multiplier
    typedef eigenX::VectorSd< dof >              VecDof;

protected:
    /// DOF index array (it's of length 1)
    typedef std::array< unsigned, dof >          IndexArray_;
    /// Vector of attaced nodes carrying DOFs
    typedef std::vector< Node* >                 VecNodePtr_;

    /// General vector
    typedef Eigen::VectorXd                      Vec_;
    /// General matrix
    typedef Eigen::MatrixXd                      Mat_;

private:
    /// Its own type
    typedef LinkPenaltyLMultiplier< NODE, DOF >               LinkPenaltyLMultiplier_;

public:
    /// Default constructor which does nothing
    LinkPenaltyLMultiplier( ) :
        nodes_( ),
        penalty_( 0. ),
        dofMult_( ),
        multiplier_( Eigen::VectorXd::Zero( dof ) ),
        increment_( Eigen::VectorXd::Zero( dof ) )
    {
        std::fill( dofMult_.begin( ), dofMult_.end( ),
                   std::numeric_limits< unsigned >::max( ) );
        return;
    }

    /// Virtual destructor to please the compiler
    virtual ~LinkPenaltyLMultiplier( )
    { 
        return;
    }

    /// Set penalty parameter
    void setPenalty( const double penalty )
    {
        penalty_ = penalty;
        return;
    }

    /// Assign numbers to the degree-of-freedom
    void numberMultDof( unsigned & counter )
    {
        for ( unsigned i = 0; i < dof; i++ ) dofMult_[ i ] = counter++;
        return;
    }

    /// Number of associated displacements DOFs
    unsigned numDofDisplacements( ) const
    {
        return Node::dof * nodes_.size( );
    }

    /// Give global DOF indices of linked DOF of stored nodes
    void getDofDisplacements( std::vector< unsigned > & dofIndices ) const
    {
        dofIndices.resize( Node::dof * nodes_.size() );
        std::vector< unsigned >::iterator iter = dofIndices.begin();
        FTL_VERIFY( this->numDofDisplacements( ) == resDisp_.size( ) );

        for ( unsigned n = 0; n < nodes_.size( ); ++n ) { 
            std::vector< unsigned > aux;
            nodes_[ n ]->copyDofArray( aux );
            // since it is a local copy
            iter = std::copy( aux.begin(), aux.end(), iter );
        }
        return;
    }

    /// Give global DOF index of Lagrange multiplier
    void getDofMultiplier( std::vector< unsigned > & dofIndices ) const
    {
        dofIndices.resize( dof );
        std::vector< unsigned >::iterator iter = dofIndices.begin();
        std::copy( dofMult_.begin( ), dofMult_.end( ), iter );
        return;
    }

    /// Return the associated Lagrange multiplier of link
    VecDof giveMultiplier( ) const
    {
        return ( multiplier_ + increment_ );
    }

    /// Set Lagrange multiplier
    void setMultiplier( const VecDof & mult )
    {
        multiplier_ = mult;
        return;
    }

    /// Update Lagrange multiplier and clear its increment
    void updateMultiplier( )
    {
        multiplier_ += increment_;
        increment_.setZero( );
        return;
    }

    /// Clear Lagrange multiplier increment
    void clearIncrement( )
    {
        increment_.clear( );
        return;
    }

    /// Set Lagrange multiplier increment
    void setIncrement( const VecDof & inc )
    {
        increment_ = inc;
        return;
    }

    /// Update Lagrange multiplier increment
    void addToIncrement( const VecDof & inc )
    {
        increment_ += inc;
        return;
    }

    /// Return negative residuum of link
    Eigen::VectorXd giveResidualDisp( ) const
    {
        return resDisp_;
    }

    /// Return negative RES contribution onto main field
    Eigen::VectorXd giveResidualMult( ) const
    {
        return resMult_;
    }

    /// Return negative residuum of link with scaled prescribed value
//    ublas::vector< double > giveResDispScaled( const double factor ) const;

    /// Return penalty matrix
    void giveMatrixDispDisp( Eigen::MatrixXd & matDispDisp ) const
    {
        matDispDisp = matDispDisp_;
    }

    /// Return link matrix containing DOF coupling
    void giveMatrixMultDisp( Eigen::MatrixXd & matMultDisp ) const
    {
        matMultDisp = matMultDisp_;
    }

    /// Return transposed link matrix containing DOF coupling
    void giveMatrixDispMult( Eigen::MatrixXd & matDispMult ) const
    { 
        matDispMult = matMultDisp_.transpose( );
    }

    /// Print contents
    std::ostream & write( std::ostream & os ) const
    {
        os << "LinkPenaltyLMultiplier" << std::endl;

        os << "    nodes=";
        std::transform( nodes_.begin( ), nodes_.end( ),
                        std::ostream_iterator< unsigned >( os, " " ),
                        std::bind( &Node::giveId, std::placeholders::_1 ) );
        os << std::endl;

        os << "    penalty=" << penalty_ << std::endl;
        
        os << "    dofMult=";
        std::copy( dofMult_.begin( ), dofMult_.end( ),
                   std::ostream_iterator< unsigned >( os, " " ) );
        os << std::endl;
        
        return os;
    }

protected:
    VecNodePtr_ nodes_;            ///< Attached nodes

    double      penalty_;          ///< penalty factor

    IndexArray_ dofMult_;          ///< index of Lagrange multiplier
    VecDof      multiplier_;       ///< Lagrange multiplier
    VecDof      increment_;        ///< Lagrange multiplier increment

    Vec_        resDisp_;          ///< Residual at nodal DOFs (displacements)
    VecDof      resMult_;          ///< Residual at Lagrange multiplier

    Mat_        matDispDisp_;      ///< derivative of #resDisp_ w.r.t. nodal DOFs (displacements)
    Mat_        matMultDisp_;      ///< derivative of #resMult_ w.r.t. nodal DOFs (displacements)
};

#endif
