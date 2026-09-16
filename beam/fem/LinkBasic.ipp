// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//
//                 Computational Structural Mechanics Lab
//                         University of Cambridge
//
// This software is copyrighted and all rights are retained by the University of
// Cambridge. All use, disclosure, and/or reproduction of any part not expressly
// authorised by the University of Cambridge is prohibited. (C) 2010
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

//! @file LinkBasic.ipp

//------------------------------------------------------------------------------
/** Accessor to the indices of the nodal degrees of freedom
 *  \param[in] dofIndices storage of dof indices
 */
template< typename NODE >
void beam::fem::LinkBasic<NODE>::
getDofIndices( std::vector<unsigned> & dofIndices ) const
{
    dofIndices.resize( nodeDofCoeffs_.size() );

    for ( unsigned linkNode = 0; linkNode < nodeDofCoeffs_.size(); ++linkNode ) {
        const Node * np = std::get< 0 >( nodeDofCoeffs_[ linkNode ] );
        const unsigned locDof = std::get< 1 >( nodeDofCoeffs_[ linkNode ] );
        std::vector<unsigned> globDofs( NODE::dof );
        np -> copyDofArray( globDofs );
        dofIndices[ linkNode ] = globDofs[ locDof ];
    }
     
    return;
}

//------------------------------------------------------------------------------
/** Accessor to the indices of the Lagrange multiplier degree-of-freedom
 *  \param[in] dofIndex  storage of dof index
 */
template< typename NODE >
void beam::fem::LinkBasic<NODE>::
getDofIndicesP( std::vector<unsigned> & dofIndex ) const
{
    dofIndex.resize( dofIndex_.size() );

    // copy to input iterator
    std::copy( dofIndex_.begin(), dofIndex_.end(), dofIndex.begin() );
    // done
    return;
}

//------------------------------------------------------------------------------
//! Return link matrix containing DOF coupling
template< typename NODE >
void beam::fem::LinkBasic<NODE>::giveCouplingMatrix( Eigen::MatrixXd & result ) const
{ 
    // number of (displacement) DOFs
    const unsigned numDof = nodeDofCoeffs_.size();

    // create coupling matrix
    for ( unsigned n = 0; n < numDof; ++n ) {
        result( n, 0 ) = std::get< 2 >( nodeDofCoeffs_[ n ] );
    }

}

//------------------------------------------------------------------------------
//! Return negative residuum of link
template< typename NODE >
Eigen::VectorXd beam::fem::LinkBasic<NODE>::giveLinkResiduum() const
{
    return this->giveLinkResiduumScaled( 1.0 );
}

//------------------------------------------------------------------------------
//! Return negative residuum of link with scaled prescribed value
template< typename NODE >
Eigen::VectorXd beam::fem::LinkBasic<NODE>::giveLinkResiduumScaled( const double factor ) const
{
    // number of (displacement) DOFs
    const unsigned numDof = nodeDofCoeffs_.size();

    // create negative residuum of link (ie constraint equation)
    VecDof res = factor * value_;
    for ( unsigned n=0; n<numDof; ++n ) {
        const Node * np = std::get< 0 >( nodeDofCoeffs_[ n ] );
        const unsigned locDof = std::get< 1 >( nodeDofCoeffs_[ n ] );
        const typename NODE::VecDof disp = np->giveDisplacements();
        res[ 0 ] -= std::get< 2 >( nodeDofCoeffs_[ n ] ) * disp[ locDof ];
    }
    return res;
}

//------------------------------------------------------------------------------
//! Return negative RHS contribution onto main field
template< typename NODE >
Eigen::VectorXd beam::fem::LinkBasic<NODE>::giveMultiplierRhs() const
{
    // number of (displacement) DOFs
    const unsigned numDof = nodeDofCoeffs_.size();

    // create coupling matrix
    Eigen::MatrixXd coupling( nodeDofCoeffs_.size(), 1 );
    this->giveCouplingMatrix( coupling );

    // create RHS contribution of Lagrange multiplier
    Eigen::VectorXd rhs( numDof );
    for ( unsigned d=0; d<numDof; ++d ) {
        const VecDof multiplier = this->giveMultiplier( );
        rhs[ d ] = - multiplier( 0 ) * coupling( d, 0 );
    }
    return rhs;
}

//------------------------------------------------------------------------------
//! Print contents
template< typename NODE >
std::ostream& beam::fem::LinkBasic<NODE>::print( std::ostream& os ) const
{
    os << "Link=" << dofIndex_[0] << std::endl;
    for ( unsigned el=0; el < nodeDofCoeffs_.size(); ++el ) {
        os << "    nodeID="  << (std::get< 0 >( nodeDofCoeffs_[ el ] ))->giveId()
           << ", DOF="   << std::get< 1 >( nodeDofCoeffs_[ el ] )
           << ", value=" << std::get< 2 >( nodeDofCoeffs_[ el ] )
           << std::endl;
    }
    os << "    value_=" << value_ << std::endl;
    os << "    LinkResiduum=" << this->giveLinkResiduum() << std::endl;
    os << std::endl;
    return os;
}

//------------------------------------------------------------------------------
//! Print connectivity 
template< typename NODE >
std::ostream& beam::fem::LinkBasic<NODE>::printDofs( std::ostream& os ) const
{
    std::vector< double > globalDofs( nodeDofCoeffs_.size() + 1 );
    std::vector< double >::iterator iter = globalDofs.begin();
    this->getDofIndices( iter );
    std::advance( iter, nodeDofCoeffs_.size() );
    this->getDofIndicesP( iter );
    for ( unsigned d=0; d<globalDofs.size(); ++d )
        os << globalDofs[ d ] << " ";
    os << std::endl;
    return os;
}
