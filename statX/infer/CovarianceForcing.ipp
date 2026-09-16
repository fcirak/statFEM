// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//
//                 Computational Structural Mechanics Lab
//                         University of Cambridge
//
// This software is copyrighted and all rights are retained by the University of
// Cambridge. All use, disclosure, and/or reproduction of any part not expressly
// authorised by the University of Cambridge is prohibited. (C) 2010
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

//! @file CovarianceForcing.ipp

//------------------------------------------------------------------------------
template <unsigned DOF, typename KERNEL>
void statX::infer::CovarianceForcing<DOF, KERNEL>::insertToRhs (
        const VectorLocal & eCovVec, const VectorUInt & dofIndices )
{
    FTL_VERIFY( dofIndices.size( ) == dof );
    for ( unsigned i = 0; i < dofIndices.size( ); i ++ ) {
        covVec_[ dofIndices[ i ] ] += eCovVec( i );
    }
    return;
}

//------------------------------------------------------------------------------
template <unsigned DOF, typename KERNEL>
void statX::infer::CovarianceForcing<DOF, KERNEL>::insertToCoordinates (
        const VecDim & nCoords )
{
    for ( unsigned i = 0; i < dof; i ++ ) {
        coords_.push_back( nCoords );
    }
    return;
}

//------------------------------------------------------------------------------
template <unsigned DOF, typename KERNEL>
void statX::infer::CovarianceForcing<DOF, KERNEL>::giveProjectedCovar (
        const HpsValueMap & hParams, const Eigen::MatrixXd & pMat,
        Eigen::MatrixXd & projCovmat ) const
{
    //! Firstly check that covVec_ and coords_ have the same size
    FTL_VERIFY( covVec_.size() == numNodes_ * dof );
    FTL_VERIFY( coords_.size() == numNodes_ * dof );
    
    //! generate triplet list of off-diagonal entries of C_f matrix
    std::vector < Eigen::Triplet<double> > tripletList;
    tripletList.reserve( ( unsigned ) 0.25 * dof * numNodes_ );

    for ( unsigned i = 0; i < numNodes_; ++ i ) {
        for ( unsigned j = i + 1; j < numNodes_; ++ j ) {
            
            //! assume independence between degrees of freedom at each node
            for ( unsigned k = 0; k < dof; ++ k ) {
                
                //! the row and column index
                const unsigned rindx = dof * i + k;
                const unsigned cindx = dof * j + k;
                
                //! take the respective coordinate and covVec
                const double phiI = covVec_[ rindx ];
                const VecDim x    = coords_[ rindx ];
                
                const double phiJ = covVec_[ cindx ];
                const VecDim y    = coords_[ cindx ];
                
                const double val = phiI * kernel_.giveValue( hParams, x, y ) * phiJ;
                
                //! cut off small values to obtain a sparse matrix
                if ( val < cutEps_ ) continue;

                tripletList.push_back(
                        Eigen::Triplet<double>( rindx, cindx, val ) );
                tripletList.push_back(
                        Eigen::Triplet<double>( cindx, rindx, val ) );
            }
        }
    }

    //! add diagonal entries to the triplet list
    for ( unsigned i = 0; i < numNodes_; ++ i ) {
        
        //! assume independence between degrees of freedom at each node
        for ( unsigned k = 0; k < dof; ++ k ) {
            
            //! the diagonal index
            const unsigned indx = dof * i + k;
            
            //! take the respective coordinate and covVec
            const double phiI = covVec_[ indx ];
            const VecDim x    = coords_[ indx ];
            
            const double val = phiI * kernel_.giveValue( hParams, x, x ) * phiI;
            
            //! cut off small values to obtain a sparse matrix
            if ( val < cutEps_ ) continue;
            
            tripletList.push_back( Eigen::Triplet<double>( indx, indx, val ) );
        }
    }

    const unsigned matSize = dof * numNodes_;
    Eigen::SparseMatrix<double> matCf( matSize, matSize );
    matCf.setFromTriplets( tripletList.begin( ), tripletList.end( ) );

    projCovmat = pMat * matCf * pMat.transpose( );

    return;
}
