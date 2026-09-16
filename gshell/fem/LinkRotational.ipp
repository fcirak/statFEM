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
//! @date   2011

//==============================================================================

//------------------------------------------------------------------------------
template< typename ELEMENT >
gshell::fem::LinkRotationalExt< ELEMENT >::LinkRotationalExt(
    Node * node,
    Element * element
    ) :
    node_( node ),
    element_( element ),
    normal_( Eigen::VectorXd::Zero( Node::dim ) )
{
    const unsigned numNodes = element_->numFunctions( );
    for ( unsigned i = 0; i< numNodes; ++i )
        nodes_.push_back( element_->giveSupportNodePtr( i ) );
    
    const unsigned i = LinkRotationalExt_::nodeLocalIndex_( element, node );
    element_->computeNormalAtNode( gshell::fem::REFERENCE, i, normal_, NULL );
    return;
}

//------------------------------------------------------------------------------
template< typename ELEMENT >
void gshell::fem::LinkRotationalExt< ELEMENT >::computeResidual( )
{
    const unsigned vId = LinkRotationalExt_::nodeLocalIndex_( element_, node_ );

    VecDim_ curNormal;
    Mat_ curNormalGradDof;
    element_->computeNormalAtNode( gshell::fem::CURRENT, vId,
                                   curNormal, &curNormalGradDof );
    const VecDim_ diffNormal = curNormal - normal_;

//    VecDim_ curTang0, curTang1;
//    Mat_ curTang0GradDof, curTang1GradDof;
//    element_->computeTangentsAtNode( gshell::fem::REFERENCE, vId, 
//                                     curTang0, curTang1,
//                                     &curTang0GradDof, &curTang1GradDof );
//    BOOST_STATIC_ASSERT( dof == 2 );
//    const VecDim_ curLag = ( multiplier_[ 0 ] * curTang0 +
//                             multiplier_[ 1 ] * curTang1 );

    const unsigned numDofDisp = this->numDofDisplacements( );
    FTL_VERIFY( curNormalGradDof.cols( ) == numDofDisp );
//    FTL_VERIFY( curTang0GradDof.cols( ) == numDofDisp );
//    FTL_VERIFY( curTang1GradDof.cols( ) == numDofDisp );

    // residuals
    resDisp_.resize( numDofDisp );
    resDisp_ = penalty_ * ( diffNormal * curNormalGradDof );
    //+
//        ublas::prod( curLag, curNormalGradDof ) +
//        multiplier_[ 0 ] * ublas::prod( diffNormal, curTang0GradDof ) +
//        multiplier_[ 1 ] * ublas::prod( diffNormal, curTang1GradDof );

//    resMult_[ 0 ] = ublas::inner_prod( diffNormal, curTang0 );
//    resMult_[ 1 ] = ublas::inner_prod( diffNormal, curTang1 );

    return;
}

//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
template< typename ELEMENT >
void gshell::fem::LinkRotationalExt< ELEMENT >::computeMatrix( )
{
    const unsigned vId = LinkRotationalExt_::nodeLocalIndex_( element_, node_ );

    VecDim_ curNormal;
    Mat_ curNormalGradDof;
    element_->computeNormalAtNode( gshell::fem::CURRENT, vId,
                                   curNormal, &curNormalGradDof );
//    const VecDim_ diffNormal = curNormal - normal_;

//    VecDim_ curTang0, curTang1;
//    Mat_ curTang0GradDof, curTang1GradDof;
//    element_->computeTangentsAtNode( gshell::fem::REFERENCE, vId, 
//                                     curTang0, curTang1,
//                                     &curTang0GradDof, &curTang1GradDof );
//    BOOST_STATIC_ASSERT( dof == 2 );
//    const VecDim_ curLag = ( multiplier_[ 0 ] * curTang0 +
//                             multiplier_[ 1 ] * curTang1 );

    const unsigned numDofDisp = this->numDofDisplacements( );
    FTL_VERIFY( curNormalGradDof.cols( ) == numDofDisp );
//    FTL_VERIFY( curTang0GradDof.cols( ) == numDofDisp );
//    FTL_VERIFY( curTang1GradDof.cols( ) == numDofDisp );

    // gradients
    matDispDisp_ = penalty_ * ( curNormalGradDof.transpose( ) * curNormalGradDof );
    //+
//        multiplier_[ 0 ] * ublas::prod( ublas::trans( curNormalGradDof ),
//                                        curTang0GradDof ) +
//        multiplier_[ 1 ] * ublas::prod( ublas::trans( curNormalGradDof ),
//                                        curTang1GradDof ) +
//        multiplier_[ 0 ] * ublas::prod( ublas::trans( curTang0GradDof ),
//                                        curNormalGradDof ) +
//        multiplier_[ 1 ] * ublas::prod( ublas::trans( curTang1GradDof ),
//                                        curNormalGradDof );

//    matMultDisp_.resize( dof, numDofDisp );
//    ublas::row( matMultDisp_, 0 ) = ublas::prod( diffNormal, curTang0GradDof ) +
//        ublas::prod( curTang0, curNormalGradDof );
//    ublas::row( matMultDisp_, 1 ) = ublas::prod( diffNormal, curTang1GradDof ) +
//        ublas::prod( curTang1, curNormalGradDof );

    return;
}

//------------------------------------------------------------------------------
template< typename ELEMENT >
void gshell::fem::LinkRotationalExt< ELEMENT >::computeResidualAndMatrix( )
{
    this->computeResidual( );
    this->computeMatrix( );
    return;
}

//------------------------------------------------------------------------------
template< typename ELEMENT >
std::ostream & gshell::fem::LinkRotationalExt< ELEMENT >::write(
    std::ostream & os
    ) const
{
    os << "LinkRotationalExt" << std::endl;
    os << "    nodeId=" << node_->giveId() << std::endl;
    os << "    normal=" << normal_ << std::endl;
    this->LinkPenaltyLMultiplier_::write( os );
    return os;
}




//==============================================================================

//------------------------------------------------------------------------------
template< typename ELEMENT >
gshell::fem::LinkRotationalInt< ELEMENT >::LinkRotationalInt(
    Node * node,
    Element * element0,
    Element * element1
    ) :
    node_( node ),
    element0_( element0 ),
    element1_( element1 ),
    normal0_( Eigen::VectorXd::Zero( Node::dim ) ),
    normal1_( Eigen::VectorXd::Zero( Node::dim ) )
{
    // add nodes of 1st element
    if ( element0_ )  this->computeNormal_( 0 );
    // add nodes of 2nd element
    if ( element1_ )  this->computeNormal_( 1 );

    return;
}

//------------------------------------------------------------------------------
template< typename ELEMENT >
void gshell::fem::LinkRotationalInt< ELEMENT >::computeResidual( )
{
    // current normal 0
    VecDim_ curNormal0;
    Mat_ curNormal0GradDof;
    const unsigned vId0 = LinkRotationalInt_::nodeLocalIndex_( element0_, node_ );
    element0_->computeNormalAtNode( gshell::fem::CURRENT, vId0,
                                    curNormal0, &curNormal0GradDof );

    // current normal 1
    VecDim_ curNormal1;
    Mat_ curNormal1GradDof;
    const unsigned vId1 = LinkRotationalInt_::nodeLocalIndex_( element1_, node_ );
    element1_->computeNormalAtNode( gshell::fem::CURRENT, vId1,
                                    curNormal1, &curNormal1GradDof );

    // difference of direction cosine
    const double diffAngle =
        curNormal1.dot( curNormal0 ) - normal1_.dot( normal0_ );

    // DOFs
    const unsigned numDofDisp = this->numDofDisplacements( );
    const unsigned numDofDisp0 = curNormal0GradDof.cols( );
    const unsigned numDofDisp1 = curNormal1GradDof.cols( );
    FTL_VERIFY( numDofDisp0 + numDofDisp1 == numDofDisp );

    // residuals
    resDisp_.resize( numDofDisp );
    resDisp_.head( numDofDisp0 ) = penalty_ * diffAngle * ( curNormal1 * curNormal0GradDof );
    resDisp_.segment( numDofDisp0, numDofDisp-numDofDisp0 ) = 
            penalty_ * diffAngle * ( curNormal0 * curNormal1GradDof );

    // resMult_ =

    return;
}

//------------------------------------------------------------------------------
template< typename ELEMENT >
void gshell::fem::LinkRotationalInt< ELEMENT >::computeMatrix( )
{
    // current normal 0
    VecDim_ curNormal0;
    Mat_ curNormal0GradDof;
    const unsigned vId0 = LinkRotationalInt_::nodeLocalIndex_( element0_, node_ );
    element0_->computeNormalAtNode( gshell::fem::CURRENT, vId0,
                                    curNormal0, &curNormal0GradDof );

    // current normal 1
    VecDim_ curNormal1;
    Mat_ curNormal1GradDof;
    const unsigned vId1 = LinkRotationalInt_::nodeLocalIndex_( element1_, node_ );
    element1_->computeNormalAtNode( gshell::fem::CURRENT, vId1,
                                    curNormal1, &curNormal1GradDof );

    // difference of direction cosine
    const double diffAngle =
        curNormal1.dot( curNormal0 ) - normal1_.dot( normal0_ );

    // DOFs
    const unsigned numDofDisp = this->numDofDisplacements( );
    const unsigned numDofDisp0 = curNormal0GradDof.cols( );
    const unsigned numDofDisp1 = curNormal1GradDof.cols( );
    FTL_VERIFY( numDofDisp0 + numDofDisp1 == numDofDisp );

    // gradients
    const Eigen::VectorXd n1ng0 = curNormal0GradDof.transpose() * curNormal1;
    const Eigen::VectorXd n0ng1 = curNormal1GradDof.transpose() * curNormal0;
    
    matDispDisp_.resize( numDofDisp, numDofDisp );
    matDispDisp_.block( 0, 0, numDofDisp0, numDofDisp0 ) = 
            penalty_ * ( n1ng0 * n1ng0.transpose( ) );
    matDispDisp_.block( 0, numDofDisp0, numDofDisp0, numDofDisp-numDofDisp0 ) =
            penalty_ * ( n1ng0 * n0ng1.transpose( ) ) + 
            penalty_ * diffAngle * ( curNormal0GradDof.transpose( ) * curNormal1GradDof );
    matDispDisp_.block( numDofDisp0, 0, numDofDisp-numDofDisp0, numDofDisp0 ) = 
            penalty_ * ( n0ng1 * n1ng0.transpose( ) ) + 
            penalty_ * diffAngle * ( curNormal1GradDof.transpose( ) * curNormal0GradDof );
    matDispDisp_.block( numDofDisp0, numDofDisp0, numDofDisp-numDofDisp0, numDofDisp-numDofDisp0 ) = 
            penalty_ * ( n0ng1 * n0ng1.transpose( ) );
   
    return;
}

//------------------------------------------------------------------------------
template< typename ELEMENT >
void gshell::fem::LinkRotationalInt< ELEMENT >::computeResidualAndMatrix( )
{
    this->computeResidual( );
    this->computeMatrix( );
    return;
}

//------------------------------------------------------------------------------
template< typename ELEMENT >
std::ostream & gshell::fem::LinkRotationalInt< ELEMENT >::write(
    std::ostream & os
    ) const
{
    os << "LinkRotationalInt" << std::endl;
    os << "    nodeId=" << node_->giveId() << std::endl;
    os << "    element0=" << element0_ << std::endl;
    os << "    element1=" << element1_ << std::endl;
    os << "    normal0=" << normal0_ << std::endl;
    os << "    normal1=" << normal1_ << std::endl;
    this->LinkPenaltyLMultiplier_::write( os );
    return os;
}

//------------------------------------------------------------------------------
template< typename ELEMENT >
void gshell::fem::LinkRotationalInt< ELEMENT >::computeNormal_(
    const unsigned elemLocId
    )
{
    Element * element = ( elemLocId == 0 ) ? element0_ : element1_;
    FTL_VERIFY( element );
    const unsigned numNodes = element->numFunctions( );
    for ( unsigned i = 0; i< numNodes; ++i )
        nodes_.push_back( element->giveSupportNodePtr( i ) );

    const unsigned i =  LinkRotationalInt_::nodeLocalIndex_( element, node_ );
    VecDim_ & normal = ( elemLocId == 0 ) ? normal0_ : normal1_;
    element->computeNormalAtNode( gshell::fem::REFERENCE, i, normal, NULL );
    return;
}
