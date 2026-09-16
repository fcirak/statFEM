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

//! @author Burkhard Bornemann
//! @date   2010

//------------------------------------------------------------------------------
//! give the dofs increasing numbers
template< typename BNODE, unsigned DOF >
void gshell::fem::NodeStatic<BNODE,DOF>::numberDOFs( unsigned & counter )
{
    for ( unsigned i = 0; i < dof; i++ ) dofIndices_[ i ] = counter++;
    return;
}

//------------------------------------------------------------------------------
//! copy the dof array given a vector
template< typename BNODE, unsigned DOF >
void gshell::fem::NodeStatic<BNODE,DOF>::copyDofArray(
    std::vector< unsigned > & dofIndices
    ) const
{ 
    dofIndices.resize( dof );
    std::copy( dofIndices_.begin(), dofIndices_.end(), dofIndices.begin() );
    return;
}

//------------------------------------------------------------------------------
template< typename BNODE, unsigned DOF >
void gshell::fem::NodeStatic<BNODE,DOF>::updateDisplacements( ) 
{ 
    disp_ += increment_;
    increment_.setZero();
    return;
}

//------------------------------------------------------------------------------
//! store the constraints
template< typename BNODE, unsigned DOF >
void gshell::fem::NodeStatic<BNODE,DOF>::storeConstraint(
    const unsigned & component,
    const double & value
    )
{
    if ( component < dof ) 
        constraints_.push_back( std::make_pair( component, value ) );
    return;
}


//------------------------------------------------------------------------------
//! pass given Constraints back to a vector
template< typename BNODE, unsigned DOF >
void gshell::fem::NodeStatic<BNODE,DOF>::giveConstraints(
    ConstraintVec_ & globalConstraints
    ) const
{
    ConstraintVec_::const_iterator cIter = constraints_.begin();
    ConstraintVec_::const_iterator cEnd  = constraints_.end( );
    for( ; cIter != cEnd; cIter ++ ) {
        const unsigned index = dofIndices_[ cIter -> first ];
        globalConstraints.push_back( std::make_pair( index,
                                                     cIter -> second ) );
    }
    return;
}

//------------------------------------------------------------------------------
//! Overwrite values in displacement vector by prescribed values
template< typename BNODE, unsigned DOF >
void gshell::fem::NodeStatic<BNODE,DOF>::copyConstraintsToIncrements( )
{
    ConstraintVec_::const_iterator cIter = constraints_.begin();
    ConstraintVec_::const_iterator cEnd  = constraints_.end( );
    for( ; cIter != cEnd; ++cIter ) {
        increment_[ cIter->first ] = cIter->second;
    }
    return;
}

//------------------------------------------------------------------------------
//! Compute non-unit unit tangent
template< typename BNODE, unsigned DOF >
void gshell::fem::NodeStatic<BNODE,DOF>::sumNonUnitTangent_(
    const enum gshell::fem::config conf,
    const unsigned dir,
    VecDim & tangent
    ) const
{
    typedef typename BasisNode::MapNPtrDouble::const_iterator MapNPtrDoubleCIter;

    const MapNPtrDoubleCIter tcBegin = tangCoeffs_[ dir ].begin();
    const MapNPtrDoubleCIter tcEnd   = tangCoeffs_[ dir ].end();
    
    tangent.setZero();
    for ( MapNPtrDoubleCIter tcIter = tcBegin; tcIter != tcEnd; ++tcIter ) {
        const NodeStatic_ * nodePtr = static_cast<const NodeStatic_*>( tcIter->first );

        VecDim coord = nodePtr->giveCoordinates();
        if ( conf == gshell::fem::REFERENCE ) {
            ; // do nothing, already set
        }
        else if ( conf == gshell::fem::CURRENT ) {
            coord += nodePtr->giveDisplacementsRange( Eigen::VectorXi::LinSpaced(dim, 0, dim) );
        }
        else if ( conf == gshell::fem::LASTCONVERGED ) {
            coord += nodePtr->giveDisplacementsRange( Eigen::VectorXi::LinSpaced(dim, 0, dim) );
            coord -= nodePtr->giveIncrement( ).head( dim );
        }
        else {
            FTL_VERIFY_DESCRIPTIVE( false, "Not impl.\n" );
        }

        const double coeff = tcIter->second;
        tangent += coeff * coord;
    }

    return;
}

//------------------------------------------------------------------------------
//! Return unit tangent
template< typename BNODE, unsigned DOF >
typename gshell::fem::NodeStatic<BNODE,DOF>::VecDim 
gshell::fem::NodeStatic<BNODE,DOF>::getTangent(
    const enum gshell::fem::config conf,
    const unsigned dir
    ) const
{
    VecDim tang;  this->sumNonUnitTangent_( conf, dir, tang );
    const double length = tang.norm( );
    return ( tang / length );
}

//------------------------------------------------------------------------------
//! Return unit normal
template< typename BNODE, unsigned DOF >
typename gshell::fem::NodeStatic<BNODE,DOF>::VecDim 
gshell::fem::NodeStatic<BNODE,DOF>::getNormal(
    const enum gshell::fem::config conf
    ) const
{
    const VecDim tang0 = this->getTangent( conf, 0 );
    const VecDim tang1 = this->getTangent( conf, 1 );
    return corlib::cross_prod( tang0, tang1 );
}

//------------------------------------------------------------------------------
template< typename BNODE, unsigned DOF >
double gshell::fem::NodeStatic<BNODE,DOF>::getTangentCoefficient(
    const enum gshell::fem::config conf,
    NodeStatic_ * nodePtr,
    const unsigned dir
    ) const
{
    const double coeff = this->BasisNode::getTangentCoefficient( nodePtr, dir );
    VecDim tang;  this->sumNonUnitTangent_( conf, dir, tang );
    const double length = tang.norm( );
    return coeff / length;
}
