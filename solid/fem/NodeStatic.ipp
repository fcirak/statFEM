// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//
//                 Computational Structural Mechanics Lab
//                         University of Cambridge
//
// This software is copyrighted and all rights are retained by the University of
// Cambridge. All use, disclosure, and/or reproduction of any part not expressly
// authorised by the University of Cambridge is prohibited. (C) 2010
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

//! @file NodeStatic.ipp

//------------------------------------------------------------------------------
//! give the dofs increasing numbers
template<typename BNODE, unsigned DOF>
void solid::fem::NodeStatic<BNODE, DOF>::numberDOFs(unsigned & counter)
{
    for ( unsigned i = 0; i < DOF; i ++ ) dofIndices_[ i ] = counter ++;
    return;
}

//------------------------------------------------------------------------------
//! copy the dof array to another vector
template<typename BNODE, unsigned DOF>
void solid::fem::NodeStatic<BNODE, DOF>::
copyDofArray( std::vector<unsigned> & out ) const
{ 
    out.resize( dof );
    std::copy( dofIndices_.begin( ), dofIndices_.end( ), out.begin() );
    return;
}

//------------------------------------------------------------------------------
//! Add increment to displacement field, clear the increment
template<typename BNODE, unsigned DOF>
void solid::fem::NodeStatic<BNODE, DOF>::updateDisplacements( ) 
{ 
    disp_ += increment_;
    increment_.clear();
    return;
}

//------------------------------------------------------------------------------
//! Store a constraint
template<typename BNODE, unsigned DOF>
void solid::fem::NodeStatic<BNODE, DOF>::
storeConstraint(const unsigned & component, const double & value)
{
    if (component < dof) 
        constraints_.push_back( std::make_pair( component, value ) );
    return;
}

//------------------------------------------------------------------------------
//! Pass given constraints back to a vector using the global dof indices
template<typename BNODE, unsigned DOF>
void solid::fem::NodeStatic<BNODE, DOF>::
giveConstraints( ConstraintVec_ & globalConstraints ) const
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






