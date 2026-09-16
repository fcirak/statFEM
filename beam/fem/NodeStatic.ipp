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
template<unsigned DIM, unsigned DOF>
void beam::fem::NodeStatic< DIM, DOF> :: numberDOFs(unsigned & counter)
{
    for ( unsigned i = 0; i < DOF; i ++ ) dofIndices_[ i ] = counter ++;
    return;
}

//------------------------------------------------------------------------------
//! copy the dof array given an iterator
template< unsigned DIM, unsigned DOF >
void beam::fem::NodeStatic< DIM, DOF >::
copyDofArray( std::vector<unsigned> & dofIndices ) const
{ 
    dofIndices.resize( dofIndices_.size() );
    std::copy( dofIndices_.begin( ), dofIndices_.end( ), dofIndices.begin() );
    return;
}

//------------------------------------------------------------------------------
template< unsigned DIM, unsigned DOF >
void beam::fem::NodeStatic< DIM, DOF > :: updateDisplacements( ) 
{ 
    disp_ += increment_;
    increment_.setZero();
    return;
}

