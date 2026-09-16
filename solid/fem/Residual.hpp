// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//
//                 Computational Structural Mechanics Lab
//                         University of Cambridge
//
// This software is copyrighted and all rights are retained by the University of
// Cambridge. All use, disclosure, and/or reproduction of any part not expressly
// authorised by the University of Cambridge is prohibited. (C) 2010
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

//! @file Residual.hpp

#ifndef solid_fem_residual_h
#define solid_fem_residual_h
//------------------------------------------------------------------------------
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_sparse.hpp>
#include <corlib/VectorAssembler.hpp>
#include <corlib/Distributor.hpp>
#include <corlib/Integrator.hpp>

//------------------------------------------------------------------------------
namespace solid{
    namespace fem{
        template<typename MESH,typename QUAD> class Residual;
        namespace ublas = boost::numeric::ublas;
    }
}

//------------------------------------------------------------------------------
template<typename MESH,typename QUAD> 
class solid::fem::Residual
{
public:
    typedef ublas::vector<double>                 DVec_;
    typedef ublas::matrix<double>                 DMat_;
    typedef ublas::compressed_matrix<double>      SparseMat_;
    typedef std::vector<unsigned>                 DofVec_;
    typedef solid::fem::Residual<MESH,QUAD>       Me_;

    //! Constructor with size, mesh pointer and constraint factor
    Residual( MESH * m );

public:
    //! pass constraints to caller
    void fillConstraintVec( std::vector< std::pair<unsigned,double> > & constraints );

    //! Distribute the computed solution
    void distributeSolution( DVec_ & X );

    //! Give the solution values to a container 
    void giveSolution( const std::vector<unsigned> & dofIndices, 
                       DVec_ & result ) const;

public:
    //! Compute the residual forces at each node
    void computeNodalResiduals( );

    //! Collect the computed residual
    void collectResidual( DVec_ & F );

    //! Let the residual be filled with forces
    void insertToRhs( const DVec_ & force, const DofVec_ & dofIndices );

protected:
    MESH *         mesh_;                  //!< access to the mesh
    QUAD           quadrature_;            //!< quadrature object
    DVec_   *      increment_;             //!< current approximate solution
    DVec_   *      residual_;              //!< vector of residual forces
};

//------------------------------------------------------------------------------
#include "Residual.ipp"
//------------------------------------------------------------------------------
#endif
