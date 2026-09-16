// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//
//                 Computational Structural Mechanics Lab
//                         University of Cambridge
//
// This software is copyrighted and all rights are retained by the University of
// Cambridge. All use, disclosure, and/or reproduction of any part not expressly
// authorised by the University of Cambridge is prohibited. (C) 2010 --
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

//! @file CovarianceForcing.hpp

#ifndef statX_infer_covarianceforcing_h
#define statX_infer_covarianceforcing_h

#include <vector>
#include <string>
#include <map>

#include <Eigen/Core>

#include <corlib/verify.hpp>

//------------------------------------------------------------------------------
namespace statX {
    namespace infer {
        template <unsigned DOF, typename KERNEL>
        class CovarianceForcing;
    }
}

//------------------------------------------------------------------------------
/** \brief Approximate covariance computed from unit load vector
 *  \details
 *
 *  Computes the approximate covariance matrix \f$ C_f \f$ and
 *  \f[
 *    (C_f)_{ij} \approx \left (\int_\Omega \phi_i(\vec x) d\Omega (\vec x) \right )
 *                       c_f(\vec x_i , \, \vec x_j)
 *                     \left (\int_\Omega \phi_j(\vec x') d\Omega (\vec x') \right )
 *  \f]
 *  gives the projected projected matrix \f$ \vec P \vec C_f \vec P^T \f$
 *
 *  \tparam DOF    degress of freedom per node
 *  \tparam KERNEL a kernel like the exponential kernel ExpKernel
 */

template <unsigned DOF, typename KERNEL>
class statX::infer::CovarianceForcing
{
public:
    //! template parameters and values
    typedef KERNEL Kernel;
    static constexpr unsigned dof = DOF;

    //! public convenience typedefs
    typedef std::map<std::string, double>      HpsValueMap;
    typedef typename Kernel::VecDim            VecDim;
    typedef typename Eigen::VectorXd           VectorLocal;
    typedef std::vector<unsigned>              VectorUInt;

private:
    //! private convenience typedefs
    typedef std::vector<VecDim> CoordVec_;
    static constexpr auto dim_ = VecDim::RowsAtCompileTime;
    typedef std::vector< std::pair<unsigned,double> >     VectorConstraint_;

public:
    //! constructor with number of nodes and cut-off threshold
    CovarianceForcing ( const unsigned numNodes, const double cutEps = 0.0 )
            : numNodes_( numNodes ), cutEps_( cutEps )
    {
        covVec_.resize( numNodes * dof );
        covVec_.setZero( );
        coords_.reserve( numNodes * dof );
        kernel_.setHypIdPre( "forcing." );
    }

    //! set preamble for hyperparameter names
    void setHypIdPre ( const std::string & preamble )
    {
        hypIdPre_ = preamble;
    }

    //! insert element RHS (eCovVec) with unit distributed load to global RHS (covVec_)
    void insertToRhs ( const VectorLocal & eCovVec, const VectorUInt & dofIndices );

    //! insert nodal coordinates (nCoords) one by one into the global vector coords_
    //! node numbering is assumed to be consistent like for the RHS
    void insertToCoordinates ( const VecDim & nCoords );


    //! give the projected RHS covariance \f$ \vec P \vec C_u \vec P^T \f$ for
    //! given hyperparameters and projection matrix P
    void giveProjectedCovar ( const HpsValueMap & params,
                              const Eigen::MatrixXd & pMat,
                              Eigen::MatrixXd & projCovmat ) const;

    //! set covVec to zero at the constrained nodes
    //! to enforce that the solution has zero variance at those nodes
    void applyConstraints( const VectorConstraint_ & constraints,
                           const double factor,
                           const double scalar )
    {
        VectorConstraint_::const_iterator cIter = constraints.begin( );
        VectorConstraint_::const_iterator cEnd  = constraints.end( );

        for ( ; cIter != cEnd; ++ cIter ) {
            const unsigned index = cIter -> first;
            const double   value = factor * (cIter -> second);
            bcIndices_.push_back( index );
            covVec_[ index ] = 0.;
        }
    }

    //! get the indices of the constrained nodes used to impose zero variance values
    void getConstrainedIndices ( VectorUInt & bcIndices )
    {
        bcIndices = bcIndices_;
        return;
    }

private:
    //! number of nodes
    const unsigned   numNodes_;

    //! cut-off threshold
    const double     cutEps_;

    //! preamble for local parameter names
    std::string      hypIdPre_;

    //! covariance kernel
    Kernel           kernel_;

    //! element RHS vector (unit distributed load) and nodal coordinates
    VectorLocal      covVec_;
    CoordVec_        coords_;
    VectorUInt       bcIndices_;

};

//------------------------------------------------------------------------------
#include "CovarianceForcing.ipp"

#endif
