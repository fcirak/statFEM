// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//
//                 Computational Structural Mechanics Lab
//                         University of Cambridge
//
// This software is copyrighted and all rights are retained by the University of
// Cambridge. All use, disclosure, and/or reproduction of any part not expressly
// authorised by the University of Cambridge is prohibited. (C) 2010
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

//! @file CouplerBeam.hpp

#ifndef beam_fem_couplerbeam_h
#define beam_fem_couplerbeam_h

#include <array>
#include <numeric>

#include <corlib/DistributeQuantity.hpp>
#include <corlib/CollectQuantity.hpp>

#include <fsi/utils/BodyForceDiscrete.hpp>
#include <beam/fem/misc.hpp>

//==============================================================================
namespace beam {
    namespace fem {
        namespace eigenX = corlib::eigenX;

        template<typename DRIVER> class CouplerBeam;
    }
}

//==============================================================================
/// CouplerBeam interface to beam drivers providing 'unified' methods
/// for coupling with other physical fields (such as fluids for instance)
///
/// \tparam DRIVER    Beam driver type
template< typename DRIVER >
class beam::fem::CouplerBeam
{
public:
    typedef DRIVER                                   Driver;

    typedef typename Driver::Mesh                    Mesh;
    typedef typename Driver::Node                    Node;
    typedef typename Driver::Element                 Element;

    typedef typename Node::VecDim                    VecDim;
    typedef std::array<unsigned,Element::numNodes>   ArrayNN;

public:
    /// Predict nodal displacement increments with Forward-Euler approach
    static
    void predictIncrementsFE( Mesh & mesh,
                              const double stepSize,
                              std::vector<VecDim> & nodeIncrements );
    /// Predict nodal displacement increments using last increment
    static
    void predictIncrementsLI( Mesh & mesh,
                              const std::vector<VecDim> & nodeIncrements );

    /// Predict nodal displacement increments to be zero
    static
    void predictIncrementsZero( Mesh & mesh,
                                std::vector<VecDim> & nodeIncrements );

    /// Copy current coordinates of beam to surface
    static
    void getSurfaceCoordinates( Mesh & mesh,
                                std::vector<VecDim> & nodeCoordinates );

    /// Copy element connectivity from beam to surface
    static
    void getSurfaceConnectivity( Mesh & mesh,
                                 std::vector<ArrayNN> & elementConnectivities );

    /// Current velocities of beam to surface with a finite difference
    static
    void getSurfaceVelocitiesFD( Mesh & mesh,
                                 const double stepSize,
                                 std::vector<VecDim> & nodeVelocities );

    /// Current velocities of beam to surface
    static
    void getSurfaceVelocities( Driver & driver, 
                               const double stepSize,
                               std::vector<VecDim> & nodeVelocities );

    /// Copy current coordinates of beam to surface
    static
    void getSurfaceNormals( Mesh & mesh,
                            std::vector<VecDim> & elementCentres,
                            std::vector<VecDim> & elementNormals,
                            const bool onePerElement = true );

    /// Copy force on surface nodes to beam nodes
    static
    void setSurfaceForceExt( Driver & driver,
                             const std::vector<VecDim> & elementTractions,
                             const bool onePerElement = true );

    //! Simply compute the residual norm
    static 
    double incrementResidual( Mesh & mesh,
                              const std::vector<VecDim> & nodeIncrementsPrev );

    /// Scale increments by #relax factor
    ///
    /// \param[in,out]  nodeIncrementsPrev    Nodal increments of previous cycle
    /// \param[in]      relax                 Relaxation parameter to scale with
    /// \return                               Norm of difference in increments
    static
    double updateIncrements( Mesh & mesh,
                             std::vector<VecDim> & nodeIncrementsPrev,
                             const double relax );

};

#include "CouplerBeam.ipp"

#endif
