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
//! @date   05/2012

#ifndef gshell_fem_coupler_h
#define gshell_fem_coupler_h

#include <numeric>
#include <array>

#include <corlib/DistributeQuantity.hpp>
#include <corlib/CollectQuantity.hpp>
#include <gshell/fem/misc.hpp>

// TODO: Activate the following once fsi module is ported. The fsi module is
//       used in Coupler.ipp as fsi::utils::BodyForceDiscreteConstant
//#include <fsi/utils/BodyForceDiscrete.hpp>

//==============================================================================
namespace gshell {
    namespace fem {

        template< typename DRIVER >
        class Coupler;

    }
}

//==============================================================================
/// Coupler interface tools for shell drivers providing 'unified' methods
/// for coupling with other physical fields (such as fluids for instance)
///
/// \tparam DRIVER    Shell driver type
template< typename DRIVER >
class gshell::fem::Coupler
{
public:
    typedef DRIVER                                   Driver;
                                                    
    typedef typename Driver::Mesh                    Mesh;
    typedef typename Driver::Node                    Node;
    typedef typename Driver::Element                 Element;
                                                    
    typedef typename Node::VecDim                    VecDim;
    typedef std::array<unsigned,Element::numNodes>   ArrayNN;

public:
    /// Predict nodal displacement increments
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

    /// Copy current coordinates of shell to surface
    static
    void getSurfaceCoordinates( const Mesh & mesh,
                                std::vector<VecDim> & nodeCoordinates );

    /// Copy element connectivity from shell to surface
    static
    void getSurfaceConnectivity( const Mesh & mesh,
                                 std::vector<ArrayNN> & elementConnectivities );

    /// Current velocities of beam to surface with inverted-backward-Euler
    static
    void getSurfaceVelocitiesFD( Mesh & mesh,
                                 const double stepSize,
                                 std::vector<VecDim>  & nodeVelocities );
    /// Copy current coordinates of shell to surface
    static
    void getSurfaceVelocities( Driver & driver,
                               const double stepSize, 
                               std::vector<VecDim>  & nodeVelocities );

    /// Copy current coordinates of shell to surface
    static
    void getSurfaceNormals( const Mesh & mesh,
                            std::vector<VecDim>  & elementCentres,
                            std::vector<VecDim>  & elementNormals );

    /// Copy force on surface nodes to shell nodes
    static
    void setSurfaceForceExt( Driver & driver,
                             const std::vector<VecDim>  & elementTractions );

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
                             std::vector<VecDim>  & nodeIncrementsPrev,
                             const double relax );

};

#include "Coupler.ipp"

#endif
