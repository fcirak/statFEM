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

//! @file   PlaneStressDecorator.hpp
//! @author Thomas Rueberg
//! @date   2010

#ifndef solid_material_planestressdecorator_h
#define solid_material_planestressdecorator_h
//------------------------------------------------------------------------------
#include <solid/material/MaterialBase.hpp>
#include <Eigen/Core>
#include <corlib/eigenX.hpp>

namespace solid{
    namespace material{
        class PlaneStressDecorator;
        
        namespace eigenX = corlib::eigenX;
    }
}


//------------------------------------------------------------------------------
/** \brief Wraps around material calls in order to achieve state of plane stress
 *
 * Masks the call to a material by computing F_33 such that P_33=0. Note that
 * the 3-direction is fixed here, whereas in future a vector normal to the plane
 * should be used (probably defaulting to (0,0,1). 
 * This class is inherited from MaterialBase in order to enforce a conforming 
 * interface.
 */
class solid::material::PlaneStressDecorator : public solid::material::MaterialBase
{
private:
    typedef solid::material::MaterialBase          MaterialBase_;

public:
    typedef eigenX::MatrixSd< 3, 3 >          Mat3x3_;
    typedef Eigen::Array< Mat3x3_, 3, 3 >     Mat3x3x3x3_;

    PlaneStressDecorator( MaterialBase * mat ) 
        : MaterialBase( mat->numInternal(), mat->planeStress() ), material_( mat ) { }

    // compute energy, second Piola-Kirchhoff stress and elasticity tensor
    virtual double strainEnergy( const Mat3x3_ & F ) const;
    virtual void FPKstress( const Mat3x3_ & F, Mat3x3_ & P ) const;
    virtual void elasticityTensor( const Mat3x3_ & F, Mat3x3x3x3_ & C ) const;


    // return value of mass density
    virtual double density() { return material_->density(); }    

    //  return the estimated velocity of the pressure wave
    virtual double estimatedWaveVelocity() const { return material_->estimatedWaveVelocity(); }

    // estimated young's modulus (used for computing stable time step)
    virtual double estimatedElasticityConst( ) { return material_->estimatedElasticityConst(); }     
    
    // initialize history vars. and perform material specific initialization
    virtual void initialization( double * history ){ return material_->initialization( history ); }


private:
    MaterialBase * material_; //<! access to the real material

    //! Newton method for P_33( F ) = 0
    void iterateUntilPlaneStress( Mat3x3_ & F ) const;
};

#include "PlaneStressDecorator.ipp"

#endif
