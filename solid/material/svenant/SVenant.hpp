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

//! @file   SVenant.hpp
//! @author Thomas Rueberg
//! @date   2010

#ifndef solid_material_svenant_h
#define solid_material_svenant_h
#include <solid/material/MaterialBase.hpp>

namespace solid {
    namespace material {
        class SVenant;
    }
}

//------------------------------------------------------------------------------
class solid::material::SVenant : public solid::material::MaterialBase 
{
private:
    SVenant( ); // use createMaterial instead

public:
    ~SVenant();
    
    double density();    
    double estimatedElasticityConst();
    double estimatedPoissonRatio();
    double estimatedWaveVelocity() const;
    void initialization(double * history);
    double strainEnergy( const Mat3x3_ & F ) const;
    void FPKstress( const Mat3x3_ & F, Mat3x3_ & P ) const;
    void elasticityTensor( const Mat3x3_ & F, Mat3x3x3x3_ & P ) const;

    static MaterialBase * createMaterial() 
    {
        return new SVenant();
    }

private:
    void computeLame( double & lambda, double & mu ) const;


    static const int numberInternalVariables_ = 1;

    double     youngModulus_;
    double     poissonRatio_;
    double     massDensity_;
};

//------------------------------------------------------------------------------
//! inlined functions
inline double solid::material::SVenant::density()
{
    return massDensity_;
}

inline double solid::material::SVenant::estimatedElasticityConst()
{
    return youngModulus_;
}

inline double solid::material::SVenant::estimatedPoissonRatio()
{
    return poissonRatio_;
}

inline void solid::material::SVenant::initialization(double * history) 
{
    return;
}

inline double solid::material::SVenant::estimatedWaveVelocity() const
{
    const double nu = poissonRatio_;
    const double  E = youngModulus_;
    // longitudinal wave velocity in infinite solid
    return std::sqrt( E / massDensity_ *(1.-nu)/(1.-2*nu)/(1.+nu) );
}

#include "SVenant.ipp"

#endif
