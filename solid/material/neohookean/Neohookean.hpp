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

//! @file   Neohookean.ipp
//! @author Thomas Rueberg
//! @date   2010

#ifndef solid_material_neohookean_h
#define solid_material_neohookean_h
#include <solid/material/MaterialBase.hpp>

#include <Eigen/Core>
#include <corlib/eigenX.hpp>

namespace solid {
    namespace material {
        class Neohookean;

        namespace eigenX = corlib::eigenX;
    }
}

//------------------------------------------------------------------------------
class solid::material::Neohookean : 
    public solid::material::MaterialBase 
{
private:
    Neohookean( ); // use createMaterial instead

    typedef eigenX::MatrixSd< 3, 3 >          Mat3x3_;
    typedef Eigen::Array< Mat3x3_, 3, 3 >     Mat3x3x3x3_;

public:
    ~Neohookean();
    
    double density();    
    double estimatedElasticityConst();
    double estimatedWaveVelocity() const;
    void initialization(double * history);
    double strainEnergy( const Mat3x3_ & F ) const;
    void FPKstress( const Mat3x3_ & F, Mat3x3_ & P ) const;
    void elasticityTensor( const Mat3x3_ & F, Mat3x3x3x3_ & C ) const;

    static MaterialBase * createMaterial() 
    {
        return new Neohookean();
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
inline double solid::material::Neohookean::density()
{
    return massDensity_;
}

inline double solid::material::Neohookean::estimatedElasticityConst()
{
    return youngModulus_;
}

inline void solid::material::Neohookean::initialization(double * history) 
{
    return;
}

inline double solid::material::Neohookean::estimatedWaveVelocity() const
{
    const double nu = poissonRatio_;
    const double  E = youngModulus_;
    return std::sqrt( E / massDensity_ *(1.-nu)/(1.-2*nu)/(1.+nu) );
}

#include "Neohookean.ipp"
#endif 
