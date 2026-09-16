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

//! @file   MaterialBase.hpp
//! @author Fehmi Cirak, Thomas Rueberg
//! @date   2010

#ifndef solid_material_materialbase_h
#define solid_material_materialbase_h

#include <iosfwd>
#include <array>

#include <corlib/PropertiesParser.hpp>
#include <Eigen/Core>
#include <corlib/eigenX.hpp>

//------------------------------------------------------------------------------
namespace solid {
    namespace material {
        class MaterialBase;

        namespace eigenX = corlib::eigenX;
    }
}

//------------------------------------------------------------------------------
/** \brief Abstract base class for material types
 */
class solid::material::MaterialBase 
{
public:
    MaterialBase( const int numInternal, const bool planeStress );
    virtual ~MaterialBase(); 
    
protected:
    typedef eigenX::MatrixSd< 3, 3 >          Mat3x3_;
    typedef Eigen::Array< Mat3x3_, 3, 3 >     Mat3x3x3x3_;

public:
    /// compute energy, second Piola-Kirchhoff stress and elasticity tensor
    virtual double strainEnergy( const Mat3x3_ & F ) const = 0;
    virtual void FPKstress( const Mat3x3_ & F, Mat3x3_ & P ) const = 0;
    virtual void elasticityTensor( const Mat3x3_ & F, Mat3x3x3x3_ & C ) const = 0;

    /// return value of mass density
    virtual double density() = 0;    

    /// return the estimated velocity of the pressure wave
    virtual double estimatedWaveVelocity() const = 0;

    /// estimated Young's modulus (used for computing stable time step)
    virtual double estimatedElasticityConst( ) = 0;

    /// initialize history vars. and perform material specific initialization
    virtual void initialization( double * history ) = 0; 

    /// @name I/O 
    //@{
    void readMaterialParams(  std::istream& is );
    void printMaterialParams( std::ostream& os );
    //@}

    /// @name accessors
    //@{
    int numInternal(  ) const { return numInternal_;}
    bool planeStress( ) const { return planeStress_;}
    //@}
    
protected:
    template <typename T>
    void registerVariable( const std::string& name, T & variable ) {
        parser_ -> registerPropertiesVar( name, variable );
    }
    
private:
    const int      numInternal_;
    const bool     planeStress_;
    
    corlib::PropertiesParser   * parser_;
};

#include "MaterialBase.ipp"
#endif

