// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//
//                 Computational Structural Mechanics Lab
//                         University of Cambridge
//
// This software is copyrighted and all rights are retained by the University of
// Cambridge. All use, disclosure, and/or reproduction of any part not expressly
// authorised by the University of Cambridge is prohibited. (C) 2010
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

//! @file MaterialParameters.hpp

#ifndef del2_fem_materialparameters_h
#define del2_fem_materialparameters_h

#include <corlib/PropertiesParser.hpp>

//==============================================================================
// forward declarations
namespace del2{
    namespace fem{
        
        class MaterialParameters;

    }
}


//==============================================================================
//! Storage of material parameters related to the computation
class del2::fem::MaterialParameters
{
public:
    /// Constructor
    MaterialParameters( std::istream & inp ) :
        conductivity_( 0.0 ), capacity_( 0.0 )
    {
        corlib::PropertiesParser * prop = new corlib::PropertiesParser;

        prop->registerPropertiesVar( "conductivity",           conductivity_ );
        prop->registerPropertiesVar( "capacity",               capacity_ );
    
        // read variables from the material.dat file
        prop->readValues( inp );
        delete prop;
    }

    double getConductivity() const { return conductivity_; }
    double getCapacity() const { return capacity_; }

protected:
    // the material parameters
    double conductivity_;
    double capacity_;
};


#endif
