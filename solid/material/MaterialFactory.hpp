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

//! @file   MaterialFactory.hpp
//! @author Fehmi Cirak, Thomas Rueberg
//! @date   2010

#ifndef solid_material_materialfactory_h
#define solid_material_materialfactory_h
//------------------------------------------------------------------------------
#include <map>
#include <string>
//------------------------------------------------------------------------------

namespace solid{
    namespace material {
        class MaterialFactory;
    }
}

class solid::material::MaterialFactory 
{ 
public:
    typedef MaterialBase* (*CreateMaterialCallBack)();
    
private:
    typedef std::map<std::string, CreateMaterialCallBack> CallBackMap_;
    
public:
    bool registerMaterial( std::string materialType, CreateMaterialCallBack cb);    
    MaterialBase* createMaterial( std::string& materialType ) const;
    
    // this is a singleton class
    static MaterialFactory * instance();
    void destroy(); 

private:
    MaterialFactory(){};
    ~MaterialFactory(){};
    static MaterialFactory *instance_;
    
    // copy constructor and assignment operator
private:
    MaterialFactory ( const MaterialFactory & mf );
    const MaterialFactory & operator=( const MaterialFactory & mf );
    
private:
    CallBackMap_    callBacks_;
};

#include "MaterialFactory.ipp"
#endif
