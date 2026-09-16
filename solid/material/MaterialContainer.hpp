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

//! @file   MaterialContainer.hpp
//! @author Fehmi Cirak, Thomas Rueberg
//! @date   2010

#ifndef solid_material_materialcontainer_h
#define solid_material_materialcontainer_h
//------------------------------------------------------------------------------
#include <vector>
#include <iosfwd>

#include <corlib/misc.hpp>
#include <corlib/verify.hpp>
//------------------------------------------------------------------------------
namespace solid{
    namespace material{
        class MaterialContainer;
        class MaterialBase;
    }
}

class solid::material::MaterialContainer 
{
public:
    // this is a singleton class
    static MaterialContainer * instance();
    void destroy();
    
    void readMaterialStream(  std::istream &is );
    void printMaterialStream( std::ostream &os );

    void addMaterial( MaterialBase *mat ) 
    {
        materials_.push_back(mat);
    }

    // acessor
    MaterialBase* getMaterial( const unsigned i ) const 
    {
        FTL_VERIFY( i < materials_.size() );
        return materials_[i];
    }
    
    // this is a singleton class
private:
    MaterialContainer(  ){};
    ~MaterialContainer( );

    static MaterialContainer *instance_;

    friend class corlib::SingletonDestroyer< MaterialContainer >;
    static corlib::SingletonDestroyer< MaterialContainer > singletonDestroyer_;

    // copy constructor and assignment operator
private:
    MaterialContainer ( const MaterialContainer & mc );
    const MaterialContainer & operator=( const MaterialContainer & mc );

private:
    typedef std::vector<MaterialBase *>   MaterialCont_;
    typedef MaterialCont_::iterator        MaterialIt_;

    MaterialCont_                          materials_;
};

#include "MaterialContainer.ipp"
#endif
