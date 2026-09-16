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

//! @file   MaterialContainer.ipp
//! @author Fehmi Cirak, Thomas Rueberg
//! @date   2010

#include <solid/material/MaterialBase.hpp>
#include <solid/material/MaterialFactory.hpp>
#include <solid/material/svenant/SVenant.hpp>
#include <solid/material/neohookean/Neohookean.hpp>

#include <corlib/misc.hpp>
#include <corlib/PropertiesParser.hpp>

#include <algorithm>
#include <iostream>
#include <typeinfo>
#include <functional>


//------------------------------------------------------------------------------
//! static pointer to itself
solid::material::MaterialContainer * 
solid::material::MaterialContainer::instance_ = NULL;

//------------------------------------------------------------------------------
//! static singleton destroyer
corlib::SingletonDestroyer<solid::material::MaterialContainer> 
solid::material::MaterialContainer::singletonDestroyer_;

//------------------------------------------------------------------------------
//! create an instance 
solid::material::MaterialContainer * 
solid::material::MaterialContainer::instance( ) 
{    
    if (!instance_) {
        instance_ = new MaterialContainer( );

        // register the createMaterial functions with the MaterialFactory
        MaterialFactory * mf = MaterialFactory::instance( );
	
        mf -> registerMaterial( "SVenant",    SVenant::createMaterial );
        mf -> registerMaterial( "Neohookean", Neohookean::createMaterial );

        singletonDestroyer_.setSingleton( instance_ );
    }
    
    return instance_;
}

//------------------------------------------------------------------------------
//! delete all materials in the container
void solid::material::MaterialContainer::destroy( )
{
    singletonDestroyer_.releaseSingleton();
    delete instance_;
}

//------------------------------------------------------------------------------
//! Destructor of materials and the factory
solid::material::MaterialContainer::~MaterialContainer( )
{
    std::for_each( materials_.begin( ), materials_.end( ), 
                   corlib::deleteFunctor( ) );
    materials_.clear();
    
    MaterialFactory * mf = MaterialFactory::instance( );
    mf -> destroy( );
}

//------------------------------------------------------------------------------
//! read materials from input stream
void solid::material::MaterialContainer::readMaterialStream( std::istream& is )
{
    MaterialFactory * mf = MaterialFactory::instance( );

    while ( is.good () ) {
        corlib::skip_comment( is );
        std::string materialType;
        is >> materialType;
	
        MaterialBase * newmat = mf -> createMaterial( materialType );
        newmat -> readMaterialParams( is );	
        materials_.push_back( newmat );
        corlib::skip_comment( is );
    }

    // initialize all materials
    //std::for_each(_materials.begin(),
    //		  _materials.end(),
    //	  std::mem_fun(&MaterialBase::initialization));
    
    return;
}

//------------------------------------------------------------------------------
//! print materials to output stream
void solid::material::MaterialContainer::printMaterialStream( std::ostream& os )
{
    MaterialIt_ it  = materials_.begin();
    MaterialIt_ ite = materials_.end();

    FTL_VERIFY( it != ite );

    while ( it != ite ) {
        os << typeid( (*it) ).name() << "\n";
        (*it) -> printMaterialParams( os );
        ++it;    
    }
    
    return;
}
