// -*- C++ -*- 
//
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//
//                                   Fehmi Cirak
//                        California Institute of Technology
//                           (C) 2005 All Rights Reserved
//
// <LicenseText>
//
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//
#include <utility>

#include <corlib/verify.hpp>

//------------------------------------------------------------------------------
solid::material::MaterialFactory * solid::material::MaterialFactory::instance_ = NULL;

//------------------------------------------------------------------------------
solid::material::MaterialFactory * solid::material::MaterialFactory::instance( ) 
{    
    if ( !instance_ ) {
        instance_ = new MaterialFactory( );
    }
    
    return instance_;
}

//------------------------------------------------------------------------------
void solid::material::MaterialFactory::destroy()
{
    callBacks_.clear();
    delete instance_;
}


bool solid::material::MaterialFactory::registerMaterial( std::string materialType, 
                                                         CreateMaterialCallBack cb ) 
{	
    return callBacks_.insert( CallBackMap_::value_type( materialType, 
                                                        cb ) ).second; 
}

//------------------------------------------------------------------------------
solid::material::MaterialBase * 
solid::material::MaterialFactory::createMaterial( std::string& materialType ) const 
{
    CallBackMap_::const_iterator it = callBacks_.find( materialType );
    
    if ( it == callBacks_.end( ) ) {
        FTL_VERIFY_DESCRIPTIVE( false, "Unknown material type" );
    }
    
    return ( it->second( ) );	
}
