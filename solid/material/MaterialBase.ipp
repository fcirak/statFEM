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

//! @file   MaterialBase.ipp
//! @author Fehmi Cirak, Thomas Rueberg
//! @date   2010

#include <corlib/verify.hpp>

//------------------------------------------------------------------------------
solid::material::MaterialBase::MaterialBase( const int numInternal, const bool planeStress) 
    : numInternal_( numInternal ), planeStress_( planeStress )
{
    parser_ = new corlib::PropertiesParser;
}

//------------------------------------------------------------------------------
solid::material::MaterialBase::~MaterialBase()
{
    FTL_VERIFY( parser_ != NULL );
    delete parser_;
}

//------------------------------------------------------------------------------
void solid::material::MaterialBase::readMaterialParams( std::istream& is ) 
{
    char c;
    is >> c;
    if( c == '{') {
        (is >> std::ws).get(c);
        while( c != '}') {
            is.putback(c);
	    
            corlib::skip_comment(is);
            (is >> std::ws).get(c); 
            if( c== '}') return;
	    
            is.putback(c);
	    
            parser_ -> readVariable( is );
	    
            (is >> std::ws).get(c); 
        }
    }		
    
    return;
}

//------------------------------------------------------------------------------
void solid::material::MaterialBase::printMaterialParams( std::ostream& os )
{
    os << '{' << "\n";
    parser_ -> printValues(os, " ", " ");
    os << " } \n";

    return;
}

