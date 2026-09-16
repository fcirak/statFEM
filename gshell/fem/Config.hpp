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
//

//! @author Burkhard Bornemann, Kosala Bandara
//! @date   2010

#ifndef gshell_fem_config_h
#define gshell_fem_config_h

//==============================================================================
// declarations
namespace gshell {
    namespace fem {

        //----------------------------------------------------------------------
        //! Configuration types
        enum config {
            REFERENCE,       //!< reference configuration
            LASTCONVERGED,   //!< last converged configuraion
            CURRENT          //!< current configuration
        };
        
    }
}

#endif
