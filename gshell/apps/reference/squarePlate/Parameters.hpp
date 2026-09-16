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

//! @author Kosala Bandara
//! @date   2010

#ifndef gshell_apps_parameters_h
#define gshell_apps_parameters_h

#include <corlib/PropertiesParser.hpp>

//------------------------------------------------------------------------------
// forward declarations
namespace app{

    class Parameters;

}

//------------------------------------------------------------------------------
//! Storage of parameters related to the computation
class app::Parameters
{
public:
    /// Constructor
    Parameters( std::istream & inp )
    {
        corlib::PropertiesParser * prop = new corlib::PropertiesParser;

        prop->registerPropertiesVar( "inputMeshFileName",      inputMeshFileName );
        prop->registerPropertiesVar( "secondaryInputFileName", secondaryInputFileName );

        prop->registerPropertiesVar( "outputFolder",           outputFolder );
        prop->registerPropertiesVar( "outputBaseFileName",     outputBaseFileName );

        prop->registerPropertiesVar( "thicknessSF",            thicknessSF );
        prop->registerPropertiesVar( "thicknessKL",            thicknessKL );

        prop->registerPropertiesVar( "materialFileName",       materialFileName );
        prop->registerPropertiesVar( "monitorInputFileName",   monitorInputFileName );

        prop->registerPropertiesVar( "animationFolder",        animationFolder );
        prop->registerPropertiesVar( "monitorOutputFileName",  monitorOutputFileName );
        prop->registerPropertiesVar( "animationBaseFileName",  animationBaseFileName );
    
        // read variables from the input.dat file
        prop->readValues( inp );
        delete prop;

    }

public: // on purpose!
    // input file names
    std::string inputMeshFileName;
    std::string secondaryInputFileName;

    // output file names
    std::string outputFolder;
    std::string outputBaseFileName;

    // mechanical shell
    double      thicknessSF;
    double      thicknessKL;
    std::string materialFileName;
    std::string monitorInputFileName;

    std::string animationFolder;
    std::string monitorOutputFileName;
    std::string animationBaseFileName;

};


#endif
