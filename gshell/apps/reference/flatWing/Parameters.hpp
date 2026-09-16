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

#ifndef gshell_app_parameters_h
#define gshell_app_parameters_h

#include <corlib/PropertiesParser.hpp>

//==============================================================================
// declarations
namespace app{

    class Parameters;

}

//==============================================================================
// definitions

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

        prop->registerPropertiesVar( "constraintsFileName",    constraintsFileName );
        prop->registerPropertiesVar( "materialFileName",       materialFileName );
        prop->registerPropertiesVar( "monitorInputFileName",   monitorInputFileName );
        nodalForcesFileName = "";
        prop->registerPropertiesVar( "nodalForcesFileName",    nodalForcesFileName );

        prop->registerPropertiesVar( "animationFolder",        animationFolder );
        prop->registerPropertiesVar( "monitorOutputFileName",  monitorOutputFileName );
        prop->registerPropertiesVar( "animationBaseFileName",  animationBaseFileName );

        prop->registerPropertiesVar( "thickness",              thickness );

        newmarkBeta = 0.25;
        prop->registerPropertiesVar( "newmarkBeta",            newmarkBeta );
        newmarkGamma = 0.5;
        prop->registerPropertiesVar( "newmarkGamma",           newmarkGamma );
        prop->registerPropertiesVar( "timeStepMax",            timeStepMax );
        prop->registerPropertiesVar( "timeStepSize",           timeStepSize );
        prop->registerPropertiesVar( "writeEveryStep",         writeEveryStep );
    
        // read variables from the input.dat file
        prop->readValues( inp );
        delete prop;

    }

public: // public on purpose!
    // input file names
    std::string inputMeshFileName;
    std::string secondaryInputFileName;

    // output file names
    std::string outputFolder;
    std::string outputBaseFileName;

    // mechanical shell
    std::string constraintsFileName;
    std::string materialFileName;
    std::string monitorInputFileName;
    std::string nodalForcesFileName;
    std::string animationFolder;
    std::string monitorOutputFileName;
    std::string animationBaseFileName;

    // geometry
    double thickness;

    // BC & load
    std::string dirichletBC;

    // time integration
    double newmarkBeta;
    double newmarkGamma;
    unsigned timeStepMax;
    double timeStepSize;
    unsigned writeEveryStep;

};


#endif
