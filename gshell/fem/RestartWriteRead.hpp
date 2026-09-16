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

//! @author Burkhard Bornemann
//! @date   04/2012

// Note:
// Functions originate from blox/apps/user/fsiTest/Restart.hpp

#ifndef gshell_fem_restartwriteread_h
#define gshell_fem_restartwriteread_h

#include <stdio.h>

//==============================================================================
namespace gshell{
    namespace fem{

        template<typename DRIVER>
        void restartWrite( const std::string & basename,
                           const unsigned step,
                           DRIVER & driver );

        template< typename GETFUN, typename DRIVER >
        void restartWriteData( std::ostream & out, GETFUN & getFun,
                               DRIVER & driver );

        template<typename DRIVER>
        void restartRead( const std::string & basename,
                          const unsigned step,
                          DRIVER & driver );

        template< typename SETFUN, typename DRIVER >
        void restartReadData( std::istream & inp, SETFUN & setFun,
                              DRIVER & driver );

        void restartRemove( const std::string & basename,
                            const unsigned step );

    }
}

//==============================================================================
//------------------------------------------------------------------------------
template<typename DRIVER>
void gshell::fem::restartWrite( const std::string & basename,
                                const unsigned step,
                                DRIVER & driver )
{
    typedef typename DRIVER::Node Node;
    typedef std::function<typename Node::VecDof(const Node*)> GetFun;
    const std::string s = std::to_string( step );

    // write disp
    const std::string dispFileName = basename + "_D_" + s + ".dat";
    std::ofstream disp( dispFileName.c_str() );
    GetFun getDisp = std::bind( &Node::giveDisplacements, std::placeholders::_1 );
    gshell::fem::restartWriteData( disp, getDisp, driver  );
    disp.close();
        
    // write veloc
    const std::string velFileName = basename + "_V_" + s + ".dat";
    std::ofstream vel( velFileName.c_str() );
    GetFun getVel = std::bind( &Node::giveVelocities, std::placeholders::_1 );
    gshell::fem::restartWriteData( vel, getVel, driver );
    vel.close();

    // write accel
    const std::string accFileName = basename + "_A_" + s + ".dat";
    std::ofstream acc( accFileName.c_str() );
    GetFun getAcc = std::bind( &Node::giveAccelerations, std::placeholders::_1 );
    gshell::fem::restartWriteData( acc, getAcc, driver );
    acc.close();

    return;
}

//------------------------------------------------------------------------------
template< typename GETFUN, typename DRIVER >
void gshell::fem::restartWriteData( std::ostream & out, GETFUN & getFun,
                                    DRIVER & driver )
{
    typedef typename DRIVER::Node Node;
    typedef std::function<unsigned( const Node* )> GetId;
    GetId getId = std::bind( &Node::giveId, std::placeholders::_1 );
    fsi::utils::DataWriter<GetId,GETFUN> wf( getId, getFun, out );
    driver.accessMesh().iterateOverNodes( wf );
}

//------------------------------------------------------------------------------
template<typename DRIVER>
void gshell::fem::restartRead( const std::string & basename,
                                   const unsigned step,
                                   DRIVER & driver )
{
    typedef typename DRIVER::Node Node;
    typedef std::function<void(Node*,const typename Node::VecDof &)> SetFun;
    const std::string s = std::to_string( step );

    // write disp
    const std::string dispFileName = basename + "_D_" + s + ".dat";
    std::ifstream disp( dispFileName.c_str() );
    SetFun setDisp = std::bind( &Node::setDisplacements, std::placeholders::_1, 
                                std::placeholders::_2 );
    gshell::fem::restartReadData( disp, setDisp, driver );
    disp.close();
        
    // write veloc
    const std::string velFileName = basename + "_V_" + s + ".dat";
    std::ifstream vel( velFileName.c_str() );
    SetFun setVel = std::bind( &Node::setVelocities, std::placeholders::_1, 
                               std::placeholders::_2 );
    gshell::fem::restartReadData( vel, setVel, driver );
    vel.close();

    // write accel
    const std::string accFileName = basename + "_A_"+ s + ".dat";
    std::ifstream acc( accFileName.c_str() );
    SetFun setAcc = std::bind( &Node::setAccelerations, std::placeholders::_1,
                               std::placeholders::_2 );
    gshell::fem::restartReadData( acc, setAcc, driver );
    acc.close();

    return;
}

//------------------------------------------------------------------------------
template< typename SETFUN, typename DRIVER >
void gshell::fem::restartReadData( std::istream & inp, SETFUN & setFun,
                                   DRIVER & driver )
{
    typedef typename DRIVER::Node Node;
    typedef std::function<unsigned( const Node* )> GetId;
    GetId getId = std::bind( &Node::giveId, std::placeholders::_1 );
    fsi::utils::DataReader<GetId,SETFUN> sr( getId, setFun, inp );
    driver.accessMesh().iterateOverNodes( sr );
}

//------------------------------------------------------------------------------
/// Remove restart files of given time step
void gshell::fem::restartRemove( const std::string & basename,
                                 const unsigned step )
{
    const std::string s = std::to_string( step );

    const std::string dispFileName = basename + "_D_" + s + ".dat";
    remove( dispFileName.c_str() );
    const std::string velFileName = basename + "_V_" + s + ".dat";
    remove( velFileName.c_str() );
    const std::string accFileName = basename + "_A_" + s + ".dat";
    remove( accFileName.c_str() );
}


#endif
