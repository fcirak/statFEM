// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//
//                 Computational Structural Mechanics Lab
//                         University of Cambridge
//
// This software is copyrighted and all rights are retained by the University of
// Cambridge. All use, disclosure, and/or reproduction of any part not expressly
// authorised by the University of Cambridge is prohibited. (C) 2010
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

//! @file NodeMixed.hpp

#ifndef nodemixed_h
#define nodemixed_h

#include <solid/fem/NodeStatic.hpp>

//------------------------------------------------------------------------------
namespace solid {
    namespace fem {
        template< unsigned DIM, unsigned DOF > class NodeMixed;
        
        namespace ublas = boost::numeric::ublas;
    }
}

//------------------------------------------------------------------------------
//! extends the SolidNode by a pressure value
template<unsigned DIM, unsigned DOF>
class solid::fem::NodeMixed : public solid::fem::NodeStatic<DIM,DOF>
{
public:
    typedef ublas::bounded_vector< double, 1 >       Vec1;

private:
    typedef solid::fem::NodeStatic<DIM,DOF>          NodeStatic_;
    typedef typename NodeStatic_::VecDof             VecDof_;
    typedef typename NodeStatic_::ZeroVec_           ZeroVec_;

    using NodeStatic_::disp_;

public:
    NodeMixed( const unsigned & id ) : NodeStatic_( id ),
                                       pressure_( 0. ),
                                       pressFlag_( false ) { }


    double givePressure(  ) const { return pressure_; }

    void setPressureFlag( ) { pressFlag_ = true; return; }

    void copyDofArrayP( std::vector<unsigned> & dofIndices ) const
    { 
        dofIndices[0] = pressId_;
        return;
    }

    bool givePressureFlag() const { return pressFlag_;}

    void storePressure( const Vec1 & p ) { pressure_ = p( 0 ); return; }

    void setPressureId( unsigned & pressCounter ) 
    { 
        if ( pressFlag_ ) pressId_ = pressCounter++; 
        return;
    }


private:
    double      pressure_;
    unsigned    pressId_;
    bool        pressFlag_;
};

#endif
