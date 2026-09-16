// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//
//                 Computational Structural Mechanics Lab
//                         University of Cambridge
//
// This software is copyrighted and all rights are retained by the University of
// Cambridge. All use, disclosure, and/or reproduction of any part not expressly
// authorised by the University of Cambridge is prohibited. (C) 2010
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

//! @file   Supports.ipp

#include <iostream>
#include "Supports.hpp"


//==============================================================================

//------------------------------------------------------------------------------
//! Constructor with given file stream inp
beam::fem::SupportsFromFile::SupportsFromFile( std::istream & inp )
{
    // get number of supports
    unsigned numSupports = 0;
    inp >> numSupports;
    // go through constraints
    for ( unsigned c = 0; c < numSupports; ++c ) {
        // get node ID
        unsigned nodeNum;
        inp >> nodeNum;
        // store support in container
        supportContainer_.insert( std::make_pair( nodeNum, 
                                                  beam::fem::Support( inp ) ) ); 
        inp.ignore( std::numeric_limits<std::streamsize>::max(), '\n' );
    }
}

//------------------------------------------------------------------------------
//! Print supports
std::ostream& beam::fem::SupportsFromFile::print( std::ostream& os ) const
{
    for ( SupportConstIter el = this -> begin(); 
          el != this -> end(); ++ el ) {
        os << "nodeID=" << el -> first << " : ";
        (el -> second).print( os );
    }
    return os;
}


//==============================================================================

//------------------------------------------------------------------------------
//! Main method to apply supports for cubic splines
template< typename LBASIC >
void beam::fem::SupportsAtNode< LBASIC, 3 >::set( std::vector< LinkBasic_ * >& links,
                                                  const beam::fem::Support supp,
                                                  const std::array< const Node *, 3 >& nodes )
{
    typedef eigenX::VectorSd<3> Vec3;

    // store constraint value directly in a Vec1
    const typename LinkBasic_::VecDof value =
        Eigen::VectorXd::Constant( LinkBasic_::dof, supp.giveValue() );
        
    switch ( supp.giveType() ) {
        // rotation type constraint
    case beam::fem::ROT :
    case beam::fem::ROTZ : {
        FTL_VERIFY( dim == 2 ); // can't do 3D yet
        // construct secant
        const VecDim coord0 = nodes[ 0 ] -> giveCoordinates();
        const VecDim coord2 = nodes[ 2 ] -> giveCoordinates();
        const VecDim secant = 0.5 * (coord2 - coord0);
        // work in 3D
        Vec3 zAxis, tang, norm;
        // out-of plane axis
        zAxis.setZero(); 
        zAxis[2] = 1.;
        // tangent vector
        tang.setZero();
        tang.head( dim ) = secant;
        // normal vector
        norm = corlib::cross_prod( zAxis, tang );
        // rotate normal vector
        eigenX::MatrixSd< 3, 3 > rot = Eigen::MatrixXd::Identity( 3, 3 );
        rot( 0, 0 ) = std::cos( value[0] );  rot( 0, 1 ) = -std::sin( value[0] );
        rot( 1, 0 ) = std::sin( value[0] );  rot( 1, 1 ) =  std::cos( value[0] );
        norm = rot * norm;
        norm /= norm.norm( );
        // Back to d-D
        VecDim normal = norm.head( dim ); 
        //normal[ 0 ] =  1./2.*coord0[ 1 ] - 1./2.*coord2[ 1 ];
        //normal[ 1 ] = -1./2.*coord0[ 0 ] + 1./2.*coord2[ 0 ];

        LinkBasic_ * linkBasic = new LinkBasic_();
        for ( unsigned d = 0; d < dim; d ++ ) {
            linkBasic->addNodeDofCoeff( nodes[ 0 ], d, -1./2. * normal[ d ] );
            linkBasic->addNodeDofCoeff( nodes[ 2 ], d,  1./2. * normal[ d ] );
        }

        //linkBasic->setValue( ublas::zero_vector< double >( LinkBasic_::dof ) );
        const double rhs = - norm.dot( tang );
        linkBasic->setValue( Eigen::VectorXd::Constant( LinkBasic_::dof, rhs ) );
        links.push_back( linkBasic );

        break;
    }
        // displacement type constraint
    case beam::fem::DISP : 
    case beam::fem::DISPX :
    case beam::fem::DISPY :
    {
        LinkBasic_ * linkBasic = new LinkBasic_();
        const Vec3 axis = supp.giveAxis();
        for ( unsigned d = 0; d < dim; d ++ ) {
            const double fac = axis[d];
            linkBasic -> addNodeDofCoeff( nodes[ 0 ], d, fac * 1./6. );
            linkBasic -> addNodeDofCoeff( nodes[ 1 ], d, fac * 2./3. );
            linkBasic -> addNodeDofCoeff( nodes[ 2 ], d, fac * 1./6. );
        }
        
        linkBasic -> setValue( value );
        links.push_back( linkBasic );
        break;
    }
    default : {
        FTL_VERIFY_DESCRIPTIVE( false, "Cannot apply chosen support type" );
    }
    }

    return;
} 
