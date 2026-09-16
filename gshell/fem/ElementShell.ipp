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
//! @date   2012

#include <algorithm>
#include <corlib/NodeBasic.hpp>
#include <corlib/misc.hpp>

//------------------------------------------------------------------------------
// Tensor indices     \alpha\beta=  11,  12,  21,  22
// C indices                        00,  01,  10,  11
// Access   2*\alpha+\beta           0,   1    2    3
// 3-Voigt C-indices                 0    1    1    2
template< typename NODE, typename SFUN >
const unsigned gshell::fem::ElementShell<NODE,SFUN>::
voigtForward[localDim][localDim] = { 
    { ShapeFun::voigtForward[0][0], ShapeFun::voigtForward[0][1] },
    { ShapeFun::voigtForward[1][0], ShapeFun::voigtForward[1][1] }
};

//------------------------------------------------------------------------------
template< typename NODE, typename SFUN >
void gshell::fem::ElementShell<NODE,SFUN>::setLimitCoefficientsAtNodes()
{
    const unsigned numFunctions = this->numFunctions();
    FTL_VERIFY( numFunctions > 0 );
    
    for ( unsigned n = 0; n < numNodes; ++n ) {
        VecNF limitCoeff( numFunctions );
        shapeFun_.evaluate( n, limitCoeff );
        if ( nodes_[ n ]->getNumberLimitCoefficients() == 0 ) {
            for ( unsigned f = 0; f < numFunctions; ++f ) {
                nodes_[ n ]->addLimitCoefficient( supportNodes_[ f ], limitCoeff( f ) );
            }
        }
    }
}

//------------------------------------------------------------------------------
template< typename NODE, typename SFUN >
void gshell::fem::ElementShell<NODE,SFUN>::setUniqueTangentCoefficientsAtNodes()
{
    const unsigned numFunctions = this->numFunctions();
    FTL_VERIFY( numFunctions > 0 );
    
    for ( unsigned n = 0; n < numNodes; ++n ) {
        MatLDimNF tangCoeff;
        tangCoeff.resize( localDim, numFunctions );
        shapeFun_.evaluateGradient( n, tangCoeff );
        for ( unsigned d = 0; d < localDim; ++d ) {
            if ( nodes_[ n ]->getNumberTangentCoefficients( d ) == 0 ) {
                for ( unsigned f = 0; f < numFunctions; ++f ) {
                    nodes_[ n ]->addTangentCoefficient( supportNodes_[ f ], d, tangCoeff( d, f ) );
                }
            }
        }
    }
}

//------------------------------------------------------------------------------
template< typename NODE, typename SFUN >
void gshell::fem::ElementShell<NODE,SFUN>::setNodePtr( const unsigned index, 
                                                       Node * node )
{
    nodes_[ index ] = node;
    return;
}

//------------------------------------------------------------------------------
template< typename NODE, typename SFUN >
std::istream & gshell::fem::ElementShell<NODE,SFUN>::readSelf(
    std::istream & inp,
    const std::vector< Node * > & nodes
    )
{
    // read node IDs and set pointers
    for ( unsigned v = 0; v < numNodes; v ++ ) {
        unsigned vId;
        inp >> vId; 
        nodes_[v] = nodes[vId];
        if ( inp.peek() == ',' ) inp.ignore();
    }
    return inp;
}

//------------------------------------------------------------------------------
template< typename NODE, typename SFUN >
void gshell::fem::ElementShell<NODE,SFUN>::setSupportNodes(
    const NodeVector & supportNodes
    )
{
    supportNodes_.insert( supportNodes_.begin( ),
                          supportNodes.begin( ), supportNodes.end( ) );
}

//------------------------------------------------------------------------------
template< typename NODE, typename SFUN >
std::ostream & gshell::fem::ElementShell<NODE,SFUN>::writeNodeIndices(
    std::ostream & out
    ) const
{    
    for ( unsigned i = 0; i < numNodes; i ++ ) 
        nodes_[ i ] -> writeId( out );
    out << std::endl;
    return out;
}

//------------------------------------------------------------------------------
template< typename NODE, typename SFUN >
std::ostream & gshell::fem::ElementShell<NODE,SFUN>::writeSupportNodeIndices(
    std::ostream & out
    ) const
{
    for ( unsigned n=0; n<supportNodes_.size( ); ++n )
        out << supportNodes_[ n ] -> writeId( out );
    out << std::endl;
    return out;
}


//------------------------------------------------------------------------------
template< typename NODE, typename SFUN >
template< typename OUT >
OUT gshell::fem::ElementShell<NODE,SFUN>::giveNodeIndices( OUT oIter ) const
{
    for( unsigned v = 0; v < numNodes; v ++ ) 
        *oIter++ = nodes_[v] -> giveId();
    return oIter;
}

//------------------------------------------------------------------------------
template< typename NODE, typename SFUN >
template< typename OP > 
OP gshell::fem::ElementShell<NODE,SFUN>::iterateOverNodes( OP op ) const
{
    return std::for_each( nodes_.begin(), nodes_.end(), op );
}

//------------------------------------------------------------------------------
template< typename NODE, typename SFUN >
template< typename OP > 
OP gshell::fem::ElementShell<NODE,SFUN>::iterateOverNodes( OP op )
{
    return std::for_each( nodes_.begin(), nodes_.end(), op );
}

//------------------------------------------------------------------------------
template<typename NODE, typename SFUN>
void gshell::fem::ElementShell<NODE,SFUN>::nodalCoordinates( MatDimNN & X ) const
{
    for ( unsigned n=0; n<nodes_.size(); ++n )
        for ( unsigned d=0; d<dim; ++d )
            X( d, n ) = nodes_[ n ] -> giveCoordinates( )[ d ];
    return;
}

//------------------------------------------------------------------------------
template<typename NODE, typename SFUN>
void gshell::fem::ElementShell<NODE,SFUN>::supportNodeCoordinates( MatDimNF & X ) const
{
    X.resize( dim, supportNodes_.size( ) );
    for ( unsigned n=0; n<supportNodes_.size( ); ++n )
        for ( unsigned d=0; d<dim; ++d )
            X( d, n ) = supportNodes_[ n ] -> giveCoordinates( )[ d ];
    return;
}

//------------------------------------------------------------------------------
template< typename NODE, typename SFUN >
void gshell::fem::ElementShell<NODE,SFUN>::sfun( const VecLDim & xi, 
                                                 VecNF & phi ) const
{
    const unsigned numFunctions = this -> numFunctions( );
    FTL_VERIFY( shapeFun_.numFunctions( ) == numFunctions );
    phi.resize( numFunctions );
    shapeFun_.evaluate( xi, phi );
    return;
}

//------------------------------------------------------------------------------
template< typename NODE, typename SFUN >
void gshell::fem::ElementShell<NODE,SFUN>::sfunGrad( const VecLDim & xi, 
                                                     MatLDimNF & dPhiDXi ) const
{
    const unsigned numFunctions = this -> numFunctions( );
    FTL_VERIFY( shapeFun_.numFunctions() == numFunctions );
    BOOST_STATIC_ASSERT( ShapeFun::localDim == localDim );
    dPhiDXi.resize( localDim, numFunctions );
    shapeFun_.evaluateGradient( xi, dPhiDXi );
    return;
}

//------------------------------------------------------------------------------
template< typename NODE, typename SFUN >
void gshell::fem::ElementShell<NODE,SFUN>::sfunHess( const VecLDim & xi,
                                                     MatSDimNF & ddPhiDDXi ) const
{
    const unsigned numFunctions = this -> numFunctions( );
    FTL_VERIFY( shapeFun_.numFunctions() == numFunctions );
    BOOST_STATIC_ASSERT( ShapeFun::sDim == sDim );
    ddPhiDDXi.resize( sDim, numFunctions );
    shapeFun_.evaluateHessian( xi, ddPhiDDXi );
    return;
}

//------------------------------------------------------------------------------
template< typename NODE, typename SFUN >
void gshell::fem::ElementShell<NODE,SFUN>::sfunGradHess( const VecLDim & xi,
                                                         VecNF & phi,
                                                         MatLDimNF & dPhiDXi,
                                                         MatSDimNF & ddPhiDDXi ) const
{
    const unsigned numFunctions = this -> numFunctions( );
    phi.resize( numFunctions );
    dPhiDXi.resize( localDim, numFunctions );
    ddPhiDDXi.resize( sDim, numFunctions );
    shapeFun_.evaluateGradHess( xi, phi, dPhiDXi, ddPhiDDXi );
    return;
}

//------------------------------------------------------------------------------
template< typename NODE, typename SFUN >
void gshell::fem::ElementShell<NODE,SFUN>::sfun( const unsigned nIndex,
                                                 VecNF & coeff ) const
{
    const unsigned numFunctions = this -> numFunctions( );
    FTL_VERIFY( numFunctions == supportNodes_.size() );
    coeff.resize( numFunctions );
    shapeFun_.evaluate( nIndex, coeff );
}

//------------------------------------------------------------------------------
template< typename NODE, typename SFUN >
void gshell::fem::ElementShell<NODE,SFUN>::sfunGrad( const unsigned nIndex,
                                                     MatLDimNF & coeff ) const
{
    const unsigned numFunctions = this -> numFunctions( );
    FTL_VERIFY( numFunctions == supportNodes_.size() );
    coeff.resize( localDim, numFunctions );
    shapeFun_.evaluateGradient( nIndex, coeff );
}

//------------------------------------------------------------------------------
template< typename NODE, typename SFUN >
void gshell::fem::ElementShell<NODE,SFUN>::geometry( const VecLDim & xi, 
                                                     VecDim  & x ) const
{
    // check xi
    FTL_VERIFY( corlib::ShapeTraits< myShape >::isInside( xi ) );
    // nodal coordinates
    MatDimNF X;  this -> supportNodeCoordinates( X ); 
    // shape function evaluation
    VecNF phi;   this -> sfun( xi, phi );
    // compute product
    x = X * phi;
    return;
}

//------------------------------------------------------------------------------
template< typename NODE, typename SFUN >
void gshell::fem::ElementShell<NODE,SFUN>::geometry( const unsigned nIndex,
                                                     VecDim  & x ) const
{
    // nodal coordinates
    MatDimNF X;  this -> supportNodeCoordinates( X ); 
    // limit surface coefficient, like shape function evaluation
    VecNF coeff;  this -> sfun( nIndex, coeff );
    // compute product
    x = X * coeff;
    return;
}

//------------------------------------------------------------------------------
template< typename NODE, typename SFUN >
double gshell::fem::ElementShell<NODE,SFUN>::jacobian( const VecLDim & xi ) const
{
    // nodal coordinates
    MatDimNF xNF;      this -> supportNodeCoordinates( xNF );
    // shape function derivatives
    MatLDimNF dPhiDXi; this -> sfunGrad( xi, dPhiDXi );
    // covariant basis
    typedef eigenX::MatrixSd<dim,localDim> MatDimLDim;
    const MatDimLDim coG = xNF * dPhiDXi.transpose( );
    // return its determinant
    return corlib::metric( coG );
}

//------------------------------------------------------------------------------
template< typename NODE, typename SFUN >
void gshell::fem::ElementShell<NODE,SFUN>::areaIntegrand( const VecLDim & xi,
                                                          const double & weight ) 
{
    surfaceArea_ += this -> jacobian( xi ) * weight;
    return;
}

//------------------------------------------------------------------------------
template< typename NODE, typename SFUN >
void gshell::fem::ElementShell<NODE,SFUN>::clearNodes_( ) 
{
    for ( unsigned i = 0; i < numNodes; i ++ ) nodes_[ i ] = NULL;
    return;
}

