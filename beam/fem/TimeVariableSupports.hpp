// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//
//                 Computational Structural Mechanics Lab
//                         University of Cambridge
//
// This software is copyrighted and all rights are retained by the University of
// Cambridge. All use, disclosure, and/or reproduction of any part not expressly
// authorised by the University of Cambridge is prohibited. (C) 2010
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

//! @file TimeVariableSupports.hpp

//! @todo   Add DISP, DISPX and ROT constraints. Inspect Supports.ipp

#ifndef beam_fem_timevariablesupports_h
#define beam_fem_timevariablesupports_h

#include <iostream>
#include <array>

#include <Eigen/Core>
#include <corlib/eigenX.hpp>

#include <beam/fem/LinkBasic.hpp>
#include <beam/fem/Supports.hpp>


//==============================================================================
namespace beam {
    namespace fem {

        namespace eigenX = corlib::eigenX;

        template< typename LINK >
        class TimeVariableSupport;

    }
}

//==============================================================================
/** \brief A Class to allow a time-varying factor to be applied to supports.
 *   The user specifies the target link to which the factor is applied. This
 *   functor is iterated over beam::fem::Mesh's Links_
 */
template< typename LINK >
class beam::fem::TimeVariableSupport
{
public:
    typedef LINK                               LinkType;

private:
    typedef typename LinkType::Node            Node_;

public:
    /// Constructor
    TimeVariableSupport( const beam::fem::supportType type,
                         const LinkType* target,
                         const double factor )
        : type_(   type ), 
          target_( target ), 
          factor_( factor )
    { }

    /// Overloaded function call operator
    void operator()( LinkType* lp )
    {
        if ( target_ == lp ) {

            switch ( type_ ) {
            case beam::fem::ROTZ :
            {
                // Extract the linked nodes
                std::array< const Node_ *, 4 > nodes;
                std::array< unsigned, 4 >      dofs;
                std::array< double, 4 >        coeffs;
                for ( unsigned i = 0; i < 4; ++i )
                    lp->getNodeDofCoeff( i, nodes[ i ], dofs[ i ], coeffs[ i ] );
                std::array< const Node_ *, 3 > nodesTriple;
                nodesTriple[ 0 ] = nodes[ 0 ];
                nodesTriple[ 1 ] = static_cast< const Node_ * >( NULL );
                nodesTriple[ 2 ] = nodes[ 1 ];

                // new support axis
                eigenX::VectorSd<3> axis; axis.setZero();
                axis[ 2 ] = 1.;

                // create the new link
                beam::fem::Support supp( type_, factor_, axis );
                beam::fem::SupportsAtNode< LinkType, 3 > suppAtNode;
                std::vector< LinkType * > links;
                suppAtNode.set( links, supp, nodesTriple );

                // get nodeDofCoeff tuple corresponding to the new link
                for ( unsigned i = 0; i < 4; ++i )
                    (links[ 0 ])->getNodeDofCoeff( i, nodes[ i ], dofs[ i ], coeffs[ i ] );
            
                // modify lp to match the new link
                for ( unsigned i = 0; i < 4; ++i )
                    lp->setNodeDofCoeff( i, nodes[ i ], dofs[ i ], coeffs[ i ]  );
                lp->setValue( (links[ 0 ])->getValue() );
            
                // delete the new link
                delete links[ 0 ];
            
                break;
            } // case

            case beam::fem::DISPY :
            {
                // Extract the linked nodes
                std::array< const Node_ *, 6 > nodes;
                std::array< unsigned, 6 >      dofs;
                std::array< double, 6 >        coeffs;
                for ( unsigned i = 0; i < 6; ++i )
                    lp->getNodeDofCoeff( i, nodes[ i ], dofs[ i ], coeffs[ i ] );
                std::array< const Node_ *, 3 > nodesTriple;
                nodesTriple[ 0 ] = nodes[ 0 ];
                nodesTriple[ 1 ] = nodes[ 1 ];
                nodesTriple[ 2 ] = nodes[ 2 ];

                // new support axis
                eigenX::VectorSd<3> axis; axis.setZero();
                axis[ 1 ] = 1.;

                // create the new link
                beam::fem::Support supp( type_, factor_, axis );
                beam::fem::SupportsAtNode< LinkType, 3 > suppAtNode;
                std::vector< LinkType * > links;
                suppAtNode.set( links, supp, nodesTriple );

                // get nodeDofCoeff tuple corresponding to the new link
                for ( unsigned i = 0; i < 6; ++i )
                    (links[ 0 ])->getNodeDofCoeff( i, nodes[ i ], dofs[ i ], coeffs[ i ] );
            
                // modify lp to match the new link
                for ( unsigned i = 0; i < 6; ++i )
                    lp->setNodeDofCoeff( i, nodes[ i ], dofs[ i ], coeffs[ i ]  );
                lp->setValue( (links[ 0 ])->getValue() );

                // delete the new link
                delete links[ 0 ];
            
                break;
            } // case

            default : {
                FTL_VERIFY_DESCRIPTIVE( false,
                                        "Cannot apply time variation to chosen support type" );
            } // default

            } // switch

        } // if
        
        return;
    }

private:
    const double                       factor_;   ///< factor to apply
    const LinkType*                    target_;   ///< apply to this link
    const enum beam::fem::supportType  type_;     ///< type of support
};

#endif
