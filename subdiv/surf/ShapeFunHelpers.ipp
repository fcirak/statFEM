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

//==============================================================================
template< >
double subdiv::surf::shapeMinDistanceToBoundary< corlib::TRIANGLE >(
    const eigenX::VectorSd< corlib::ShapeTraits< corlib::TRIANGLE >::dim > & xi
    )
{
    const unsigned numVertices = corlib::ShapeTraits< corlib::TRIANGLE >::numVertices;
    const std::array< double, numVertices > allDist =
        { { xi(0), xi(1), 1.-xi(0)-xi(1) } };
    const double minDistance = *std::min_element( allDist.begin(), allDist.end() );
    return minDistance;
}

template< >
double subdiv::surf::shapeMinDistanceToBoundary< corlib::QUADRILATERAL >(
    const eigenX::VectorSd< corlib::ShapeTraits< corlib::QUADRILATERAL >::dim > & xi
    )
{
    const unsigned numVertices = corlib::ShapeTraits< corlib::QUADRILATERAL >::numVertices;
    const std::array< double, numVertices > allDist =
        { { xi(0), xi(1), 1.-xi(0), 1.-xi(1) } };
    const double minDistance = *std::min_element( allDist.begin(), allDist.end() );
    return minDistance;
}


//==============================================================================
template< typename FACET >
void subdiv::surf::ShapeChildCoordinate< corlib::TRIANGLE, FACET >::operator()(
    Facet * facet,
    VecLDim & xi,
    MatLDimLDim & jac,
    unsigned & child
    )
{
    static_assert( myShape == FACET::myShape );

    // parametric co-ordinate short-cuts
    const double u = xi( 0 );
    const double v = xi( 1 );
    const double w = 1. - u - v;

    if ( w > .5 ) {
        child = 0;
        xi( 0 ) = 2.*u;
        xi( 1 ) = 2.*v;
        jac *= 2.;        
    }
    else if ( u > .5 ) {
        child = 1;
        xi( 0 ) = 2.*u - 1.;
        xi( 1 ) = 2.*v;
        jac *= 2.;
    }
    else if ( v > .5 ) {
        child = 2;
        xi( 0 ) = 2.*u;
        xi( 1 ) = 2.*v - 1.;
        jac *= 2.;
    }
    else {
        child = 3;
        xi( 0 ) = 1. - 2.*w;
        xi( 1 ) = 1. - 2.*u;
        MatLDimLDim trm;
        trm( 0, 0 ) =  2.;  trm( 0, 1 ) = 2.;
        trm( 1, 0 ) = -2.;  trm( 1, 1 ) = 0.;
        jac = trm * jac;
    }

    return;
}

template< typename FACET >
void subdiv::surf::ShapeChildCoordinate< corlib::QUADRILATERAL, FACET >::operator()(
    Facet * facet,
    VecLDim & xi,
    MatLDimLDim & jac,
    unsigned & child
    )
{
    static_assert( myShape == FACET::myShape );

    const double u = xi( 0 );
    const double v = xi( 1 );

    if ( ( u <= 0.5 ) and ( v <= 0.5 ) )  {
        child = 0;
    }
    else if ( ( u > 0.5 ) and ( v <= 0.5 ) ) {
        child = 1;
    }
    else if ( ( u > 0.5 ) and ( v > 0.5 ) ) {
        child = 2;
    }
    else if ( ( u <= 0.5 ) and ( v > 0.5) ) {
        child = 3;
    }

    // compute new co-ordinate and jacobian
    switch ( child ) {
        case 0 :   xi( 0 ) = 2.*u;        xi( 1 ) = 2.*v;        break;
        case 1 :   xi( 0 ) = 2.*u - 1.;   xi( 1 ) = 2.*v;        break;
        case 2 :   xi( 0 ) = 2.*u - 1.;   xi( 1 ) = 2.*v - 1.;   break;
        case 3 :   xi( 0 ) = 2.*u;        xi( 1 ) = 2.*v - 1.;   break;
        default :  FTL_VERIFY( false );
    }
    jac *= 2.;

    FTL_VERIFY_DESCRIPTIVE( corlib::fuzzyGreaterEqual(xi(0),0) and corlib::fuzzyLessEqual(xi(0),1), "xi(0) out of range [0,1]");
    FTL_VERIFY_DESCRIPTIVE( corlib::fuzzyGreaterEqual(xi(0),0) and corlib::fuzzyLessEqual(xi(0),1), "xi(1) out of range [0,1]");

    return;
}
