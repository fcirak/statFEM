// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//
//                 Computational Structural Mechanics Lab
//                         University of Cambridge
//
// This software is copyrighted and all rights are retained by the University of
// Cambridge. All use, disclosure, and/or reproduction of any part not expressly
// authorised by the University of Cambridge is prohibited. (C) 2010
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

//! @file CouplerBeam.ipp

//------------------------------------------------------------------------------
/** Predict nodal increments with a first-order Taylor series:
 *  \f[
 *      (x^{(0)}_{n+1} - x_n) = \Delta t * v_n
 *  \f]
 *  Steps done by this function:
 *  -  set nodal increments to this value in nodes of the beam mesh
 *  -  copy the new nodal increments to the given array
 */
template< typename DRIVER >
void beam::fem::CouplerBeam<DRIVER>::
predictIncrementsFE( Mesh & mesh,
                     const double stepSize,
                     std::vector<VecDim> & nodeIncrements )
{
    // set nodal increments in the beam driver
    std::function<void(Node*)> predictIncr =
        std::bind( &Node::setIncrement, std::placeholders::_1,
                     std::bind( beam::fem::Multiply<VecDim,double,VecDim>(),
                                  stepSize, std::bind( &Node::giveVelocities, std::placeholders::_1 ) ) );
    mesh.iterateOverNodes( predictIncr );

    // copy to interface
    nodeIncrements.resize( mesh.numNodes() );
    std::transform( mesh.nodesBegin(),
                    mesh.nodesEnd(),
                    nodeIncrements.begin(),
                    std::bind( &Node::giveIncrement, std::placeholders::_1 ) );
}

//------------------------------------------------------------------------------
//! Pass the increments given in a vector to the nodes
template< typename DRIVER >
void beam::fem::CouplerBeam<DRIVER>::
predictIncrementsLI( Mesh & mesh,
                     const std::vector<VecDim> & nodeIncrements )
{
    if ( not nodeIncrements.empty() ) {
        corlib::DistributeQuantity<Node,VecDim>
            distributeI( std::bind( &Node::setIncrement, std::placeholders::_1, std::placeholders::_2 ) );
        typename std::vector<VecDim>::const_iterator increIter = nodeIncrements.begin();
        mesh.iterateOverNodes( std::bind( distributeI, std::placeholders::_1,
                                            std::ref( increIter ) ) );
    }
}

//------------------------------------------------------------------------------
// Predict nodal displacement increments to be zero
template< typename DRIVER >
void beam::fem::CouplerBeam<DRIVER>::
predictIncrementsZero( Mesh & mesh,
                       std::vector<VecDim> & nodeIncrements )
{
    mesh.iterateOverNodes( std::bind( &Node::clearIncrement, std::placeholders::_1 ) );
    VecDim zero; zero.clear();
    nodeIncrements.resize( mesh.numNodes(), zero );
}

//------------------------------------------------------------------------------
// Copy current coordinates of shell to surface
template< typename DRIVER >
void beam::fem::CouplerBeam<DRIVER>::
getSurfaceCoordinates( Mesh & mesh,
                       std::vector<VecDim> & nodeCoordinates )
{
    nodeCoordinates.resize( mesh.numNodes() );
    std::transform( mesh.nodesBegin(),
                    mesh.nodesEnd(),
                    nodeCoordinates.begin(),
                    std::bind( std::plus<VecDim>(),
                                 std::bind( &Node::giveCoordinates, std::placeholders::_1 ),
                                 std::bind( &Node::giveDisplacements, std::placeholders::_1 ) ) );
}

//------------------------------------------------------------------------------
// Copy element connectivity from shell to surface
template< typename DRIVER >
void beam::fem::CouplerBeam<DRIVER>::
getSurfaceConnectivity( Mesh & mesh,
                        std::vector<ArrayNN> & elementConnectivities )
{
    elementConnectivities.resize( mesh.numElements() );
    typename Driver::Mesh::ElementConstIterator eIter =
        mesh.elementsBegin();
    typename std::vector<ArrayNN>::iterator ecIter = elementConnectivities.begin();
    for ( ; eIter != mesh.elementsEnd(); ++eIter, ++ecIter ) {
        typename Element::NodeConstIterator nIter = (*eIter)->nodesBegin();
        typename ArrayNN::iterator niIter = ecIter->begin();
        for ( ; nIter != (*eIter)->nodesEnd(); ++nIter, ++niIter ) {
            *niIter = (*nIter)->giveId();
        }
    }
}

//------------------------------------------------------------------------------
// Current surface velocity as a finite difference
template< typename DRIVER >
void beam::fem::CouplerBeam<DRIVER>::
getSurfaceVelocitiesFD( Mesh & mesh,
                        const double stepSize,
                        std::vector<VecDim> & nodeVelocities )
{
    nodeVelocities.resize( mesh.numNodes() );
    typename std::vector<VecDim>::iterator velocIter = nodeVelocities.begin();

    // extract velocity v_{n+1} = 1/\Delta{t} * \Delta{d}_{n+1}^{<k>}
    corlib::CollectQuantity<Node,VecDim>
        collectV( std::bind( beam::fem::Multiply<VecDim,double,VecDim>(),
                               1. / stepSize,
                               std::bind( &Node::giveIncrement, std::placeholders::_1 ) ) );
    
    mesh.iterateOverNodes( std::bind( collectV, std::placeholders::_1, std::ref( velocIter ) ) );
}

//------------------------------------------------------------------------------
// Copy current velocites
template< typename DRIVER >
void beam::fem::CouplerBeam<DRIVER>::
getSurfaceVelocities( Driver & driver, 
                      const double stepSize,
                      std::vector<VecDim> & nodeVelocities )
{
    // store velocities v_{n}
    corlib::CollectQuantity<Node,VecDim>
        collectV( std::bind( &Node::giveVelocities, std::placeholders::_1 ) );
    std::vector<VecDim> veloc( driver.accessMesh().numNodes() );
    typename std::vector<VecDim>::iterator velocIter = veloc.begin();
    driver.accessMesh().iterateOverNodes( std::bind( collectV, std::placeholders::_1,
                                                       std::ref( velocIter ) ) );

    // store accelerations a_{n}
    corlib::CollectQuantity<Node,VecDim>
        collectA( std::bind( &Node::giveAccelerations, std::placeholders::_1 ) );
    std::vector<VecDim> accel( driver.accessMesh().numNodes() );
    typename std::vector<VecDim>::iterator accelIter = accel.begin();
    driver.accessMesh().iterateOverNodes( std::bind( collectA, std::placeholders::_1,
                                                        std::ref( accelIter ) ) );

    // update ... this will overwrite v_n := v_{n+1}, a_n := a_{n+1}
    driver.dynamicUpdate( stepSize );

    // extract velocity v_{n+1}
    nodeVelocities.resize( driver.accessMesh().numNodes() );
    std::transform( driver.accessMesh().nodesBegin(),
                    driver.accessMesh().nodesEnd(),
                    nodeVelocities.begin(),
                    std::bind( &Node::giveVelocities, std::placeholders::_1 ) );

    // reset previously saved velocities v_{n}
    corlib::DistributeQuantity<Node,VecDim>
        distributeV( std::bind( &Node::setVelocities, std::placeholders::_1, std::placeholders::_2 ) );
    velocIter = veloc.begin();
    driver.accessMesh().iterateOverNodes( std::bind( distributeV, std::placeholders::_1,
                                                       std::ref( velocIter ) ) );

    // reset previously saved accelerations a_{n}
    corlib::DistributeQuantity<Node,VecDim>
        distributeA( std::bind( &Node::setAccelerations, std::placeholders::_1, std::placeholders::_2 ) );
    accelIter = accel.begin();
    driver.accessMesh().iterateOverNodes( std::bind( distributeA, std::placeholders::_1,
                                                        std::ref( accelIter ) ) );
}

//------------------------------------------------------------------------------
// Copy current coordinates of shell to surface
template< typename DRIVER >
void beam::fem::CouplerBeam<DRIVER>::
getSurfaceNormals( Mesh & mesh,
                   std::vector<VecDim> & elementCentres,
                   std::vector<VecDim> & elementNormals,
                   const bool onePerElement )
{
    // one evalution point per element
    if ( onePerElement ) {
        const typename Element::VecLDim centre =
            corlib::ShapeTraits<Element::myShape>::centroid();
    
        elementCentres.resize( mesh.numElements() );
        std::transform( mesh.elementsBegin(),
                        mesh.elementsEnd(),
                        elementCentres.begin(),
                        std::bind( std::plus<VecDim>(),
                                     std::bind( &Element::giveCoordinate, std::placeholders::_1, centre ),
                                     std::bind( &Element::giveDisplacement, std::placeholders::_1, centre ) ) );

        elementNormals.resize( mesh.numElements() );
        std::transform( mesh.elementsBegin(),
                        mesh.elementsEnd(),
                        elementNormals.begin(),
                        std::bind( &Element::giveScaledNormal, std::placeholders::_1, centre ) );
    }
    // multiple evalution points per element using the quadrature rule
    else {
        typename Driver::Quadrature quadrature;
        const unsigned numCentres = mesh.numElements()
            * Driver::Quadrature::numPoints;
        elementCentres.resize( numCentres );
        elementNormals.resize( numCentres );
        typename std::vector<VecDim>::iterator ecIter = elementCentres.begin();
        typename std::vector<VecDim>::iterator enIter = elementNormals.begin();
        
        typename Driver::Mesh::ElementConstIterator eIter
            = mesh.elementsBegin();
        for ( ; eIter != mesh.elementsEnd() ; ++eIter ) {
            typename Driver::Quadrature::QuadIter qIter = quadrature.begin();
            for ( ; qIter != quadrature.end(); ++qIter ) {
                *ecIter = (*eIter)->giveCoordinate( qIter->second )
                    + (*eIter)->giveDisplacement( qIter->second );
                *enIter = (*eIter)->giveScaledNormal( qIter->second );
                ++ecIter;
                ++enIter;
            }
        }
    }
}

//------------------------------------------------------------------------------
// Copy force on surface nodes to shell nodes
template< typename DRIVER >
void beam::fem::CouplerBeam<DRIVER>::
setSurfaceForceExt( Driver & driver,
                    const std::vector<VecDim> & elementTractions,
                    const bool onePerElement )
{
    typedef typename std::vector<VecDim> VecVecDim;

    // one evalution point per element
    if ( onePerElement ) {
        FTL_VERIFY( elementTractions.size() == driver.accessMesh().numElements() );

        typedef fsi::utils::BodyForceDiscreteConstant<Element,VecVecDim> PCBodyForce;
        typedef typename PCBodyForce::LoadFun                            LoadFun;
        PCBodyForce pcbf( std::bind( &Element::template bodyForce<LoadFun>,
                                       std::placeholders::_1, std::placeholders::_2, 
                                       std::placeholders::_3, std::placeholders::_4, 1. ),
                          elementTractions );

        driver.setExternalTractionFun( std::bind( pcbf, std::placeholders::_1, std::placeholders::_2, 
                                                  std::placeholders::_3 ) );
    }
    // multiple evalution points per element using the quadrature rule
    else {
        FTL_VERIFY( elementTractions.size() == ( driver.accessMesh().numElements()
                                                 * Driver::Quadrature::numPoints ) );

        typedef fsi::utils::BodyForceDiscrete<Element,VecVecDim>    PCBodyForce;
        typedef typename PCBodyForce::LoadFun                       LoadFun;
        PCBodyForce pcbf( std::bind( &Element::template bodyForce<LoadFun>,
                                     std::placeholders::_1, std::placeholders::_2, 
                                     std::placeholders::_3, std::placeholders::_4, 1. ),
                          elementTractions,
                          Driver::Quadrature::numPoints );

        driver.setExternalTractionFun( std::bind( pcbf, std::placeholders::_1, std::placeholders::_2, 
                                                  std::placeholders::_3 ) );
    }

}

//------------------------------------------------------------------------------
/** Return the vector-norm of the difference between previous and current increments
 *  \f[
 *      \| (x_{n+1}^{(k+1)} - x_{n}) - (x_{n+1}^{(k)} - x_{n}) \|_2
 *  \f]
 */
template< typename DRIVER >
double beam::fem::CouplerBeam<DRIVER>::
incrementResidual( Mesh & mesh,
                   const std::vector<VecDim>  & nodeIncrementsPrev )
{
    // residual displacements (increment of displacement increments)
    const double normedResIncrSq =
        std::inner_product( mesh.nodesBegin(),
                            mesh.nodesEnd(),
                            nodeIncrementsPrev.begin(), 0.,
                            std::plus<double>(),
                            std::bind( &beam::fem::norm2Squared<VecDim>,
                                         std::bind( std::minus<VecDim>(),
                                                      std::bind( &Node::giveIncrement, 
                                                                 std::placeholders::_1 ), std::placeholders::_2 ) ) );
                        
    return std::sqrt( normedResIncrSq );
}


//------------------------------------------------------------------------------
// Update scaled increments by #relax factor
template< typename DRIVER >
double beam::fem::CouplerBeam<DRIVER>::
updateIncrements( Mesh & mesh,
                  std::vector<VecDim>  & nodeIncrementsPrev,
                  const double relax )
{
    if ( nodeIncrementsPrev.empty() ) {
        VecDim zero; zero.clear();
        nodeIncrementsPrev.resize( mesh.numNodes(), zero );
    }
    else {
        FTL_VERIFY( nodeIncrementsPrev.size() == mesh.numNodes() );
    }

    // incr = relax*incr + (1-relax)*incrPrev
    corlib::DistributeQuantity<Node,VecDim>
        scaleIncr( std::bind( &Node::setIncrement, std::placeholders::_1,
                                std::bind( std::plus<VecDim>(),
                                             std::bind( beam::fem::Multiply<VecDim,double,VecDim>(),
                                                          relax,
                                                          std::bind( &Node::giveIncrement, std::placeholders::_1 ) ),
                                             std::bind( beam::fem::Multiply<VecDim,double,VecDim>(),
                                                          1.-relax,
                                                          _2 ) ) ) );
    typename std::vector<VecDim>::const_iterator nodeIncrIter = nodeIncrementsPrev.begin();
    mesh.iterateOverNodes( std::bind( scaleIncr, std::placeholders::_1, std::ref( nodeIncrIter ) ) );
    
    // residual displacements (increment of displacement increments)
    const double normedResIncrSq =
        std::inner_product( mesh.nodesBegin(),
                            mesh.nodesEnd(),
                            nodeIncrementsPrev.begin(), 0.,
                            std::plus<double>(),
                            std::bind( &beam::fem::norm2Squared<VecDim>,
                                         std::bind( std::minus<VecDim>(),
                                                      std::bind( &Node::giveIncrement, std::placeholders::_1 ), std::placeholders::_2 ) ) );
                        
    // store increment
    std::transform( mesh.nodesBegin(),
                    mesh.nodesEnd(),
                    nodeIncrementsPrev.begin(),
                    std::bind( &Node::giveIncrement, std::placeholders::_1 ) );

    // clear the multiplier increments
    typedef typename Driver::Link Link;
    mesh.iterateOverLinks( std::bind( &Link::clearIncrement, std::placeholders::_1 ) );

    // return norm
    return std::sqrt( normedResIncrSq );
}
