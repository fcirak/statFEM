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

//! @author Burkhard Bornemann
//! @date   2012

//------------------------------------------------------------------------------
// Copy current coordinates of shell to surface
template<typename DRIVER>
void gshell::fem::Coupler<DRIVER>::
predictIncrementsFE( Mesh & mesh,
                     const double stepSize, 
                     std::vector<VecDim> & nodeIncrements )
{
    // predict increments \Delta{d}_{n}^0 := \Delta{t}_n * v_n
    std::function<void(Node*)> predictIncr =
        std::bind( &Node::setIncrement, std::placeholders::_1,
                     std::bind( gshell::fem::Multiply<VecDim,double,VecDim>(),
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
// Predict nodal displacement increments using last increment
template<typename DRIVER>
void gshell::fem::Coupler<DRIVER>::
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
template<typename DRIVER>
void gshell::fem::Coupler<DRIVER>::
predictIncrementsZero( Mesh & mesh,
                       std::vector<VecDim> & nodeIncrements )
{
    mesh.iterateOverNodes( std::bind( &Node::clearIncrement, std::placeholders::_1 ) );
    VecDim zero; zero.clear();
    nodeIncrements.resize( mesh.numNodes(), zero );
}

//------------------------------------------------------------------------------
// Copy current coordinates of shell to surface
template<typename DRIVER>
void gshell::fem::Coupler<DRIVER>::
getSurfaceCoordinates( const Mesh & mesh,
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
void gshell::fem::Coupler<DRIVER>::
getSurfaceConnectivity( const Mesh & mesh,
                        std::vector<ArrayNN>& elementConnectivities )
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
template<typename DRIVER>
void gshell::fem::Coupler<DRIVER>::
getSurfaceVelocitiesFD( Mesh & mesh,
                        const double stepSize,
                        std::vector<VecDim> & nodeVelocities )
{
    nodeVelocities.resize( mesh.numNodes() );
    typename std::vector<VecDim>::iterator velocIter = nodeVelocities.begin();
    
    // extract velocity v_{n+1} = 1/\Delta{t} * \Delta{d}_{n}^{<k>}
    corlib::CollectQuantity<Node,VecDim>
        collectV( std::bind( gshell::fem::Multiply<VecDim,double,VecDim>(),
                               1. / stepSize,
                               std::bind( &Node::giveIncrement, std::placeholders::_1 ) ) );
    
    mesh.iterateOverNodes( std::bind( collectV, std::placeholders::_1, std::ref( velocIter ) ) );
}

//------------------------------------------------------------------------------
// Copy current coordinates of shell to surface
template<typename DRIVER>
void gshell::fem::Coupler<DRIVER>::
getSurfaceVelocities( Driver & driver,
                      const double stepSize, 
                      std::vector<VecDim> & nodeVelocities )
{
    typedef typename std::vector<VecDim>::iterator     VecVecDimIter;

    // store velocities v_{n}
    corlib::CollectQuantity<Node,VecDim>
        collectV( std::bind( &Node::giveVelocities, std::placeholders::_1 ) );
    std::vector<VecDim> veloc( driver.accessMesh().numNodes() );
    VecVecDimIter velocIter = veloc.begin();
    driver.accessMesh().iterateOverNodes( std::bind( collectV, std::placeholders::_1,
                                                        std::ref( velocIter ) ) );

    // store accelerations a_{n}
    corlib::CollectQuantity<Node,VecDim>
        collectA( std::bind( &Node::giveAccelerations, std::placeholders::_1 ) );
    std::vector<VecDim> accel( driver.accessMesh().numNodes() );
    VecVecDimIter accelIter = accel.begin();
    driver.accessMesh().iterateOverNodes( std::bind( collectA, std::placeholders::_1,
                                                        std::ref( accelIter ) ) );
    
    // update ... this will overwrite d_n := d_{n+1}, \Delta d_n := 0,
    //                                v_n := v_{n+1}, a_n := a_{n+1}
    driver.dynamicUpdate( stepSize );

    // extract velocity v_{n} = v_{n+1}
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
template<typename DRIVER>
void gshell::fem::Coupler<DRIVER>::
getSurfaceNormals( const Mesh & mesh,
                   std::vector<VecDim> & elementCentres,
                   std::vector<VecDim> & elementNormals )
{
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

//------------------------------------------------------------------------------
// Copy force on surface nodes to shell nodes
template<typename DRIVER>
void gshell::fem::Coupler<DRIVER>::
setSurfaceForceExt( Driver & driver,
                    const std::vector<VecDim> & elementTractions )
{
    FTL_VERIFY( elementTractions.size() == driver.accessMesh().numElements() );

    typedef std::vector<VecDim>                           VecVecDim;
    typedef fsi::utils::BodyForceDiscreteConstant<Element,VecVecDim> PCBodyForce;

    typedef typename PCBodyForce::LoadFun                           LoadFun;
    PCBodyForce pcbf( std::bind( &Element::template bodyForce<LoadFun>,
                                 std::placeholders::_1, std::placeholders::_2, 
                                 std::placeholders::_3, std::placeholders::_4, 1. ),
                      elementTractions );

    driver.setExternalTractionFun( std::bind( pcbf, std::placeholders::_1, 
                                              std::placeholders::_2, std::placeholders::_3 ) );
}

//------------------------------------------------------------------------------
/** Return the vector-norm of the difference between previous and current increments
 *  \f[
 *      \| (x_{n+1}^{(k+1)} - x_{n}) - (x_{n+1}^{(k)} - x_{n}) \|_2
 *  \f]
 */
template< typename DRIVER >
double gshell::fem::Coupler<DRIVER>::
incrementResidual( Mesh & mesh,
                   const std::vector<VecDim> & nodeIncrementsPrev )
{
    // residual displacements (increment of displacement increments)
    const double normedResIncrSq =
        std::inner_product( mesh.nodesBegin(),
                            mesh.nodesEnd(),
                            nodeIncrementsPrev.begin(), 0.,
                            std::plus<double>(),
                            std::bind( &gshell::fem::norm2Squared<VecDim>,
                                         std::bind( std::minus<VecDim>(),
                                                      std::bind( &Node::giveIncrement, 
                                                                 std::placeholders::_1 ), 
                                                                 std::placeholders::_2 ) ) );
                        
    return std::sqrt( normedResIncrSq );
}


//------------------------------------------------------------------------------
// Update scaled increments by #relax factor
template< typename DRIVER >
double gshell::fem::Coupler<DRIVER>::
updateIncrements( Mesh & mesh,
                  std::vector<VecDim> & nodeIncrementsPrev,
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
                                             std::bind( gshell::fem::Multiply<VecDim,double,VecDim>(),
                                                          relax,
                                                          std::bind( &Node::giveIncrement, std::placeholders::_1 ) ),
                                             std::bind( gshell::fem::Multiply<VecDim,double,VecDim>(),
                                                          1.-relax,
                                                          std::placeholders::_2 ) ) ) );
    typename std::vector<VecDim>::const_iterator nodeIncrIter = nodeIncrementsPrev.begin();
    mesh.iterateOverNodes( std::bind( scaleIncr, std::placeholders::_1, std::ref( nodeIncrIter ) ) );
    
    // residual displacements (increment of displacement increments)
    const double normedResIncrSq =
        std::inner_product( mesh.nodesBegin(),
                            mesh.nodesEnd(),
                            nodeIncrementsPrev.begin(), 0.,
                            std::plus<double>(),
                            std::bind( &gshell::fem::norm2Squared<VecDim>,
                                         std::bind( std::minus<VecDim>(),
                                                      std::bind( &Node::giveIncrement, 
                                                                 std::placeholders::_1 ), 
                                                                 std::placeholders::_2 ) ) );
                        
    // store increment
    std::transform( mesh.nodesBegin(),
                    mesh.nodesEnd(),
                    nodeIncrementsPrev.begin(),
                    std::bind( &Node::giveIncrement, std::placeholders::_1 ) );

    // return norm
    return std::sqrt( normedResIncrSq );
}
