// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//
//                 Computational Structural Mechanics Lab
//                         University of Cambridge
//
// This software is copyrighted and all rights are retained by the University of
// Cambridge. All use, disclosure, and/or reproduction of any part not expressly
// authorised by the University of Cambridge is prohibited. (C) 2010
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

//! @file LevelSet.ipp

//------------------------------------------------------------------------------
template<unsigned DIM>
void tools::lset::LevelSet<DIM>::
perform(const bool onesided,
        const enum mode mode)
{
    //! Reset the states of the grid points
    std::for_each( gridPoints_.begin(), gridPoints_.end(), 
                   boost::bind( &Point<dim>::initialise, _1 ) );

    if(verbose_) std::cout << "perform_()" << std::endl;

    // preconditions
    FTL_VERIFY_DESCRIPTIVE(latticeSet_, "Background grid not set.\n");
    FTL_VERIFY_DESCRIPTIVE(meshSet_, "Embedded mesh not set.\n");

    // compute AABB for the whole mesh
    const Box_ meshAABB =
        detail_::computeAABB<dim, typename VecVecDim::iterator>
        (coordsMesh_.begin(), coordsMesh_.end(), AABBinc_);

    // find ids of grid points that lie in the mesh AABB
    std::vector<int> pointIdsWithinMeshAABB;
    VecPointsCIter_ begin = gridPoints_.begin();
    VecPointsCIter_ end   = gridPoints_.end();
    for (VecPointsCIter_ iter = begin; iter != end; ++ iter) {
        if(detail_::withinBBox<dim>(iter->getCoordinates(), meshAABB))
            pointIdsWithinMeshAABB.push_back(iter->getId());
    }

    // compute the CPT
    if ((mode == BRUTEFORCE) or (mode == DEFAULT)) {
        computeCPTWithBruteForce_(pointIdsWithinMeshAABB, onesided);
    } else if (mode == BBOX) {
        computeCPTWithElementBBoxes_(pointIdsWithinMeshAABB, onesided);
    } else {
        FTL_VERIFY_DESCRIPTIVE(false, "Cannot handle chosen CPT mode.\n");
    }

    // fill unused points with a dummy CPT
    floodFill_(onesided);

    return;
}

//------------------------------------------------------------------------------
template<unsigned DIM>
void tools::lset::LevelSet<DIM>::
setSurface( const VecVecDim & coordinates, const VecSimplex & connectivity )
{

    if (verbose_) {
        std::cout << "setSurface_()"  << std::endl;
        std::cout << "  nNodes = "    << coordinates.size()  << std::endl;
        std::cout << "  nElements = " << connectivity.size() << std::endl;
    }

    coordsMesh_    = coordinates;
    connecMesh_    = connectivity;

    meshSet_ = true;

    return;
}

//------------------------------------------------------------------------------
template<unsigned DIM>
void tools::lset::LevelSet<DIM>::
setLattice(const int numPointsX,
           const int numPointsY,
           const int numPointsZ,
           const VecVecDim & coordinates)
{
    if (verbose_) std::cout << "setLattice_()" << std::endl;

    // check dimensions
    if (dim < 3) FTL_VERIFY(numPointsZ == 1);
    if (dim < 2) FTL_VERIFY(numPointsY == 1);
    FTL_VERIFY(numPointsX);

    // store number of cells and points per direction
    if (dim > 2) numPoints_[2] = numPointsZ;
    if (dim > 1) numPoints_[1] = numPointsY;
    numPoints_[0] = numPointsX;

    // check number of points
    std::size_t product = 1;
    for (unsigned d = 0; d < dim; ++d) product *= numPoints_[d];
    FTL_VERIFY(product == coordinates.size());
    
    // convert raw coordinates into vector of grid points
    const unsigned numPoints = static_cast<unsigned>(coordinates.size());
    gridPoints_.resize(numPoints);
    for (unsigned i = 0; i < numPoints; i ++) {
        gridPoints_[i] = Point<dim>(coordinates[i], i);
    }

    latticeSet_ = true;

    return;
}

//------------------------------------------------------------------------------
template<unsigned DIM>
void tools::lset::LevelSet<DIM>::
getLevelSet(std::vector<double> & levelSet) const
{
    if (verbose_) std::cout << "getLevelSet_()" << std::endl;

    levelSet.resize(gridPoints_.size());
    std::transform(gridPoints_.begin(), gridPoints_.end(), levelSet.begin(),
                   boost::bind(&Point<dim>::getLevelSet, _1));
}

//------------------------------------------------------------------------------
template<unsigned DIM>
void tools::lset::LevelSet<DIM>::
getCPT(VecVecDim & cpt) const
{
    if (verbose_) std::cout << "getCPT_()" << std::endl;

    cpt.resize(gridPoints_.size());
    std::transform(gridPoints_.begin(), gridPoints_.end(), cpt.begin(),
                   boost::bind(&Point<dim>::getCPT, _1));
}

//------------------------------------------------------------------------------
template<unsigned DIM>
void tools::lset::LevelSet<DIM>::
getClosestElementIds(std::vector<int> & closestElements) const
{
    if (verbose_) std::cout << "getClosestElementIds_()" << std::endl;

    closestElements.resize(gridPoints_.size());
    std::transform(gridPoints_.begin(), gridPoints_.end(),
                   closestElements.begin(),
                   boost::bind(&Point<dim>::getClosestElementId, _1));
}

//------------------------------------------------------------------------------
/**
 * \brief Implementation of closest point transform using BRUTE FORCE.
 *
 * For the points whose point ids are stored in the pointIdsWithinMeshAABB
 * vector, the method computes the cpt data for each element of the embedded
 * mesh. This method should have time complexity O(MN) where M is the number of
 * _elements_ in the embedded mesh and N is the number _nodes_ in the embedding
 * grid.
 *
 * \sa computeCPTWithElementBBoxes_()
 */
template<unsigned DIM>
void tools::lset::LevelSet<DIM>::
computeCPTWithBruteForce_(const std::vector<int> & pointIdsWithinMeshAABB,
                          const bool onesided)
{
    if (verbose_) {
        std::cout << "computeCPTWithBruteForce_()" << std::endl;
        std::cout << "  " << pointIdsWithinMeshAABB.size()
                  << " points within incremented Mesh AABB" << std::endl;
    }

    VecIntCIter_ pointIdBegin = pointIdsWithinMeshAABB.begin();
    VecIntCIter_ pointIdEnd   = pointIdsWithinMeshAABB.end();
    for (VecIntCIter_ idp = pointIdBegin; idp != pointIdEnd; ++ idp){
        for (unsigned e = 0; e < connecMesh_.size(); e ++){

            // extract element e
            const ArrVecDim_ verts = this -> extractSimplex_(e);

            // compute closest point
            const VecDim closestPoint =
                detail_::closestPointToSimplex<dim>(gridPoints_[*idp].getCoordinates(), verts);

            // compute distance to closest point
            const double dist =
                ublas::norm_2(gridPoints_[*idp].getCoordinates() - closestPoint);
            // compute orientation if signed distance is required
            double orientDet = 1.0;
            if(!onesided)
                orientDet = detail_::orientPointWithSimplex<dim>(gridPoints_[*idp].getCoordinates(), verts);
            // update the cpt tuple of the point in the global points vector
            gridPoints_[*idp].updateCPT(dist, e, closestPoint, orientDet);
        }
    }

    return;
}

//------------------------------------------------------------------------------
/**
 * \brief Implementation of the closest point transform using an ELEMENT
 * BOUNDING BOX approach.
 *
 * First the grid points that are within in the mesh AABB are copied into
 * DIM vectors which are sorted in each coordinate direction in turn. For
 * each element in the embedded mesh an inflated AABB is computed, the cpt
 * data for points in these boxes only is computed. The points to be tested
 * are found using a binary search of the sorted arrays to find slices whose
 * intersection is the element AABB. The slice containing the fewest points
 * is found and we iterate over these points computing the cpt data for the
 * points which lie in the element AABB. The complexity of this method is...
 *
 * \sa computeCPTWithBruteForce_()
 */
template<unsigned DIM>
void tools::lset::LevelSet<DIM>::
computeCPTWithElementBBoxes_(const std::vector<int> & pointIdsWithinMeshAABB,
                             const bool onesided)
{
    if(verbose_) std::cout << "computeCPTWithElementBBoxes_()" << std::endl;

    // extract DIM copies of the grid points within the mesh AABB and sort in each direction
    boost::array<VecPoints_, dim> gridPointsSorted;
    for (unsigned d = 0; d < dim; ++ d) {
        gridPointsSorted[d].resize(pointIdsWithinMeshAABB.size());
        for (unsigned i = 0; i < pointIdsWithinMeshAABB.size(); ++ i)
            gridPointsSorted[d][i] = gridPoints_[ pointIdsWithinMeshAABB[i] ];

        std::sort(gridPointsSorted[d].begin(),
                  gridPointsSorted[d].end(),
                  boost::bind(detail_::comparePointCoords<dim>, _1, _2, d));
    }

    // iterate over the mesh elements
    for (unsigned e = 0; e < connecMesh_.size(); ++e) {

        // extract element e
        const ArrVecDim_ verts = this -> extractSimplex_(e);

        // compute (inflated) element AABB
        Box_ bbox =
            detail_::computeAABB<dim, typename ArrVecDim_::const_iterator>
            (verts.begin(), verts.end(), AABBinc_);
        Point<dim> low(bbox[0], 0);
        Point<dim> upp(bbox[1], 0);

        // identify points in the ranges of interest and find dim with fewest points
        boost::array<VecPointsIter_, dim> itup, itlo;
        unsigned dimMinPoints = std::numeric_limits<unsigned>::max();
        for (unsigned d = 0; d < dim; ++ d) {
            itup[d] = std::upper_bound(gridPointsSorted[d].begin(),
                                       gridPointsSorted[d].end(),
                                       upp, boost::bind(detail_::comparePointCoords<dim>, _1, _2, d));
            itlo[d] = std::lower_bound(gridPointsSorted[d].begin(),
                                       gridPointsSorted[d].end(),
                                       low, boost::bind(detail_::comparePointCoords<dim>, _1, _2, d));
            unsigned numPoints = static_cast<unsigned>(itup[d] - itlo[d]);
            if (numPoints < dimMinPoints) dimMinPoints = d;
        }

        // compute LevelSet for range of interest
        VecPointsIter_ testStart = itlo[dimMinPoints];
        VecPointsIter_ testEnd   = itup[dimMinPoints];
        for (VecPointsIter_ iter = testStart; iter != testEnd; ++ iter) {
            // check if point is within the inflated element AABB
            if(detail_::withinBBox<dim>(iter->getCoordinates(), bbox)) {

                // compute closest point
                const VecDim closestPoint =
                    detail_::closestPointToSimplex<dim>(iter->getCoordinates(), verts);

                // compute distance to closest point
                const double dist =
                    ublas::norm_2(iter->getCoordinates() - closestPoint);

                // compute orientation if signed distance is required
                double orientDet = 1.0;
                if(!onesided)
                    orientDet = detail_::orientPointWithSimplex<dim>(iter->getCoordinates(), verts);

                // update the cpt tuple of the point in the global points vector
                gridPoints_[iter->getId()].updateCPT(dist, e, closestPoint, orientDet);
            }
        }
    }

    return;
}

//------------------------------------------------------------------------------
/** Extract vertex coordinates for an element with the specified id */
template<unsigned DIM>
typename tools::lset::LevelSet<DIM>::ArrVecDim_
tools::lset::LevelSet<DIM>::
extractSimplex_(const unsigned e) const
{
    // vector of element vertex coordinates
    ArrVecDim_ verts;

    // loop over verticies of element
    for (unsigned v = 0; v < vpemesh; ++ v) {
        // get global vertex index
        const unsigned i = connecMesh_[e][v];

        // insert coordinate of element's vertex
        verts[v] = coordsMesh_[i]; 
    }

    return verts;
}

//------------------------------------------------------------------------------
/** 
 * Fill the values of the level set at unused points (i.e. point.active_ ==
 * false) with the maximum computed distance and set the sign of the distance
 */
template<unsigned DIM>
void tools::lset::LevelSet<DIM>::
floodFill_(const bool onesided)
{
    if (verbose_) std::cout << "floodFill_()" << std::endl;

    // find maximum distance function amongst active points
    double maxComputedDist = 0.0;
    VecPointsIter_ begin = gridPoints_.begin();
    VecPointsIter_ end = gridPoints_.end();
    for (VecPointsCIter_ iter = begin; iter != end; ++iter) {
        if (iter->isActive() and (iter->getDistance() > maxComputedDist))
            maxComputedDist = iter->getDistance();
    }
    if (verbose_) 
        std::cout << "  Maximum distance amongst points computed"
                  << " in perform() method: " << maxComputedDist << std::endl;

    // loop over seed points
    unsigned numSeeds = 1;
    for (unsigned d = 1; d < dim; ++d) numSeeds *= numPoints_[d];
    for (unsigned s = 0; s < numSeeds; ++s) {
        // get point index of seed point
        const unsigned idx = s * numPoints_[0];

        // is seed point already computed?
        const bool seedComputed = gridPoints_[idx].isActive();

        // if needed find the cpt for the seed point
        if (not seedComputed) {
            VecDim coords = gridPoints_[idx].getCoordinates();
            for (unsigned e = 0; e < connecMesh_.size(); e ++) {
                ArrVecDim_ verts     = this -> extractSimplex_(e);
                VecDim closestPoint = detail_::closestPointToSimplex<dim>(coords, verts);
                double dist          = ublas::norm_2(coords - closestPoint);
                double orientDet = 1.0;
                if (!onesided) orientDet = detail_::orientPointWithSimplex<dim>(coords, verts);
                gridPoints_[idx].updateCPT(dist, e, closestPoint, orientDet);
            }
        }

        // flood fill points along x-direction from seed
        double prevSign = gridPoints_[idx].getSign();
        for (int i = 0; i < numPoints_[0]; ++i) {
            if (not gridPoints_[idx + i].isActive()) {
                (gridPoints_[idx + i]).setDistance(maxComputedDist);
                (gridPoints_[idx + i]).setSign(prevSign);
            } else {
                prevSign = (gridPoints_[idx + i]).getSign();
            }
        }

        // clean up artefacts of filling process
        if (not seedComputed) {
            double seedSign = gridPoints_[idx].getSign();
            gridPoints_[idx].initialise();
            gridPoints_[idx].setDistance(maxComputedDist);
            gridPoints_[idx].setSign(seedSign);
        }
    }

    return;
}

//------------------------------------------------------------------------------
template<unsigned DIM>
typename tools::lset::Point<DIM>::VecDim
tools::lset::Point<DIM>::getCPT() const
{
    if(active_)
        return closestPoint_;
    else{
        VecDim zeroVector;
        zeroVector.clear();
        return zeroVector;
    }
}

//------------------------------------------------------------------------------
template<unsigned DIM>
void tools::lset::Point<DIM>::
updateCPT(const double dist,
          const int elementId,
          const VecDim & point,
          const double orientDet)
{
    // set activity flag since point has been used
    active_ = true;

    if (detail_::combinedFuzzyEqual(dist, distance_)) {
        // if distances match chose the face with larger orientation determinant
        distance_         = dist;
        closestElementId_ = std::fabs(orientDet) > std::fabs(orientDet_) ? elementId : closestElementId_;
        closestPoint_     = point;
        orientDet_        = std::fabs(orientDet) > std::fabs(orientDet_) ? orientDet : orientDet_;
    } else if (dist < distance_) {
        distance_         = dist;
        closestElementId_ = elementId;
        closestPoint_     = point;
        orientDet_        = orientDet;
    }

    return;
}

//--------------------------------------------------------------------------
template<>
boost::numeric::ublas::bounded_vector<double,2>
tools::lset::detail_::closestPointToSimplex<2>(const ublas::bounded_vector<double,2>                 & point,
                                               const boost::array<ublas::bounded_vector<double,2>,2> & verts)
{
    const ublas::bounded_vector<double,2> ab = verts[1] - verts[0];
    // Project p onto ab, computing parameterized position d(t) = a + t*(b - a)
    double t = ublas::inner_prod(point - verts[0], ab) / ublas::inner_prod(ab, ab);
    // If outside segment, clamp t (and therefore d) to the closest endpoint
    if (t < 0.0) t = 0.0;
    if (t > 1.0) t = 1.0;
    // Compute projected position from the clamped t
    return verts[0] + t * ab;
}

//--------------------------------------------------------------------------
template<>
boost::numeric::ublas::bounded_vector<double,3>
tools::lset::detail_::closestPointToSimplex<3>(const ublas::bounded_vector<double, 3>                & point,
                                               const boost::array<ublas::bounded_vector<double,3>,3> & verts)
{
    // for convenience
    typedef ublas::bounded_vector<double, 3> VecDim;

    // Check if P in vertex region outside A
    const VecDim ab = verts[1] - verts[0];
    const VecDim ac = verts[2] - verts[0];
    const VecDim ap = point    - verts[0];
    const double d1 = ublas::inner_prod(ab, ap);
    const double d2 = ublas::inner_prod(ac, ap);
    if (d1 <= 0.0 and d2 <= 0.0)
        return verts[0]; // barycentric coordinates (1,0,0)

    // Check if P in vertex region outside B
    const VecDim bp = point - verts[1];
    const double d3 = ublas::inner_prod(ab, bp);
    const double d4 = ublas::inner_prod(ac, bp);
    if (d3 >= 0.0 and d4 <= d3)
        return verts[1]; // barycentric coordinates (0,1,0)

    // Check if P in edge region of AB, if so return projection of P onto AB
    double vc = d1 * d4 - d3 * d2;
    if (vc <= 0.0 and d1 >= 0.0 and d3 <= 0.0) {
        const double v = d1 / (d1 - d3);
        return verts[0] + v * ab; // barycentric coordinates (1-v,v,0)
    }

    // Check if P in vertex region outside C
    const VecDim cp = point - verts[2];
    const double d5 = ublas::inner_prod(ab, cp);
    const double d6 = ublas::inner_prod(ac, cp);
    if (d6 >= 0.0 and d5 <= d6)
        return verts[2]; // barycentric coordinates (0,0,1)

    // Check if P in edge region of AC, if so return projection of P onto AC
    const double vb = d5*d2 - d1*d6;
    if (vb <= 0.0 and d2 >= 0.0 and d6 <= 0.0) {
        const double w = d2 / (d2 - d6);
        return verts[0] + w * ac; // barycentric coordinates (1-w,0,w)
    }

    // Check if P in edge region of BC, if so return projection of P onto BC
    const double va = d3*d6 - d5*d4;
    if (va <= 0.0 and (d4 - d3) >= 0.0 and (d5 - d6) >= 0.0) {
        const double w = (d4 - d3) / ((d4 - d3) + (d5 - d6));
        return verts[1] + w * (verts[2] - verts[1]); // barycentric coordinates (0,1-w,w)
    }

    // P inside face region. Compute Q through its barycentric coordinates (u,v,w)
    const double denom = 1.0 / (va + vb + vc);
    const double v = vb * denom;
    const double w = vc * denom;
    return verts[0] + ab * v + ac * w; // = u*a + v*b + w*c, u = va * denom = 1.0 - v - w
}

//--------------------------------------------------------------------------
template<>
double
tools::lset::detail_::orientPointWithSimplex<2>(const ublas::bounded_vector<double,2>                 & point,
                                                const boost::array<ublas::bounded_vector<double,2>,2> & verts)
{
    // Compute which side p lies on
    const double det =
        (verts[0](0) - point(0)) * (verts[1](1) - point(1)) -
        (verts[1](0) - point(0)) * (verts[0](1) - point(1));

    return det;
}

//--------------------------------------------------------------------------
template<>
double
tools::lset::detail_::orientPointWithSimplex<3>(const ublas::bounded_vector<double,3>                 & point,
                                                const boost::array<ublas::bounded_vector<double,3>,3> & verts)
{
    // Compute which side p lies on
    typedef ublas::bounded_vector<double, 3> VecDim;
    const VecDim pa = verts[0] - point;
    const VecDim pb = verts[1] - point;
    const VecDim pc = verts[2] - point;
    const VecDim temp = corlib::cross_prod(pb, pc);
    const double det = ublas::inner_prod(pa, temp);

    return det;
}

//------------------------------------------------------------------------------
template<unsigned DIM>
bool tools::lset::detail_::
comparePointCoords(const Point<DIM> & p0,
                   const Point<DIM> & p1,
                   unsigned dir)
{
    ublas::bounded_vector<double,DIM> p0Coords, p1Coords;
    p0Coords = p0.getCoordinates();
    p1Coords = p1.getCoordinates();

    return p0Coords(dir) < p1Coords(dir);
}

//------------------------------------------------------------------------------
template<unsigned DIM>
bool tools::lset::detail_::
withinBBox(const ublas::bounded_vector<double, DIM>                & vec,
           const boost::array<ublas::bounded_vector<double,DIM>,2> & bbox)
{
    bool inBBox = true;
    for (unsigned d = 0; d < DIM; ++ d) {
        if(vec(d) < bbox[0](d) or vec(d) > bbox[1](d))
            inBBox = false;
    }
    
    return inBBox;
}

//------------------------------------------------------------------------------
template <unsigned DIM, typename ITER>
boost::array<boost::numeric::ublas::bounded_vector<double,DIM>, 2>
tools::lset::detail_::
computeAABB(const ITER first, const ITER last, const double increment)
{
    typedef boost::array<ublas::bounded_vector<double,DIM>,2> Box;
    typedef ublas::bounded_vector<double,DIM> VecDim;

    Box bbox;
    for (unsigned d = 0; d < DIM; ++ d) {
        const VecDim max = *std::max_element(first, last,
                                             boost::bind(compareCoords<DIM>, _1, _2, d));
        const VecDim min = *std::min_element(first, last,
                                             boost::bind(compareCoords<DIM>, _1, _2, d));
        bbox[0](d) = min(d) - increment;
        bbox[1](d) = max(d) + increment;
    }

    return bbox;
}
