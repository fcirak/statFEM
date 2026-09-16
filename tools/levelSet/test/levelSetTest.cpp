// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//
//                 Computational Structural Mechanics Lab
//                         University of Cambridge
//
// This software is copyrighted and all rights are retained by the University of
// Cambridge. All use, disclosure, and/or reproduction of any part not expressly
// authorised by the University of Cambridge is prohibited. (C) 2010
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

//! @file testLevelSet.hpp

//------------------------------------------------------------------------------
//! System includes
#include <fstream>
#include <string>

//! Boost includes
#include <boost/lexical_cast.hpp>

//! Corlib includes
#include <corlib/verify.hpp>
#include <corlib/SmfHead.hpp>
#include <corlib/DistributeQuantity.hpp>
#include <corlib/CollectQuantity.hpp>

//! Level set tool
#include <tools/levelSet/LevelSet.hpp>

//! Testing functions
#include "levelSetTest.hpp"

//------------------------------------------------------------------------------
//! Deduce grid dimension from given sgf input stream
unsigned readDimFromSGF(std::istream & sgf)
{
    corlib::SmfHead sgfHead;
    sgfHead.read(sgf);
    const unsigned numPointsPerCell = sgfHead.giveElementNumPoints();

    sgf.clear();
    sgf.seekg(0);

    if      (numPointsPerCell == 2) return 1;
    else if (numPointsPerCell == 4) return 2;
    else if (numPointsPerCell == 8) return 3;
    else FTL_VERIFY(false);
    return 0;
}

//------------------------------------------------------------------------------
//! Deduce surface mesh embedding dimension from given smf input stream
unsigned readDimFromSMF(std::istream & smf)
{
    corlib::SmfHead smfHead;
    smfHead.read(smf);
    const corlib::shape elemShape = smfHead.giveElementShape();

    smf.clear();
    smf.seekg(0);

    if      (elemShape == corlib::LINE)     return 2;
    else if (elemShape == corlib::TRIANGLE) return 3;
    else FTL_VERIFY(false);
    return 0;
}

/** Generate grid and surface, compute level set, write VTS file */
template<unsigned DIM, unsigned DEGREE>
struct GenerateCompute
{
    static void apply(std::istream & sgf,
                      std::istream & smf,
                      std::ostream & vts,
                      const bool isUnsigned,
                      const double maxLevelSet,
                      const bool bruteForce)
    {
        // construct Grid
        typedef lset::Point<DIM>  Point;
        typedef lset::Grid<Point> Grid;
        Grid grid(sgf);

        // construct Surface
        typedef lset::Surface<DIM> Surface;
        Surface surface(smf);

        // constuct LevelSet
        const bool verbose = true;
        tools::lset::LevelSet<DIM> levelSet(maxLevelSet, verbose);

        // compute level set
        {
            // collect grid coordinates and set grid
            std::vector<typename Grid::VecDim> gridCoordinates;
            corlib::CollectQuantity<Point,typename Grid::VecDim>
                collectGridCoordinates(boost::bind(&Point::getCoordinates, _1));
            grid.iterateOverPoints(boost::bind(collectGridCoordinates, _1,
                                               std::back_inserter(gridCoordinates)));

            boost::numeric::ublas::bounded_vector<unsigned,3> numCells;
            numCells = grid.getNumCells();

            boost::numeric::ublas::bounded_vector<unsigned,3> numPoints;
            numPoints[0] = numCells[0] + DEGREE;
            numPoints[1] = numCells[1] + DEGREE;
            numPoints[2] = DIM > 2 ? numCells[2] + DEGREE: 1;
            levelSet.setLattice(numPoints[0], numPoints[1], numPoints[2], gridCoordinates);

            // collect surface coordinates and connectivity and set surface
            std::vector<typename Surface::VecDim>
                surfaceCoordinates(surface.coordinatesBegin(), surface.coordinatesEnd());
            std::vector<typename Surface::IndexSimplex>
                surfaceConnectivity(surface.connectivityBegin(), surface.connectivityEnd());
            levelSet.setSurface(surfaceCoordinates, surfaceConnectivity);

            // compute level set
            if (bruteForce)
                levelSet.perform(isUnsigned, tools::lset::LevelSet<DIM>::BRUTEFORCE);
            else
                levelSet.perform(isUnsigned, tools::lset::LevelSet<DIM>::BBOX);
        }

        // apply to grid
        {
            // retrieve distance values from level set object
            std::vector<double> distances;
            levelSet.getLevelSet(distances);
            corlib::DistributeQuantity<Point,double> passLevelSet(&Point::setLevelSet);
            grid.iterateOverPoints(boost::bind(passLevelSet, _1, distances.begin()));
        }

        {
            // retrieve distance values from level set object
            std::vector<typename Grid::VecDim> cps;
            levelSet.getCPT(cps);
            corlib::DistributeQuantity<Point,typename Grid::VecDim> passCP(&Point::setCP);
            grid.iterateOverPoints(boost::bind(passCP, _1, cps.begin()));
        }

        {
            // retrieve distance values from level set object
            std::vector<int> faceIds;
            levelSet.getClosestElementIds(faceIds);
            corlib::DistributeQuantity<Point,int> passFaceId(&Point::setFaceId);
            grid.iterateOverPoints(boost::bind(passFaceId, _1, faceIds.begin()));
        }

        // write VTS file
        grid.writeVTS(vts);
    }
};

//------------------------------------------------------------------------------
int main(int argc, char * argv[])
{
    if(argc != 3 and argc != 4 and argc != 5 and argc != 6){
        std::cout << "Usage " << argv[0]
                  << " Grid.sgf Surface.smf BruteForceFlag(0,1) UnsignedFlag(0,1) MaxLevelSet"
                  << std::endl
                  << std::endl;
        exit(0);
    }

    // read grid file
    const std::string gridFile = boost::lexical_cast<std::string>(argv[1]);
    std::ifstream sgf(gridFile.c_str());
    FTL_VERIFY(sgf.is_open());
    const unsigned sgfDim = readDimFromSGF(sgf);

    // read surface file
    const std::string surfaceFile = boost::lexical_cast<std::string>(argv[2]);
    std::ifstream smf(surfaceFile.c_str());
    FTL_VERIFY(smf.is_open());
    const unsigned smfDim = readDimFromSMF(smf);

    // read other options
    const bool bruteForce = boost::lexical_cast<bool>(argv[3]);
    const bool isUnsigned = boost::lexical_cast<bool>(argv[4]);
    const double maxLevelSet = boost::lexical_cast<double>(argv[5]);
    FTL_VERIFY( sgfDim == smfDim );

    // open vts file for grid output
    const std::string gridName = gridFile.substr(0, gridFile.find(".sgf"));
    const std::string surfName = surfaceFile.substr(0, surfaceFile.find(".smf"));
    const std::string modeName = bruteForce ? "BF" : "EBB";
    const std::string vtsFile = gridName + "_" + surfName + "_" + modeName + ".vts";
    std::ofstream vts(vtsFile.c_str());

    if (sgfDim == 2)
        GenerateCompute<2,1>::apply(sgf, smf, vts, isUnsigned, maxLevelSet, bruteForce);
    else if (sgfDim == 3)
        GenerateCompute<3,1>::apply(sgf, smf, vts, isUnsigned, maxLevelSet, bruteForce);
    else FTL_VERIFY(false);

    sgf.close();
    smf.close();
    vts.close();

    return 0;
}
