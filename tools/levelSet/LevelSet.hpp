// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//
//                 Computational Structural Mechanics Lab
//                         University of Cambridge
//
// This software is copyrighted and all rights are retained by the University of
// Cambridge. All use, disclosure, and/or reproduction of any part not expressly
// authorised by the University of Cambridge is prohibited. (C) 2010
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

//! @file LevelSet.hpp

#ifndef tools_levelset_h
#define tools_levelset_h

// standard library includes
#include <iostream>
#include <fstream>
#include <iterator>
#include <vector>
#include <limits>
#include <algorithm>
#include <cmath>
#include <utility>
#include <bitset>
// for use with floating point comparison
#include <float.h>
#include <math.h>

// boost includes
#include <boost/algorithm/string.hpp>
#include <boost/tuple/tuple.hpp>
#include <boost/function.hpp>
#include <boost/bind.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/array.hpp>
// #include <boost/timer.hpp>

// corlib includes
#include <corlib/verify.hpp>
#include <corlib/fuzzyEqual.hpp>
#include <corlib/linalg.hpp>

namespace tools {
    namespace lset {
        namespace ublas = boost::numeric::ublas;

        template<unsigned DIM> class LevelSet;
        template<unsigned DIM> class Point;

        namespace detail_ {

            // Computational Geometry. These functions are implemented as
            // shown in: 'Real-Time Collision Detection' by Christer Ericsion
            //------------------------------------------------------------------
            /** Return the closest point on the given simplex to the given point */
            template<unsigned DIM>
            ublas::bounded_vector<double,DIM>
            closestPointToSimplex( const ublas::bounded_vector<double,DIM> &,
                                   const boost::array<ublas::bounded_vector<double,DIM>,DIM> &);

            // \cond SKIPDOX
            // specialization for line-segment
            template<>
            ublas::bounded_vector<double,2>
            closestPointToSimplex<2>( const ublas::bounded_vector<double,2> &,
                                      const boost::array<ublas::bounded_vector<double,2>,2> &);
            /// specialization for triangle
            template<>
            ublas::bounded_vector<double,3>
            closestPointToSimplex<3>( const ublas::bounded_vector<double,3> &,
                                      const boost::array<ublas::bounded_vector<double,3>,3> &);
            // \endcond

            /** Orient points w.r.t simplices using determinant predicates, return determinant */
            template<unsigned DIM>
            double
            orientPointWithSimplex( const ublas::bounded_vector<double,DIM> &,
                                    const boost::array<ublas::bounded_vector<double,DIM>,DIM> &);
            // \cond SKIPDOX
            // specialization for line segment
            template<>
            double
            orientPointWithSimplex<2>( const ublas::bounded_vector<double,2> &,
                                       const boost::array<ublas::bounded_vector<double,2>,2> &);
            // specialization for triangle
            template<>
            double
            orientPointWithSimplex<3>( const ublas::bounded_vector<double,3> &,
                                       const boost::array<ublas::bounded_vector<double,3>,3> &);
            // \endcond

            // Other helper functions
            //------------------------------------------------------------------
            /** Compare two coordinates in a given direction */
            template<unsigned DIM>
            bool compareCoords(const ublas::bounded_vector<double, DIM> & c0,
                               const ublas::bounded_vector<double, DIM> & c1,
                               const unsigned dir)
            {return c0(dir) < c1(dir);}

            /** Compares the coordinates of two points in a given direction */
            template<unsigned DIM>
            bool comparePointCoords(const Point<DIM> &,
                                    const Point<DIM> &,
                                    const unsigned);

            /** Check if given coordinate lies in given bounding box */
            template<unsigned DIM>
            bool withinBBox(const ublas::bounded_vector<double, DIM> &,
                            const boost::array<ublas::bounded_vector<double,DIM>,2> &);

            /** Return incremented axis-aligned bounding-box for given point cloud */
            template <unsigned DIM, typename ITER>
            boost::array<ublas::bounded_vector<double,DIM>, 2>
            computeAABB(const ITER first, const ITER last, const double increment = 0.0);

            /** Compute the sign of a number */
            int sgn(double a) {return a < 0.0 ? -1 : 1;}

            /**
             * \brief A robust floating point comparison
             * \detail See C.Ericsion 'Real-Time Collision Detection' Chap. 11 'Numerical Robustness'.
             * This test:
             * - is better than a relative test: works well if a, b < 1.0
             * - is better than an absolute test: distance between adjacent numbers depends on size
             * - has a good choice for tol: approx half the available digits of prescision will match
             */
            bool combinedFuzzyEqual(double a, double b, double tol=sqrt(DBL_EPSILON))
            {
                return std::fabs(a - b) <= tol * std::max(std::max(std::fabs(a), std::fabs(b)), 1.0);
                // faster alternative
                //return std::fabs(a - b) <= tol * (std::fabs(a) + std::fabs(b) + 1.0);}
            }
        } // detail_
    } // lset
} // tools

//------------------------------------------------------------------------------
template<unsigned DIM>
class tools::lset::LevelSet
{
public:
    static const unsigned dim     = DIM;
    static const unsigned embdim  = dim    - 1;
    static const unsigned vpemesh = embdim + 1;

    typedef ublas::bounded_vector<double, dim> VecDim;
    typedef boost::array<unsigned, vpemesh>    Simplex;
    typedef std::vector<Simplex>               VecSimplex;
    typedef std::vector<VecDim>                VecVecDim;

    // mode how Level Set is computed
    enum mode {
        DEFAULT,    //!< The default mode is set to BRUTEFORCE
        BRUTEFORCE, //!< This algorithm iterates over all the background mesh points
                    //!< and then over all the embedded mesh elements computing the CPT
        BBOX,       //!< This algorithm iterates over the embedded mesh elements and
                    //!< computes the CPT only in a (inflated) AABB around each element
    };

private:
    // container typeefs
    typedef ublas::bounded_vector<int, dim> VecIntDim_;
    typedef std::vector<Point<dim> >        VecPoints_;
    typedef boost::array<VecDim, dim>       ArrVecDim_;

    // for legibility
    typedef boost::array<VecDim, 2> Box_;

    // iterator typedefs
    typedef typename std::vector<int>::const_iterator VecIntCIter_;
    typedef typename VecPoints_::iterator             VecPointsIter_;
    typedef typename VecPoints_::const_iterator       VecPointsCIter_;

public:
    /**
     * Constructor given the maximal distance of interest and optionally a
     * verbosity flag.
     *
     * \param[in] inflationLength Length increase to apply to the AABBs
     *                            This should be in the order of the largest cell size
     * \param[in] verbose         Verbosity Flag
     */
    LevelSet(double inflationLength,
             bool verbose = false)
        : AABBinc_(inflationLength),
          verbose_(verbose)
    {
        this -> initialise();
    }

    //! Set all members to an initial state
    void initialise() 
    {
        coordsMesh_.clear();
        connecMesh_.clear();
        gridPoints_.clear();
        for (unsigned d = 0; d < dim; ++d) numPoints_[d] = 0;
        latticeSet_ = false;
        meshSet_    = false;
    }

    /**
     * Compute the signed closest point transform using the specified method
     * \param[in] onesided  Flag indicating if level set is signed
     * \param[in] mode      Method used to compute the levelset
     */
    void perform(const bool onesided,
                 const enum mode mode = DEFAULT);

    /** Set the representation of the immersed surface */
    void setSurface( const VecVecDim  & coordinates,
                     const VecSimplex & connectivity);

    /** Set the representation of the background grid, sgf format is assumed */
    void setLattice(const int numPointsX,
                    const int numPointsY,
                    const int numPointsZ,
                    const VecVecDim & coordinates);

    /** Pass level set from the grid points to the given container */
    void getLevelSet(std::vector<double> & levelSet) const;

    /** Pass closest point vectors the to given container */
    void getCPT(VecVecDim & cpt) const;

    /** Pass closest element ID of grid points to given container */
    void getClosestElementIds(std::vector<int> & closestElements) const;

private:
    // Compute CPT for all points in mesh bounding box
    void computeCPTWithBruteForce_(const std::vector<int> & pointIdsWithinMeshAABB,
                                   const bool onesided);

    // Compute CPT for points within element bounding boxes to reduce work
    void computeCPTWithElementBBoxes_(const std::vector<int> & pointIdsWithinMeshAABB,
                                      const bool onesided);

    // Extract vertex coordinates for an element with the specified id
    ArrVecDim_ extractSimplex_(const unsigned) const;

    // Assign a dummy level set to points not computed in perform()
    void floodFill_(const bool onesided);

private:
    // LevelSet parameters
    const double AABBinc_;  //!< AABB inflation factor
    const bool   verbose_;  //!< Verbosity flag

    // Embedded mesh and grid
    VecVecDim  coordsMesh_;    //!< coordinates of the embedded mesh nodes
    VecSimplex connecMesh_;    //!< connectivity of the embedded mesh

    // Grid coordinates as points
    VecPoints_ gridPoints_;                   //!< Vector of grid points which store the cpt dat
    VecIntDim_ numPoints_;                    //!< Number of grid points in each direction

    // Sanity checking
    bool latticeSet_;
    bool meshSet_;
};

//------------------------------------------------------------------------------
/**
 * The Point class has a coordinate, an id and stores all the CPT infomation (i.e.
 * the distance from the embedded mesh, the sign of the distance, the coordinates
 * of the closest point on the embedded mesh in. The Point is also copy-constructable.
 */
template<unsigned DIM>
class tools::lset::Point
{
public:
    typedef ublas::bounded_vector<double, DIM> VecDim;

    /** Default constructor */
    Point()
        : coord_(ublas::zero_vector<double>(DIM)),
          id_(-1)
    {
        this -> initialise();
    }

    /** Constructor given a coordinate and id */
    Point(VecDim coordIn, unsigned idIn)
        : coord_(coordIn),
          id_(idIn)
    {
        this -> initialise();
    }

    //! Make sure, the private data are invalidated
    void initialise()
    {
        active_           = false;
        distance_         = std::numeric_limits<double>::max();
        orientDet_        = std::numeric_limits<double>::max();
        closestElementId_ = -1; 
        closestPoint_     = std::numeric_limits<double>::max() * 
            ublas::scalar_vector<double>(DIM);
    }

    /** Update the CPT and set the activity flag */
    void updateCPT(const double dist,
                   const int elementId,
                   const VecDim & point,
                   const double orientDet);

    /** Set sign of this point. This is used in the LevelSet flood fill method */
    void setSign(const double sign) {orientDet_ = sign;}

    void setDistance(const double dist) {distance_ = dist;}

    //! @name Access data members
    //@{
    VecDim getCoordinates()        const {return coord_;}
    int    getId()                 const {return id_;}
    bool   isActive()              const {return active_;}
    double getDistance()           const {return distance_;}
    double getSign()               const {return detail_::sgn(orientDet_);}
    int    getClosestElementId()   const {return closestElementId_;}
    /** Position vector of the closest point */
    VecDim getClosestPoint()       const {return closestPoint_;}
    /** Relative position vector of closest point i.e. closestPoint - point */
    VecDim getCPT() const;
    double getLevelSet()           const {return detail_::sgn(orientDet_) * distance_;}
    //@}

private:
    VecDim coord_;
    int    id_;
    bool   active_; //!< Indicates if the point has been used in the CPT computation
    double distance_;
    double orientDet_;
    int    closestElementId_;
    VecDim closestPoint_;
};

#include "LevelSet.ipp"
#endif
