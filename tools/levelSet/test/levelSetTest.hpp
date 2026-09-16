// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//
//                 Computational Structural Mechanics Lab
//                         University of Cambridge
//
// This software is copyrighted and all rights are retained by the University of
// Cambridge. All use, disclosure, and/or reproduction of any part not expressly
// authorised by the University of Cambridge is prohibited. (C) 2010
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

//! @file levelSetTest.hpp

#ifndef lset_levelsettest_h
#define lset_levelsettest_h

//! std includes
#include <istream>
#include <vector>

//! boost includes
#include <boost/array.hpp>
#include <boost/utility.hpp>
#include <boost/numeric/ublas/vector.hpp>

//! corlib includes
#include <corlib/verify.hpp>
#include <corlib/misc.hpp>
#include <corlib/Accessor.hpp>
#include <corlib/VTUwriter.hpp> // to access the free write functions
#include <corlib/SmfHead.hpp>

namespace lset {
    namespace ublas = boost::numeric::ublas;

    template<unsigned DIM> class Point;
    template<typename POINT> class Grid;
    template<unsigned DIM> class Surface;

    namespace detail_ {
        /** Construct VecInt3 from given components and padded with zeros */
        template<unsigned DIM>
        ublas::bounded_vector<int,3>
        makeMultiIndex_0(const int i1=0, const int i2=0, const int i3=0)
        {
            boost::array<int,3> aux = {{ i1, i2, i3}};
            ublas::bounded_vector<int,3> i;
            for (unsigned d = 0; d < 3; ++d) i[d] = 0;
            for (unsigned d = 0; d < DIM; ++d) i[d] = aux[d];
            return i;
        }

        /** Construct VecDim from Vec3 */
        template<unsigned DIM>
        ublas::bounded_vector<double,DIM>
        makeMultiIndex(const ublas::bounded_vector<double,3> & i)
        {
            ublas::bounded_vector<double,DIM> idx;
            for (unsigned d = 0; d < DIM; ++d) idx[d] = i[d];
            return idx;
        }
    }
}

//------------------------------------------------------------------------------
template <unsigned DIM>
class lset::Point
{
public:
    static const unsigned dim = DIM;
    typedef ublas::bounded_vector<double,dim> VecDim;
    typedef ublas::bounded_vector<int,dim>    VecIntDim;

    Point(VecIntDim idx, VecDim coords)
    {
        // invalidate point data
        coords_ = coords;
        index_  = idx;

        // invalidate other data
        levelSet_ = -std::numeric_limits<float>::max();
        faceId_   = std::numeric_limits<int>::max();
        cp_       = ublas::scalar_vector<double>(dim, std::numeric_limits<double>::max());
    }

    void setCoordinates(const VecDim &x) {coords_    = x;}
    void setId(const VecIntDim &idx)     {index_     = idx;}
    void setLevelSet(const double ls)    {levelSet_ = ls;}
    void setFaceId(const int f)          {faceId_   = f;}
    void setCP(const VecDim cp)          {cp_       = cp;}

    VecDim    getCoordinates() const {return coords_;}
    VecIntDim getId()          const {return index_;}
    double    getLevelSet()    const {return levelSet_;}
    int       getFaceId()      const {return faceId_;}
    VecDim    getCP()          const {return cp_;}

private:
    VecDim    coords_;   //!< Storage of coordinates
    VecIntDim index_;    //!< Storage of multi-index
    double    levelSet_; //!< Value of the control point's level set function
    int       faceId_;
    VecDim    cp_;
};

//------------------------------------------------------------------------------
template <typename POINT>
class lset::Grid : boost::noncopyable
{
public:
    static const unsigned dim = POINT::dim;

    typedef POINT Point;
    typedef typename std::vector<Point*> PointsVec;
    typedef typename Point::VecIntDim VecIntDim;

    typedef ublas::bounded_vector<int,3> VecInt3;
    typedef ublas::bounded_vector<double,3> Vec3;
    typedef ublas::bounded_vector<double,dim> VecDim;

public:
    Grid(std::istream &sgf)
    {
        this->readHeader_(sgf);
        this->readPoints_(sgf);
    };

    ~Grid()
    {
        this->iterateOverPoints(corlib::deleteFunctor());
    };

    template<typename OP>
    OP iterateOverPoints(OP op)
    {
        return corlib::for_each_if(points_.begin(), points_.end(),
                                   op, corlib::ValidPointer());
    };

    /** Write VTS file */
    void writeVTS(std::ostream & vts)
    {
        this->writeCoords_(vts);
        this->openPointData_(vts);
        this->writePointQuantity_(vts, corlib::accessorFun(&Point::getLevelSet, "LevelSet"));
        this->writePointQuantity_(vts, corlib::accessorFun(&Point::getFaceId, "FaceId"));
        this->writePointQuantity_(vts, corlib::accessorFun(&Point::getCP, "CP"));
        this->closePointData_(vts);
        this->finalize(vts);
    };

    VecInt3 getNumCells() {return numCells_;};

private:
    //! Functions to construct the grid
    void readHeader_(std::istream &);
    void readPoints_(std::istream &);

    //! Functions to write VTS file
    void writeCoords_(std::ostream &);
    void openPointData_(std::ostream &vts) {vts << "<PointData>" << std::endl;};
    template<typename PACC> 
    void writePointQuantity_(std::ostream &, PACC, const bool = false, const bool = false);
    void closePointData_(std::ostream &vts) {vts << "</PointData>" << std::endl;};
    void finalize(std::ostream &);

private:
    boost::array<int, 6> extents_;
    VecInt3 numCells_;
    PointsVec points_;
};

//------------------------------------------------------------------------------
template<unsigned DIM>
class lset::Surface : public boost::noncopyable
{
public:
    static const unsigned dim = DIM;
    static const unsigned numNodesPerElement = dim;

    typedef boost::numeric::ublas::bounded_vector<double,dim> VecDim;
    typedef boost::array<unsigned,numNodesPerElement>         IndexSimplex;
    typedef boost::array<VecDim,numNodesPerElement>           Simplex;

    typedef typename std::vector<VecDim>::iterator       CoordinateIterator;
    typedef typename std::vector<IndexSimplex>::iterator IndexSimplexIterator;

public:
    /** Constructor via SMF file */
    Surface(std::istream & smf)
    {
        //! Read and validate the header of the file
        corlib::SmfHead smfHead;
        smfHead.read(smf);

        //! number of nodes and elements
        unsigned numNodes, numElements;
        smf >> numNodes >> numElements;
        smf.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

        //! Read in the coordinates
        coordinates_.resize(numNodes);
        for (unsigned n = 0; n < numNodes; n ++) {
            // read x
            VecDim x;
            for (unsigned d = 0; d < dim; d ++) smf >> x[d];
            smf.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
            coordinates_[n] = x;
        }

        //! Read in the connectivity
        connectivity_.resize(numElements);
        for (unsigned e = 0; e < numElements; e ++) {
            // read simplex
            IndexSimplex s;
            for (unsigned v = 0; v < numNodesPerElement; v ++) smf >> s[v];
            smf.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
            connectivity_[e] = s;
        }
    }

    CoordinateIterator coordinatesBegin() { return coordinates_.begin(); }
    CoordinateIterator coordinatesEnd()   { return coordinates_.end(); }
    IndexSimplexIterator connectivityBegin() { return connectivity_.begin(); }
    IndexSimplexIterator connectivityEnd()   { return connectivity_.end(); }

private:
    std::vector<VecDim>       coordinates_;
    std::vector<IndexSimplex> connectivity_;
};

//------------------------------------------------------------------------------
template <typename POINT>
void lset::Grid<POINT>::readHeader_(std::istream &sgf)
{
    // read header
    corlib::SmfHead smfHead;
    smfHead.read(sgf);

    // read extents (ix_min, ix_max, iy_min, iy_max, iz_min, iz_max)
    for (unsigned e = 0; e < 6; e ++) sgf >> extents_[e];
    sgf.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

    // numbers of cells per direction (pad with zeros)
    numCells_ = detail_::makeMultiIndex_0<dim>(extents_[1] - extents_[0],
                                               extents_[3] - extents_[2],
                                               extents_[5] - extents_[4]);

    std::cout << numCells_[0] << " " << numCells_[1] << " " << numCells_[2] << std::endl;

    FTL_VERIFY(extents_[0] == 0 && extents_[2] == 0 && extents_[4] == 0);

    // check dimensions
    if (dim < 3) FTL_VERIFY((extents_[5] - extents_[4]) == 0);
    if (dim < 2) FTL_VERIFY((extents_[3] - extents_[2]) == 0);
    FTL_VERIFY((extents_[1] - extents_[0]));
}

//------------------------------------------------------------------------------
template <typename POINT>
void lset::Grid<POINT>::readPoints_(std::istream &sgf)
{
    // loop over points
    for (int ix = 0; ix <= numCells_[0]; ++ix) {
        for (int iy = 0; iy <= numCells_[1]; ++iy) {
            for (int iz = 0; iz <= numCells_[2]; ++iz) {
                // get index and coordinates
                VecIntDim idx = detail_::makeMultiIndex_0<dim>(ix, iy, iz);
                Vec3 coords;
                for (unsigned d = 0; d < 3; ++d) sgf >> coords[d];

                // store point
                Point *p = new Point(idx, detail_::makeMultiIndex<dim>(coords));
                points_.push_back(p);
            }
        }
    }
}

//------------------------------------------------------------------------------
template <typename POINT>
void lset::Grid<POINT>::writeCoords_(std::ostream & vts)
{
    vts << "<?xml version=\"1.0\"?>" << std::endl
        << "<VTKFile type=\"StructuredGrid\" version="
        << "\"0.1\" byte_order=\"LittleEndian\">" << std::endl
        << "  <StructuredGrid WholeExtent=\"";

    std::copy(extents_.begin(), extents_.end(),
              std::ostream_iterator<unsigned>(vts, " "));

    vts << "\">" << std::endl
        << "    <Piece Extent=\"";

    std::copy(extents_.begin(), extents_.end(),
              std::ostream_iterator<unsigned>(vts, " "));

    vts << "\">" << std::endl;

    //! Write nodal coordinates
    vts << "      <Points>" << std::endl;
    this -> writePointQuantity_(vts, corlib::accessorFun(&Point::getCoordinates,
                                                         "Coordinates"), true);
    vts << "      </Points>" << std::endl;
}

//------------------------------------------------------------------------------
template<typename POINT>
template<typename PACC>
void lset::Grid<POINT>::writePointQuantity_(std::ostream & vts,
                                            PACC pointAcc,
                                            const bool coords,
                                            const bool noXML)
{
    // get information about quantity
    const std::string quantityName = pointAcc.name();
    const unsigned numComponents = PACC::numComponents;
    typedef typename PACC::ValueType ValueType;
    const std::string nameOfType = corlib::detail_::VtkDataType<ValueType>::name();

    //-- Dirty hack follows here, because Paraview does not understand 2D --//
    unsigned numVTUcomp = (numComponents == 2 ? 3 : numComponents);
    if(coords) numVTUcomp = 3;

    if (not noXML) {
        vts << "        <DataArray type=\"" << nameOfType
            << "\" NumberOfComponents=\"" << numVTUcomp
            << "\" Name=\""               << quantityName
            << "\" format=\""             << "ascii"
            << "\">"
            << std :: endl;
    }

    // go through all points
    for (unsigned n = 0; n < points_.size(); ++n) {
        // obtain nodal quantity
        typename PACC::ReturnType quantity = pointAcc(points_[n]);
        //----------------------------------------------------------------------
        // write the quantity
        if (coords or (numComponents == 2)) {
            // write a vec3 and fill with zeros
            ublas::bounded_vector<double,3> q;
            corlib::detail_::padWithZeros(q, quantity);
            corlib::detail_::writeQuantity(vts, q);
        }
        else corlib::detail_::writeQuantity(vts, quantity);

        vts << std::endl;
    }
    vts << "      </DataArray>" << std :: endl;
}

//------------------------------------------------------------------------------
template <typename POINT>
void lset::Grid<POINT>::finalize(std::ostream & vts)
{
    vts << "    </Piece>" << std::endl
        << "  </StructuredGrid>" << std::endl
        << "</VTKFile>" << std::endl;
    return;
}


#endif
