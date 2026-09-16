// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//
//                 Computational Structural Mechanics Lab
//                         University of Cambridge
//
// This software is copyrighted and all rights are retained by the University of
// Cambridge. All use, disclosure, and/or reproduction of any part not expressly
// authorised by the University of Cambridge is prohibited. (C) 2010
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

//! @file MeshSizeInfo.hpp

#ifndef corlib_meshsizeinfo_h
#define corlib_meshsizeinfo_h
//------------------------------------------------------------------------------
#include <vector>
#include <algorithm>
#include <numeric>

#include <Eigen/Core>

#include <corlib/verify.hpp>
#include <corlib/Shape.hpp>

//------------------------------------------------------------------------------
namespace corlib {
    template<typename ELEMENT> class MeshSizeInfo;
}

//------------------------------------------------------------------------------
/** \brief Functor for mesh size statistics
 *  \details Collects the mesh size parameters of all elements and allows to
 *  query the average or write some statistics about the mesh.
 *  \tparam ELEMENT Type of element
 */
template<typename ELEMENT>
class corlib::MeshSizeInfo
    : public boost::function<void( const ELEMENT* ) >
{
public:
    static const corlib::shape elemShape   = ELEMENT::myShape;
    static const unsigned      numVertices = corlib::ShapeTraits<elemShape>::numVertices;

    //! Empty constructor
    MeshSizeInfo( ) { };
    //! Destructor to clear the data
    ~MeshSizeInfo( ) { meshSizes_.clear( ); }

    /** \brief Operation is to store the mesh size parameters in a container
     *  \param ep  Pointer to the element it is operating on  */
    void operator() ( const ELEMENT* ep );

    //! Write some statistice to a stream
    std::ostream & writeStats( std::ostream & out ) const;

    //! Return the maximum element size
    double giveMaximum() const;

    //! Return the minimum element size
    double giveMinimum() const;

    //! Return the average mesh size
    double giveAverage() const;

private:
    std::vector<double> meshSizes_; //!< Container of mesh sizes
};

//------------------------------------------------------------------------------
/** Function call operator which stores the element's size indicator. Note
 *  that this is based on the element's corner vertices only.
 */
template<typename ELEMENT>
void corlib::MeshSizeInfo<ELEMENT>::operator()( const ELEMENT * ep )
{
    // get corner nodes as shape's vertices
    std::array<Eigen::VectorXd, numVertices> verts;

    typename ELEMENT::NodeConstIterator itnb = ep -> nodesBegin();
    typename ELEMENT::NodeConstIterator itne = ep -> nodesEnd();

    for ( unsigned v = 0; (v < numVertices) and ( itnb != itne ); v ++, itnb ++ ) 
    {
        verts[v] = (*itnb) -> giveCoordinates();        
    }

    // store mesh size indicator
    meshSizes_.push_back( corlib::elementSize<elemShape>( verts ) );
}

//------------------------------------------------------------------------------
/** Write statistics like maximum, minimum, and average of the stored mesh sizes
 *  to a given stream
 *  \param[in,out] out  Output stream
 *  \retval        out  Output stream
 */
template<typename ELEMENT>
std::ostream & corlib::MeshSizeInfo<ELEMENT>::writeStats( std::ostream & out ) const
{
    out << "--------------------------------------------------" << std::endl
        << " Mesh size info: " << std :: endl
        << " Minimal element size: "
        << this -> giveMinimum( )
        << std::endl
        << " Maximal element size: "
        << this -> giveMaximum( )
        << std::endl
        << " Average element size: "
        << this -> giveAverage( )
        << std :: endl
        << "--------------------------------------------------" << std::endl;

    return out;
}

//------------------------------------------------------------------------------
/** Return the maximum mesh size
 *  \retval h_max  The maximum mesh size
 */
template<typename ELEMENT>
double corlib::MeshSizeInfo<ELEMENT>::giveMaximum() const
{
    FTL_VERIFY( meshSizes_.begin( ) != meshSizes_.end( ) );

    return *std::max_element( meshSizes_.begin( ), meshSizes_.end( ) );
}

//------------------------------------------------------------------------------
/** Return the minimum mesh size
 *  \retval h_min  The minimum mesh size
 */
template<typename ELEMENT>
double corlib::MeshSizeInfo<ELEMENT>::giveMinimum() const
{
    FTL_VERIFY( meshSizes_.begin( ) != meshSizes_.end( ) );

    return *std::min_element( meshSizes_.begin( ), meshSizes_.end( ) );
}

//------------------------------------------------------------------------------
/** Compute the average mesh size: \f$ \bar h = 1/N \sum_{i=1}^N h_i \f$,
 *  using std::accumulate.
 *  \retval h  The average mesh size
 */
template<typename ELEMENT>
double corlib::MeshSizeInfo<ELEMENT>::giveAverage() const
{
    FTL_VERIFY( meshSizes_.begin( ) != meshSizes_.end( ) );

    return std::accumulate( meshSizes_.begin( ),
                            meshSizes_.end( ), 0. ) / meshSizes_.size( );
}


#endif
