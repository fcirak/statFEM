//
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//
//                 Computational Structural Mechanics Lab
//                         University of Cambridge
//
// This software is copyrighted and all rights are retained by the CSMLab. All
// use, disclosure, and/or reproduction of any part not expressly authorized by
// F Cirak is prohibited. (C) 2012.
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//

//! @file   partitioners/utils/misc.hpp
//! @author Burkhard Bornemann, Matija Kecman
//! @date   01/2012

#ifndef tools_input_partitioners_utils_misc_h
#define tools_input_partitioners_utils_misc_h

#include <boost/numeric/ublas/vector.hpp>

#include <vector>
#include <algorithm>

//==============================================================================
namespace tools{
    namespace input{
        namespace partitioners{
            namespace utils{

                namespace ublas = boost::numeric::ublas;

                // nested for_each loop for iterating on a container of containers
                template< typename OUTERINPITER, typename INNERFUN >
                INNERFUN for_each_2d( OUTERINPITER outerFirst, OUTERINPITER outerLast,
                                      INNERFUN innerFun )
                {
                    for ( ; outerFirst != outerLast; ++outerFirst )
                        innerFun = std::for_each( outerFirst->begin(), outerFirst->end(), innerFun );
                    return innerFun;
                }

                //--------------------------------------------------------------
                /** Computes a vector of the coordinates of the element centers
                 *  of the mesh described by the given coordinates and connectivity.
                 *
                 *  \param[in]  coords  Coordinates of the given mesh
                 *  \param[in]  connec  Connectivity of the given mesh
                 *  \param[out] centers Centers of the elements in the mesh
                 */
                void computeElementCenters( const std::vector<std::vector<unsigned> >           & connec, 
                                            const std::vector<ublas::bounded_vector<double,3> > & coords,
                                            std::vector<ublas::bounded_vector<double,3> >       & centers )
                {
                    typedef std::vector<std::vector<unsigned> >::size_type VecVecU_sz;
                    typedef ublas::bounded_vector<double,3>::size_type VecDim3_sz;

                    const VecVecU_sz numElements        = connec.size();
                    const VecDim3_sz numNodesPerElement = connec[0].size();

                    // compute the centres of the elements
                    for ( VecVecU_sz e = 0; e < numElements; ++e ) {

                        ublas::bounded_vector<double,3> center;
                        center.clear();

                        for( VecDim3_sz i = 0; i < numNodesPerElement; ++i ) {
                            const unsigned nodeId = connec[e][i];
                            center += coords[nodeId];
                        }

                        center /= static_cast<double>( numNodesPerElement );
                        centers.push_back( center );
                    }

                    return;
                }

            }
        }
    }
}

#endif
