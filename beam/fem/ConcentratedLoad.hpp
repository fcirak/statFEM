// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//
//                 Computational Structural Mechanics Lab
//                         University of Cambridge
//
// This software is copyrighted and all rights are retained by the University of
// Cambridge. All use, disclosure, and/or reproduction of any part not expressly
// authorised by the University of Cambridge is prohibited. (C) 2010
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

//! @file ConcentratedLoad.hpp

#ifndef beam_fem_concentratedload_h
#define beam_fem_concentratedload_h

//------------------------------------------------------------------------------
#include <iostream>
#include <tuple>

#include <corlib/eigenX.hpp>
#include <corlib/verify.hpp>

//------------------------------------------------------------------------------
namespace beam {
    namespace fem {

        template< typename ELEMENT >
        class ConcentratedLoad;

        namespace eigenX = corlib::eigenX;

    }
}

//------------------------------------------------------------------------------
/** \brief Management of nodally concentrated loads (forces and torques)
 *
 *  \details This object reads the nodally given concentrated loads
 *  described in a file. The concentrated loads are forces and torques/moments.
 *  The found concentrated loads are attached to to an element of which
 *  the loaded node is a vertex. The loaded elements togethe with load
 *  specifications are stored in #loads_.
 *  The object provides an operator() which can be parsed by an iteration
 *  of all mesh elements.
 *
 *  \tparam  ELEMENT  Type of element
 */
template< typename ELEMENT >
class beam::fem::ConcentratedLoad
{
public:
    typedef ELEMENT                                      ElementType;
    typedef typename ELEMENT::Node                       NodeType;

private:
    static const unsigned localDim_ = ElementType::localDim;

private:
    typedef ElementType*                                 ElementPtr_;
    typedef typename ElementType::VecLDim                VecLDim_;
    typedef eigenX::VectorSd<3>                          Vec3_;
    typedef std::tuple< VecLDim_, Vec3_, Vec3_ >         Load_;
    typedef std::map< ElementType*, Load_ >              LoadMap_;
    typedef typename LoadMap_::iterator                  LoadMapIter_;
    typedef typename LoadMap_::const_iterator            LoadMapConstIter_;


public:
    //! Constructor given the function object representing the body force
    //!
    //! The input file is parsed. Nodal concentrated loads are found.
    //! The concentrated loads are attached to an element of which
    //! the loaded node is a vertex.
    //!
    //! \tparam          MESH   The finite element mesh type
    //! \param[in,out]   inp    input stream of file
    //! \param[in]       mesh   Finite element mesh
    template< typename MESH >
    ConcentratedLoad( std::istream& inp,
                      const MESH& mesh );

    //! Skip comment lines
    //!
    //! \param[in,out]   inp    input stream of file
    void skipCommentLines( std::istream& inp ) const;

    //! Set (global) load factor
    void setFactor( const double& factor )
    {
        factor_ = factor;
        return;
    }

    //! Basic operation to add concentrated loads
    void operator()( ELEMENT * ep );

    //! Print concentrated loads
    std::ostream& print( std::ostream& os ) const;

private:
    //! character indicating a comment line in input file
    const char    commentChar_;

    //! concentrated loads
    LoadMap_      loads_;

    //! Multiplier for the forces
    double        factor_;

};


//------------------------------------------------------------------------------
#include "ConcentratedLoad.ipp"


//------------------------------------------------------------------------------
#endif
