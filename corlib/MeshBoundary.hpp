// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//
//                 Computational Structural Mechanics Lab
//                         University of Cambridge
//
// This software is copyrighted and all rights are retained by the University of
// Cambridge. All use, disclosure, and/or reproduction of any part not expressly
// authorised by the University of Cambridge is prohibited. (C) 2010
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

//! @file MeshBoundary.hpp

#ifndef corlib_meshboundary_h
#define corlib_meshboundary_h

//------------------------------------------------------------------------------
#include <corlib/Shape.hpp>
#include <corlib/misc.hpp>

#include <boost/static_assert.hpp>
#include <boost/type_traits/is_same.hpp>

#include <vector>
#include <iostream>
#include <algorithm>

namespace corlib{
    template<typename SELEMENT> class MeshBoundary;
    
}

//------------------------------------------------------------------------------
/** \brief Boundary mesh (or a part of it) 
 *  \details This object stores the boundary elements of a mesh. 
 *  The stored boundary elements can be a subset of the entire boundary mesh. 
 *  \tparam SELEMENT  Type of the stored boundary element
 */
template<typename SELEMENT>
class corlib::MeshBoundary
{
public:
    typedef SELEMENT SurfElement;

    //! Constructor with access to the mesh and input stream
    template<typename MESH>
    MeshBoundary( std::istream & inp, MESH & mesh ) 
    {

        // some sanity checks for the right types
        {
            typedef typename MESH::Element DomainElement;

            FTL_STATIC_ASSERT_MSG( (boost::is_same<typename SurfElement::Node,
                                                   typename MESH::Node>::value),
                                   "Node types do not mathc");

            FTL_STATIC_ASSERT_MSG( (SurfElement::myShape == 
                                    corlib::ShapeTraits<DomainElement::myShape>::faceShape ),
                                   "Shapes do not match" );
            
            FTL_STATIC_ASSERT_MSG( (SurfElement::ShapeFun::sfc ==
                                      DomainElement::ShapeFun::sfc), 
                                     "Shape function types do not match" );

            FTL_STATIC_ASSERT_MSG( (SurfElement::ShapeFun::degree ==
                                    DomainElement::ShapeFun::degree),
                                   "Shape function degrees do not match" );
            
        }

        // get number of faces, i.e., surface elements
        unsigned numFaces; 
        inp >> numFaces;
        inp.ignore( std::numeric_limits<std::streamsize>::max(), '\n' );
        
        // get nodes from mesh
        std::vector< typename MESH::Node * > nodes( mesh.numNodes() );
        std::copy( mesh.nodesBegin(), mesh.nodesEnd(), nodes.begin() );

        // read surface element connectivity
        surfaceElements_.reserve( numFaces );

        for ( unsigned f = 0; f < numFaces; f ++ ) {
            SurfElement * surfElement = new SurfElement;
            surfElement -> readSelf( inp, nodes );
            inp.ignore( std::numeric_limits<std::streamsize>::max(), '\n' );
            surfaceElements_.push_back( surfElement );
        }

    }    

    //! delete surface elements
    ~MeshBoundary() {
        this -> iterateOverElements( corlib::deleteFunctor( ) );
    }
        
    
    //! Apply op to all elements
    template< typename OP >
    OP iterateOverElements( OP op ) 
    {
        typedef typename std::vector<SurfElement*>::iterator ElemIter;
        ElemIter begin = surfaceElements_.begin( );
        ElemIter end   = surfaceElements_.end( );
        return std::for_each( begin, end, op );
    }


private:
    //! Elements which form a surface mesh
    std::vector< SurfElement *> surfaceElements_;
};


#endif

