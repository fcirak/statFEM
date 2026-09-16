// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//
//                 Computational Structural Mechanics Lab
//                         University of Cambridge
//
// This software is copyrighted and all rights are retained by the University of
// Cambridge. All use, disclosure, and/or reproduction of any part not expressly
// authorised by the University of Cambridge is prohibited. (C) 2010
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

//! @file ShapeToVtk.hpp

#ifndef corlib_shapetovtk_h
#define corlib_shapetovtk_h

//------------------------------------------------------------------------------
// headers
#include <string>
#include <cassert>

#include <corlib/Shape.hpp>


namespace corlib {

    //--------------------------------------------------------------------------
    //! gives the vtkType of the combination (SHAPE,numPoints) 
    template< shape SHAPE, unsigned numPoints > struct VtkType;
    //! \cond SKIPDOX
    template<> struct VtkType< POINT,         1 > { static const unsigned vtk =  1; };
    template<> struct VtkType< LINE,          2 > { static const unsigned vtk =  3; };
    template<> struct VtkType< LINE,          3 > { static const unsigned vtk = 21; };
    template<> struct VtkType< TRIANGLE,      3 > { static const unsigned vtk =  5; };
    template<> struct VtkType< TRIANGLE,      6 > { static const unsigned vtk = 22; };
    template<> struct VtkType< QUADRILATERAL, 4 > { static const unsigned vtk =  9; };
    template<> struct VtkType< QUADRILATERAL, 8 > { static const unsigned vtk = 23; };
    template<> struct VtkType< QUADRILATERAL, 9 > { static const unsigned vtk = 23; };
    template<> struct VtkType< TETRAHEDRON,   4 > { static const unsigned vtk = 10; };
    template<> struct VtkType< TETRAHEDRON,  10 > { static const unsigned vtk = 24; };
    template<> struct VtkType< HEXAHEDRON,    8 > { static const unsigned vtk = 12; };
    template<> struct VtkType< HEXAHEDRON,   20 > { static const unsigned vtk = 25; };
    template<> struct VtkType< HEXAHEDRON,   27 > { static const unsigned vtk = 25; };
    template<> struct VtkType< PIXEL,         4 > { static const unsigned vtk =  8; };
    template<> struct VtkType< VOXEL,         8 > { static const unsigned vtk = 11; };

    // work around for multi-noded pixel / voxel elements
    template<unsigned numPoints> struct VtkType< PIXEL, numPoints> { static const unsigned vtk = 0; };
    template<unsigned numPoints> struct VtkType< VOXEL, numPoints> { static const unsigned vtk = 0; };
    //! \endcond

    //--------------------------------------------------------------------------
    //! gives the Shape based on vtkType
    template< unsigned VTKTYPE > struct VtkTypeShape;
    //! \cond SKIPDOX
    template<> struct VtkTypeShape<  1 > { static const enum shape shape =         POINT; };
    template<> struct VtkTypeShape<  3 > { static const enum shape shape =          LINE; };
    template<> struct VtkTypeShape< 21 > { static const enum shape shape =          LINE; };
    template<> struct VtkTypeShape<  5 > { static const enum shape shape =      TRIANGLE; };
    template<> struct VtkTypeShape< 22 > { static const enum shape shape =      TRIANGLE; };
    template<> struct VtkTypeShape<  9 > { static const enum shape shape = QUADRILATERAL; };
    template<> struct VtkTypeShape< 23 > { static const enum shape shape = QUADRILATERAL; };
    template<> struct VtkTypeShape< 10 > { static const enum shape shape =   TETRAHEDRON; };
    template<> struct VtkTypeShape< 24 > { static const enum shape shape =   TETRAHEDRON; };
    template<> struct VtkTypeShape< 12 > { static const enum shape shape =    HEXAHEDRON; };
    template<> struct VtkTypeShape< 25 > { static const enum shape shape =    HEXAHEDRON; };
    template<> struct VtkTypeShape<  8 > { static const enum shape shape =         PIXEL; };
    template<> struct VtkTypeShape< 11 > { static const enum shape shape =         VOXEL; }; 
    //! \endcond

    //--------------------------------------------------------------------------
    //! 
    template< unsigned VTKTYPE > struct VtkTypeNodeNum;
    //! \cond SKIPDOX

    // line element with two nodes
    template<> struct VtkTypeNodeNum<3>  { static const unsigned np = 2; };

    // line element with three nodes
    template<> struct VtkTypeNodeNum<21> { static const unsigned np = 3; };

    // triangle with three nodes
    template<> struct VtkTypeNodeNum<5>  { static const unsigned np = 3; };

    // triangle with six nodes
    template<> struct VtkTypeNodeNum<22> { static const unsigned np = 6; };

    // quadrilateral with four nodes
    template<> struct VtkTypeNodeNum<9>  { static const unsigned np = 4; };

    // quadrilateral with eight nodes
    template<> struct VtkTypeNodeNum<23> { static const unsigned np = 8; };

    // tetrahedron with four nodes
    template<> struct VtkTypeNodeNum<10> { static const unsigned np = 4; };

    // tetrahedron with ten nodes
    template<> struct VtkTypeNodeNum<24> { static const unsigned np = 10;};

    // hexahedron with eight nodes
    template<> struct VtkTypeNodeNum<12> { static const unsigned np =  8;};

    // hexahedron with twenty nodes
    template<> struct VtkTypeNodeNum<25> { static const unsigned np = 20;};

    // pixel with four nodes
    template<> struct VtkTypeNodeNum<8>  { static const unsigned np =  4; };

    // voxel with eight nodes
    template<> struct VtkTypeNodeNum<11> { static const unsigned np =  8; };

    // work around for multi-noded pixel / voxel elements
    template<> struct VtkTypeNodeNum<0> { static const unsigned np = 0; };

    //! \endcond

    //--------------------------------------------------------------------------
    //! Determine VtkType of given element shape \p eleShape with \p numPoints 
    //! points attached
    //!
    //! \param[in]   eleShape   element shape
    //! \param[in]   numPoints  number of element points
    //! \return                 VtkType number (default is 0)
    unsigned vtkTypeOfShape( const enum shape & eleShape,
                             const unsigned & numPoints )
    {
        // initialise output
        unsigned vtkTypeNum = 0;

        switch ( eleShape ) {
        case UNDEFINED : 
            FTL_VERIFY_DESCRIPTIVE( false, "Cannot convert undefined shape" );
            break;
        case LINE :
            switch ( numPoints ) {
            case 2  : vtkTypeNum = VtkType< LINE, 2 >::vtk;           break;
            case 3  : vtkTypeNum = VtkType< LINE, 3 >::vtk;           break;
            default : FTL_VERIFY_DESCRIPTIVE( false, "Cannot find VTK type of shape with given number of points" );
            }
            break;
        case TRIANGLE :
            switch ( numPoints ) {
            case 3  : vtkTypeNum = VtkType< TRIANGLE, 3 >::vtk;       break;
            case 6  : vtkTypeNum = VtkType< TRIANGLE, 6 >::vtk;       break;
            default : FTL_VERIFY_DESCRIPTIVE( false, "Cannot find VTK type of shape with given number of points" );
            }
            break;
        case QUADRILATERAL :
            switch ( numPoints ) {
            case 4  : vtkTypeNum = VtkType< QUADRILATERAL, 4 >::vtk;  break;
            case 8  : vtkTypeNum = VtkType< QUADRILATERAL, 8 >::vtk;  break;
            case 9  : vtkTypeNum = VtkType< QUADRILATERAL, 9 >::vtk;  break;
            default : FTL_VERIFY_DESCRIPTIVE( false, "Cannot find VTK type of shape with given number of points" );
            }
            break;
        case TETRAHEDRON :
            switch ( numPoints ) {
            case 4  : vtkTypeNum = VtkType< TETRAHEDRON, 4 >::vtk;    break;
            case 10 : vtkTypeNum = VtkType< TETRAHEDRON, 10 >::vtk;   break;
            default : FTL_VERIFY_DESCRIPTIVE( false, "Cannot find VTK type of shape with given number of points" );
            }
            break;
        case HEXAHEDRON :
            switch ( numPoints ) {
            case 8  : vtkTypeNum = VtkType< HEXAHEDRON, 8 >::vtk;     break;
            case 20 : vtkTypeNum = VtkType< HEXAHEDRON, 20 >::vtk;    break;
            case 27 : vtkTypeNum = VtkType< HEXAHEDRON, 27 >::vtk;    break;
            default : FTL_VERIFY_DESCRIPTIVE( false, "Cannot find VTK type of shape with given number of points" );
            }
            break;
        case PIXEL :
            switch ( numPoints ) {
            case 4  : vtkTypeNum = VtkType< PIXEL, 4 >::vtk;          break;
            default : FTL_VERIFY_DESCRIPTIVE( false, "Cannot find VTK type of shape with given number of points" );
            }
            break;
        case VOXEL :
            switch ( numPoints ) {
            case 8  : vtkTypeNum = VtkType< VOXEL, 8 >::vtk;          break;
            default : FTL_VERIFY_DESCRIPTIVE( false, "Cannot find VTK type of shape with given number of points" );
            }
            break;
        default : 
            FTL_VERIFY_DESCRIPTIVE( false, "Cannot find VTK type of shape" );
        }

        return vtkTypeNum;
    }
    
    //--------------------------------------------------------------------------
    //! Determine number of nodes of given VtkType 
    //!
    //! \param[in]   vtkType    VtkType number
    //! \return                 number of element points (default is 0)
    unsigned vtkTypeNumberOfNodes( const unsigned vtkType )
    {
        // initialise output
        unsigned numNodes = 0;
	
        switch ( vtkType ) {
        case  3 : numNodes = VtkTypeNodeNum<  3 >::np; break;
        case 21 : numNodes = VtkTypeNodeNum< 21 >::np; break;
        case  5 : numNodes = VtkTypeNodeNum<  5 >::np; break;
        case 22 : numNodes = VtkTypeNodeNum< 22 >::np; break;
        case  9 : numNodes = VtkTypeNodeNum<  9 >::np; break;
        case 23 : numNodes = VtkTypeNodeNum< 23 >::np; break;
        case 10 : numNodes = VtkTypeNodeNum< 10 >::np; break;
        case 24 : numNodes = VtkTypeNodeNum< 24 >::np; break;
        case 12 : numNodes = VtkTypeNodeNum< 12 >::np; break;
        case 25 : numNodes = VtkTypeNodeNum< 25 >::np; break;
        default : FTL_VERIFY_DESCRIPTIVE( false, "VTK type not known." );
        }

        return numNodes;
    }

}


#endif
