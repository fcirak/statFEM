// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//
//                 Computational Structural Mechanics Lab
//                         University of Cambridge
//
// This software is copyrighted and all rights are retained by the University of
// Cambridge. All use, disclosure, and/or reproduction of any part not expressly
// authorised by the University of Cambridge is prohibited. (C) 2010
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

//! @file   HighResVTU.hpp

#ifndef beam_fem_highresvtu_h
#define beam_fem_highresvtu_h
//------------------------------------------------------------------------------
//! std    includes
#include <ostream>
#include <algorithm>
#include <iterator>
#include <string>
//! Eigen includes
#include <Eigen/Core>
//! corlib includes
#include <corlib/eigenX.hpp>
#include <corlib/EvaluateField.hpp>
#include <corlib/Accessor.hpp>
#include <corlib/VTUwriter.hpp> // to access the free write functions

//------------------------------------------------------------------------------
namespace beam{
    namespace fem{
        template<typename MESH>  class HighResVTU;

        namespace eigenX = corlib::eigenX;

        namespace detail_{

            //------------------------------------------------------------------
            /** Workaround for element accessors which do not depend on the 
             *  coordinate.
             *  \tparam EACC  Accessor to element quantity
             *  \tparam ARITY Function arity of the accessor
             */
            template<typename EACC, unsigned ARITY = EACC::arity>
            struct ElementAccessPolicy;

            //! \cond SKIPDOX
            template<typename EACC> // no dependence on local coordinate
            struct ElementAccessPolicy<EACC,1> 
            {
                // ignore second argument
                template<typename DUMMY>
                typename EACC::result_type operator()( EACC & cacc,
                                                       typename EACC::argument_type c,
                                                       DUMMY & dummy ) const
                {
                    return cacc( c );
                }
            };

            template<typename EACC> // depends on local coordinate
            struct ElementAccessPolicy<EACC,2>
            {
                // use both arguments
                typename EACC::result_type operator()( EACC& cacc, 
                                                       typename EACC::arg1_type c,
                                                       typename EACC::arg2_type s ) const
                {
                    return cacc( c, s );
                }
            };
            //! \endcond
        }
    }
}

//------------------------------------------------------------------------------
/** \brief Writer for unstructured meshes in VTU-format with high resolution
 *  
 *  \details The output data is completely obtained by evaluation of the 
 *  desired quantities (including geometry) at local coordinates of each cell.
 *  In the standard setup this is done at the element's origin and for 
 *  completeness at the right-extremal points for the last points in the mesh.
 *  The element data are obtained by directly accessing the elements.
 *
 *  This object has the functionality to sample the output at a higher rate.
 *  This means that one can write out a mesh which is finer than the computational
 *  mesh. This feature is useful for mesh refinement or better resolution in case
 *  of higher-order splines.
 *  \tparam MESH  Type of mesh to post-process
 */
template<typename MESH>  
class beam::fem::HighResVTU
{
public:
    typedef MESH                                           Mesh;
    static const unsigned dim = Mesh::Element::localDim;

    //! Cstor with stream and mesh reference
    HighResVTU( Mesh * meshPtr, 
                std::ostream & vtu, 
                const int resolution = 1 ) 
        : meshPtr_(     meshPtr ),
          vtu_(         vtu ),
          resolution_(  resolution )
    { }

    //--------------------------------------------------------------------------
    //! Write the XML-header and extents array
    void writeMesh()
    {
        const unsigned nElements = resolution_ * (meshPtr_ -> numElements());
        const bool     isClosed  = meshPtr_ -> hasClosedTopology();
        const unsigned nNodes    = nElements + (isClosed ? 0 : 1);

        // write header of VTU file
        vtu_ << "<?xml version=\"1.0\"?>" << std :: endl
             << "<VTKFile type=\"UnstructuredGrid\" byte_order=\"LittleEndian\">" << std :: endl
             << "  <UnstructuredGrid>" << std :: endl
             << "  <Piece NumberOfPoints=\"" << nNodes 
             << "\"  NumberOfCells=\"" << nElements 
             << "\">" << std :: endl;
        this -> writePoints();

        this -> writeCells( isClosed );
        
        return;
    }

    //--------------------------------------------------------------------------
    //! Write footer
    void finalize()
    {
        vtu_ << "    </Piece>" << std::endl
             << "  </UnstructuredGrid>" << std::endl
             << "</VTKFile>" << std::endl;
        return;
    }

private:
    //--------------------------------------------------------------------------
    //! Write nodal coordinates
    void writePoints()
    {
        vtu_ << "      <Points>" << std::endl;
        this -> writePointQuantity( corlib::accessorFun( &Mesh::Node::giveCoordinates, 
                                                         "Coordinates"), true );
        vtu_ << "      </Points>" << std::endl;
        return;
    }

    //---------------------------------------------------------------------------
    //! Write connectivity, hard-wired to sequential line elements
    void writeCells( const bool isClosed ) 
    {
        const unsigned nElements = resolution_ * (meshPtr_ -> numElements());

        vtu_ << "      <Cells>\n";
        // connectivity
        {
            vtu_ << "        <DataArray type=\"Int32\" Name=\"connectivity\" "
                 << "NumberOfComponents=\"1\" format=\"ascii\">\n";
            unsigned n = 0;
            for ( ; n < nElements - 1; n ++ ) {
                vtu_ << n << " " << n+1 << "\n";
            }
            vtu_ << n << " " << (isClosed ? 0 : n + 1) << "\n" 
                 << "        </DataArray>\n";
        }

        // offsets
        {
            vtu_ << "        <DataArray type=\"Int32\" Name=\"offsets\" "
                 << "NumberOfComponents=\"1\" format=\"ascii\">\n";
            for ( unsigned n  = 0; n < nElements; n ++ )
                vtu_ << 2 * (n+1) << "\n";
            vtu_ << "        </DataArray>\n";
        }
        
        // vtk element type
        {
            vtu_ << "        <DataArray type=\"Int32\" Name=\"types\" "
                 << "NumberOfComponents=\"1\" format=\"ascii\">\n";
            for ( unsigned n = 0; n < nElements; n ++ ) 
                vtu_ << 3 << "\n";
            vtu_ << "        </DataArray>\n";
        }

        vtu_ << "      </Cells>\n";
    }

public:
    //--------------------------------------------------------------------------
    //! this function writes the final tags to close the file
    void openPointData( )
    {
        vtu_ << "      <PointData>" << std :: endl;
        return;
    }

    //--------------------------------------------------------------------------
    //! this function writes the final tags to close the file
    void closePointData( )
    {
        vtu_ << "      </PointData>" << std :: endl;
        return;
    }

    //--------------------------------------------------------------------------
    /** \brief Computation of point data by evaluation of the field in elements
     *  \details For any field 'u' (including the geometry), the approximation
     *  \f[
     *         u_h(\xi) = \sum B_k(\xi) u_k
     *  \f]
     *  is assumed within every element. This approximation is here evaluated
     *  for every element at sub-points as desired.
     *  \tparam    PACC      Type of point accessor
     *  \param[in] pointAcc  Accessor for the point quantity
     *  \param[in] coords    Flag for writing coordinates
     */
    template<typename PACC>
    void writePointQuantity( PACC pointAcc, 
                             const bool coords = false )
    {
        // get name and amount of quantities 
        const std::string quantityName  = pointAcc.name();
        const unsigned    numComponents = PACC::numComponents;

        typedef typename PACC::ValueType ValueType;
        const std::string nameOfType = corlib::detail_::VtkDataType<ValueType>::name();

        //-- Dirty hack follows here, because Paraview does not understand 2D --//
        unsigned numVTUcomp = ( numComponents == 2 ? 3 : numComponents );
        if( coords ) numVTUcomp = 3;

        vtu_ << "        <DataArray type=\"" << nameOfType 
             << "\" NumberOfComponents=\"" << numVTUcomp 
             << "\" Name=\""               << quantityName 
             << "\" format=\""             << "ascii"
             << "\">"
             << std :: endl;

        const unsigned nElementsCoarse = meshPtr_ -> numElements();
        
        // all cells in z-direction
        for ( unsigned c = 0; c <= nElementsCoarse; c++ ) {
            
            // all sub-cells in z-direction
            const int numSubElements = ( c == nElementsCoarse ? 1 : resolution_ );
            for ( unsigned f = 0; f < numSubElements; f++ ) {
                
                // check index limits to get the right-most cell-pointer
                unsigned i = c;
                if ( c == nElementsCoarse ) i = c - 1;
                    
                // evaluation point in sub-cell mesh
                typename Mesh::Element::VecLDim xi; 
                xi[0] = static_cast<double>( f ) / static_cast<double>( resolution_ );
                                    
                // in the right-most limit
                if ( c == nElementsCoarse ) xi[0] = 1.;
                
                // get datum evaluated at xi
                corlib::EvaluateField<PACC,typename Mesh::Element> ef( pointAcc );
                const typename PACC::result_type quantity = 
                    ef( (meshPtr_ -> getElementPointer( i )), xi );
                    
                // work-around for coordinates and 2D
                if ( ( coords ) or (numComponents == 2) ) {
                    // write a vec3 and fill with zeros
                    eigenX::VectorSd<3> q;
                    corlib::detail_::padWithZeros(  q, quantity );
                    corlib::detail_::writeQuantity( vtu_, q );
                }
                else corlib::detail_::writeQuantity( vtu_, quantity );
                vtu_ << "\n";
            }
        }
        vtu_ << "        </DataArray>" << std::endl;
        return;
    }

    //--------------------------------------------------------------------------
    //! this function writes the final tags to close the file
    void openCellData( )
    {
        vtu_ << "      <CellData>" << std :: endl;
        return;
    }

    //--------------------------------------------------------------------------
    //! this function writes the final tags to close the file
    void closeCellData( )
    {
        vtu_ << "      </CellData>" << std :: endl;
        return;
    }

    //--------------------------------------------------------------------------
    /** Evaluate an element datum represented by a functor and a name. 
     *  Corresponding to the resolution_ variable, the cell quantity will be 
     *  evaluated at the cell's parametric coordinates
     *  \f[
     *         \xi = (i+1) / (r+1)
     *  \f]
     *  where \f$ i \f$ is the index of a sub-mesh element and \f$ r \f$ the
     *  order of resultion. E.g.: \f$ r=1 \to \xi = 1/2 \f$ and 
     *  \f$ r=2 \to \xi \in \{ 1/3, 2/3 \} \f$. 
     *  \tparam    EACC    Type of Element quantity Accessor
     *  \param[in] cellAcc Functor to access the cell datum
     *  \param[in] name    Name for the output datum
     */
    template<typename EACC>
    void writeCellQuantity( EACC & cellAcc, const std::string & name )
    {
        typedef typename EACC::result_type                       Datum;
        typedef typename corlib::detail_::AccessorTraits<Datum>  AccTraits;

        const unsigned    numComponents = AccTraits::size;
        typedef typename AccTraits::ValueType ValueType;
        const std::string nameOfType = corlib::detail_::VtkDataType<ValueType>::name();

        //-- Dirty hack follows here, because Paraview does not understand 2D --//
        const unsigned numVTUcomp = ( numComponents == 2 ? 3 : numComponents );

        vtu_ << "      <DataArray type=\"" << nameOfType 
             << "\" NumberOfComponents=\"" << numVTUcomp 
             << "\" Name=\""               << name 
             << "\" format=\""             << "ascii"
             << "\">"
             << std :: endl;
        
        // numbers of cells per direction
        const unsigned nElementsCoarse = meshPtr_ -> numElements();

        // all cells in z-direction
        for ( unsigned f = 0; f < resolution_* nElementsCoarse; f ++ ) {
            
            const unsigned i = f / resolution_; // actual cell index

            const unsigned s = f - i * resolution_; // sub-cell index

            // evaluation point in sub-cell mesh
            typename Mesh::Element::VecLDim xi;
            xi[0] = 
                static_cast<double>( s + 1 ) / 
                static_cast<double>( resolution_ + 1);
                    
            // let cell do the evaluation
            // use policy in order to handle accessors which don't depend on xi
            Datum quantity = 
                detail_::ElementAccessPolicy<EACC>()( cellAcc, 
                                                      meshPtr_ -> getElementPointer( i ),
                                                      xi );


            // work-around for 2D
            if ( (numComponents == 2) ) {
                // write a vec3 and fill with zeros
                eigenX::VectorSd<3> q;
                corlib::detail_::padWithZeros(  q, quantity );
                corlib::detail_::writeQuantity( vtu_, q );
            }
            else corlib::detail_::writeQuantity( vtu_, quantity );
            vtu_ << std::endl;
        }
        vtu_ << "        </DataArray>" << std::endl;
    }

private:
    //--------------------------------------------------------------------------
    Mesh *              meshPtr_; //!< Access to mesh
    std::ostream &          vtu_; //!< Output stream
    const int        resolution_; //!< Element-resolution
};

#endif
