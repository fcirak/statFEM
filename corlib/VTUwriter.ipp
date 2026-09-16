// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//
//                 Computational Structural Mechanics Lab
//                         University of Cambridge
//
// This software is copyrighted and all rights are retained by the University of
// Cambridge. All use, disclosure, and/or reproduction of any part not expressly
// authorised by the University of Cambridge is prohibited. (C) 2010
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

//! @file VTUwriter.ipp
#include <ostream>
#include <functional>
#include <boost/functional.hpp>

#include <corlib/Shape.hpp>
#include <corlib/verify.hpp>
#include <corlib/Accessor.hpp>
#include <corlib/eigenX.hpp>

//------------------------------------------------------------------------------
//! this function writes the header and the mesh 
template<typename MESH> 
void corlib::VTUwriter<MESH>::writeMesh( )  
{
    const unsigned nNodes    = this -> setNumberOfNodes();
    const unsigned nElements = this -> setNumberOfElements();

    // write header of VTU file
    vtu_ << "<?xml version=\"1.0\"?>" << std :: endl
         << "<VTKFile type=\"UnstructuredGrid\" byte_order=\"LittleEndian\">" << std :: endl
         << "  <UnstructuredGrid>" << std :: endl
         << "  <Piece NumberOfPoints=\"" << nNodes 
         << "\"  NumberOfCells=\"" << nElements 
         << "\">" << std :: endl;

    // write nodal coordinates
    vtu_ << "    <Points>" << std :: endl;

    this -> writeNodalQuantity( corlib::accessorFun( &MESH::Node::giveCoordinates, "Coordinates"),
                                true );

    vtu_  << "    </Points>" << std :: endl;

    // write cells
    vtu_ << "    <Cells>" << std :: endl;
    
    // prepare map for global2local translation
    std::vector<unsigned> globalIndices;
    unsigned numNodes = meshPtr_ -> numNodes();
    for ( unsigned i = 0; i < numNodes; ++i ) {
        globalIndices.push_back( ( meshPtr_ -> getNodePointer( i ) ) -> giveId() );
    }
    // make a map from the vector
    std::map<unsigned,unsigned> global2LocalMap;
    for ( unsigned ind = 0; ind < globalIndices.size(); ++ind ) {
        global2LocalMap.insert( std::make_pair( globalIndices[ind], ind ) );
    }
    globalIndices.clear();

    corlib::detail_::ElementConnectivity<typename MESH::Element> ec( &global2LocalMap );
    this -> writeElementQuantity( ec );
    global2LocalMap.clear();
   
    corlib::detail_::ElementOffset<typename MESH::Element> eo;
    this -> writeElementQuantity( eo );

    corlib::detail_::VtkTypeFunctor<typename MESH::Element> ev;
    this -> writeElementQuantity( ev );
    
    vtu_ << "    </Cells>" << std :: endl;

    return;
}

//------------------------------------------------------------------------------
//! set the number of nodes (overloadable)
template< typename MESH > 
unsigned corlib::VTUwriter<MESH>::setNumberOfNodes()
{
    return meshPtr_ -> numNodes();
}

//------------------------------------------------------------------------------
//! set the number of nodes (overloadable)
template< typename MESH > 
unsigned corlib::VTUwriter<MESH>::setNumberOfElements()
{
    return meshPtr_ -> numElements();
}


//------------------------------------------------------------------------------
//! this function writes the final tags to close the file
template< typename MESH > 
void corlib::VTUwriter<MESH>::finalize( )
{
    vtu_ << "  </Piece>" << std :: endl
         << "  </UnstructuredGrid>" << std :: endl;
    vtu_ << "</VTKFile>" << std :: endl;
    return;
}

//------------------------------------------------------------------------------
//! this function writes the final tags to close the file
template< typename MESH > 
void corlib::VTUwriter<MESH>::openPointData( )
{
    vtu_ << "    <PointData>" << std :: endl;
    return;
}

//------------------------------------------------------------------------------
//! this function writes the final tags to close the file
template< typename MESH > 
void corlib::VTUwriter<MESH>::closePointData( )
{
    vtu_ << "    </PointData>" << std :: endl;
    return;
}


//------------------------------------------------------------------------------
//! this function writes any nodal quantity accessed by a Node functor
template< typename MESH > 
template< typename NFUN >
void corlib::VTUwriter<MESH>::writeNodalQuantity( NFUN nodeFunctor, const bool coords )
{
    // get name and amount of quantities 
    const std::string quantityName  = nodeFunctor.name();
    const unsigned    numComponents = NFUN::numComponents;

    typedef typename NFUN::ValueType ValueType;
    const std::string nameOfType = corlib::detail_::VtkDataType<ValueType>::name();

    //-- Dirty hack follows here, because Paraview does not understand 2D --//
    unsigned numVTUcomp = ( numComponents == 2 ? 3 : numComponents );
    if( coords ) numVTUcomp = 3;
    //---------------------------------------------------------------------//


    vtu_ << "      <DataArray type=\"" << nameOfType 
         << "\" NumberOfComponents=\"" << numVTUcomp 
         << "\" Name=\""               << quantityName 
         << "\" format=\""             << dataFormat_ 
         << "\">"
         << std :: endl;

    const unsigned nNodes = this -> setNumberOfNodes();

    // go through all nodes of the mesh
    for ( unsigned n = 0; n < nNodes; n ++ ) {
        // obtain node pointer
        const typename MESH::Node* nodePtr = meshPtr_ -> getNodePointer( n );
        // obtain nodal quantity
        typename NFUN::ReturnType quantity = nodeFunctor( nodePtr );
        //----------------------------------------------------------------------
        // write the quantity
        if ( coords or (numComponents == 2) ) {
            // write a vec3 and fill with zeros
            eigenX::VectorSd<3> q; 
            detail_::padWithZeros( q, quantity );
            detail_::writeQuantity( vtu_, q );
        }
        else detail_::writeQuantity( vtu_, quantity );

        vtu_ << std::endl;
    }
    vtu_ << "      </DataArray>" << std :: endl;

    return;
}

//------------------------------------------------------------------------------
//! this function writes the final tags to close the file
template< typename MESH > 
void corlib::VTUwriter<MESH>::openCellData( )
{
    vtu_ << "    <CellData>" << std :: endl;
    return;
}

//------------------------------------------------------------------------------
//! this function writes the final tags to close the file
template< typename MESH > 
void corlib::VTUwriter<MESH>::closeCellData( )
{
    vtu_ << "    </CellData>" << std :: endl;
    return;
}


//------------------------------------------------------------------------------
//! this function writes any element quantity accessed by a Element functor
template< typename MESH > 
template< typename EFUN >
void corlib::VTUwriter<MESH>::writeElementQuantity( EFUN elementFunctor )
{
    // get name and amount of quantities 
    const std::string quantityName  = elementFunctor.name();
    //const std::string quantityName  = elementFunctor.name();
    const unsigned    numComponents = EFUN::numComponents;

    typedef typename EFUN::ValueType ValueType;
    const std::string nameOfType = corlib::detail_::VtkDataType<ValueType>::name();

    //-- Dirty hack follows here, because Paraview does not understand 2D --//
    const unsigned    numVTUcomp = ( numComponents == 2 ? 3 : numComponents );
    
    vtu_ << "      <DataArray type=\"" << nameOfType 
         << "\" NumberOfComponents=\"" << numVTUcomp 
         << "\" Name=\""               << quantityName 
         << "\" format=\""             << dataFormat_ 
         << "\">"
         << std :: endl;

    const unsigned nElements = this -> setNumberOfElements();

    // expect some Eigen::Matrix type of quantity
    for( unsigned e = 0; e < nElements; e ++ ) {
        // obtain element quantity
        const typename MESH::Element * elPtr = meshPtr_ -> getElementPointer( e );
        typename EFUN::ReturnType quantity = elementFunctor( elPtr );
        if ( numComponents ==  2 ) {
            eigenX::VectorSd<3> q; 
            detail_::padWithZeros( q, quantity );
            detail_::writeQuantity( vtu_, q );
        }
        else detail_::writeQuantity( vtu_, quantity );
        vtu_ << std::endl;
    }
    
    vtu_ << "      </DataArray>" << std :: endl;

    return;
}


//------------------------------------------------------------------------------
// Write primitive quantities
namespace corlib{
    namespace detail_{
        std::ostream & writeQuantity( std::ostream & out, 
                                      const eigenX::VectorSd<1> & quantity )
        {
            writeQuantity( out, quantity[0] );
            return out;
        }

        std::ostream & writeQuantity( std::ostream & out, 
                                      const eigenX::VectorSd<2> & quantity )
        {
            for ( unsigned d1 = 0; d1 < 2; d1 ++ ) {
                writeQuantity(out, quantity(d1));
            }
            return out;
        }
    
        std::ostream & writeQuantity( std::ostream & out, 
                                      const eigenX::VectorSd<3> & quantity )
        {
            for ( unsigned d1 = 0; d1 < 3; d1 ++ ) {
                writeQuantity(out, quantity(d1));
            }
            return out;
        }
    
        std::ostream & writeQuantity( std::ostream & out, 
                                      const eigenX::MatrixSd<3> & quantity )
        {
            for ( unsigned d1 = 0; d1 < 3; d1 ++ ) {
                for ( unsigned d2 = 0; d2 < 3; d2 ++ ) {
                    writeQuantity( out, quantity( d1, d2 ) );
                }
            }
            return out;
        }

        std::ostream & writeQuantity( std::ostream & out, 
                                      const std::vector<unsigned> & quantity )
        {
            std::copy( quantity.begin(), quantity.end(), std::ostream_iterator<unsigned>( out, " " ) );
            return out;
        }

        std::ostream & writeQuantity( std::ostream & out, 
                                      const double & quantity )
        {
            out << static_cast<float>( quantity ) << " ";
            return out;
        }

        std::ostream & writeQuantity( std::ostream & out, 
                                      const unsigned int & quantity )
        {
            out << quantity << " ";
            return out;
        }

        std::ostream & writeQuantity( std::ostream & out, 
                                      const int & quantity )
        {
            out << quantity << " ";
            return out;
        }

        std::ostream & writeQuantity( std::ostream & out, 
                                      const bool & quantity )
        {
            out << std::noboolalpha << quantity << " ";
            return out;
        }
    }
}

//==============================================================================//
// Stuff for writing a pvd-container file                                       //
//==============================================================================//

//------------------------------------------------------------------------------
//! set up a name for the output file given the current time
std::string corlib::VTUanim::snapshotName( const double time, const unsigned stepNum )
{
    corlib::UniqueFilename filenameGen( 4 );

    // generate the filename
    std::string snapshot = filenameGen( basename_, stepNum, suffix_ );
    // store the name for the final pvd-file
    animation_.push_back( std::make_pair( time, snapshot ) );
    return snapshot;
}

//------------------------------------------------------------------------------
//! write the file containing the timestamps and snapshot file names
std::ostream & corlib::VTUanim::writeAnimationFile( std::ostream & vta ) const
{
    // write header
    vta << "<?xml version=\"1.0\"?>" << std::endl
        << "<VTKFile type=\"Collection\" version=\"0.1\" "
        << "byte_order=\"LittleEndian\">" << std::endl
        << "  <Collection>" << std::endl;

    // go through vector of snapshots
    std::vector< Snapshot_ > :: const_iterator siter = animation_.begin( );
    std::vector< Snapshot_ > :: const_iterator send  = animation_.end( );
    for( ; siter != send; ++ siter ) {
        vta << "    <DataSet timestep=\"" << siter -> first 
            << "\" file=\"" << siter -> second << "\"/>" << std::endl;
    }
        
    // finish file
    vta << "  </Collection>" << std::endl
        << "</VTKFile>" << std::endl;

    return vta;
}

//------------------------------------------------------------------------------
//! Store filenames together with current time and number of part in container
void corlib::VTUparaAnim::storeSnapshots( const double time, const unsigned stepNum )
{
    corlib::UniqueFilename filenameGen( 3, 4 );
    
    for ( unsigned p = 0; p < numProcesses_; p ++ ) {
        std::string snapshot = filenameGen( basename_, p, stepNum, suffix_ );
        paraAnimation_.push_back( boost::make_tuple( time, p, snapshot ) );
    }
    return;
}

//------------------------------------------------------------------------------
//! Store filenames together with current time and number of part in container
//! In addition to previous, get actual basename as argument - this allows embedding VTU files into 
//! directories based on on timesteps.
void corlib::VTUparaAnim::storeSnapshots( std::string basename, const double time, const unsigned stepNum )
{
    corlib::UniqueFilename filenameGen( 3, 4 );
    
    for ( unsigned p = 0; p < numProcesses_; p ++ ) {
        std::string snapshot = filenameGen( basename, p, stepNum, suffix_ );
        paraAnimation_.push_back( boost::make_tuple( time, p, snapshot ) );
    }
    return;
}

//------------------------------------------------------------------------------
//! write the file containing the timestamps and snapshot file names
std::ostream & corlib::VTUparaAnim::writeAnimationFile( std::ostream & vta ) const
{
    // write header
    vta << "<?xml version=\"1.0\"?>" << std::endl
        << "<VTKFile type=\"Collection\" version=\"0.1\" "
        << "byte_order=\"LittleEndian\">" << std::endl
        << "  <Collection>" << std::endl;

    // go through vector of snapshots
    std::vector< ParaSnapshot_ > :: const_iterator siter = paraAnimation_.begin( );
    std::vector< ParaSnapshot_ > :: const_iterator send  = paraAnimation_.end( );
    for( ; siter != send; ++ siter ) {
        vta << "    <DataSet timestep=\"" << siter -> get<0>()
            << "\" part=\"" << siter -> get<1>() 
            << "\" file=\"" << siter -> get<2>() << "\"/>" << std::endl;
    }
        
    // finish file
    vta << "  </Collection>" << std::endl
        << "</VTKFile>" << std::endl;

    return vta;
}

//------------------------------------------------------------------------------
//! Simply write the pieces for a static computation
std::ostream & corlib::VTUparaAnim::writeGroupFile( std::ostream & vtg ) const
{
    // write header
    vtg << "<?xml version=\"1.0\"?>" << std::endl
        << "<VTKFile type=\"Collection\" version=\"0.1\" "
        << "byte_order=\"LittleEndian\">" << std::endl
        << "  <Collection>" << std::endl;
    
    for ( unsigned p = 0; p < numProcesses_; p++ ) {
        std::string name = corlib::UniqueFilename(3)( basename_, p, suffix_ );
        vtg << "    <DataSet  part=\"" << p
            << "\" file=\"" << name << "\"/>" << std::endl;
    }

    // finish file
    vtg << "  </Collection>" << std::endl
        << "</VTKFile>" << std::endl;

    return vtg;
}
