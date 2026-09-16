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

//! @author Burkhard Bornemann
//! @date   06/2012

#ifndef aux_parameters_h
#define aux_parameters_h

#include <corlib/PropertiesParser.hpp>

//------------------------------------------------------------------------------
namespace aux{

    // make ublas namespace available in namespace
    namespace ublas = boost::numeric::ublas;

    class GridParameters;

}

//==============================================================================
/** Read and store input parameters */
class aux::GridParameters
{
public:
    typedef ublas::bounded_vector<double,3>          VecD3;
    typedef ublas::bounded_vector<int,3>        VecI3;

public:
    GridParameters( std::istream & inp )
    {
        corlib::PropertiesParser * prop = new corlib::PropertiesParser;
            
        dim = 0;
        prop -> registerPropertiesVar( "dim",               dim );

        bbmin = ublas::scalar_vector<double>( 3, std::numeric_limits<double>::max() );
        prop -> registerPropertiesVar( "bbminX",            bbmin[0] );
        prop -> registerPropertiesVar( "bbminY",            bbmin[1] );
        prop -> registerPropertiesVar( "bbminZ",            bbmin[2] );

        bbmax = ublas::scalar_vector<double>( 3, -std::numeric_limits<double>::max() );
        prop -> registerPropertiesVar( "bbmaxX",            bbmax[0] );
        prop -> registerPropertiesVar( "bbmaxY",            bbmax[1] );
        prop -> registerPropertiesVar( "bbmaxZ",            bbmax[2] );

        numCells.clear();
        prop -> registerPropertiesVar( "numCellsX",         numCells[0] );
        prop -> registerPropertiesVar( "numCellsY",         numCells[1] );
        prop -> registerPropertiesVar( "numCellsZ",         numCells[2] );
            
        prop -> registerPropertiesVar( "numParts",          numParts );
        overlapLayers = 0;
        prop -> registerPropertiesVar( "overlapLayers",     overlapLayers );

        writeVTU = false;
        prop -> registerPropertiesVar( "writeVTU",          writeVTU );
            
        // read variables from the input.dat file
        prop -> readValues( inp );
        delete prop;

        cellSizes.clear();
        for ( unsigned d = 0; d < dim; ++d )
            cellSizes[d] = ( bbmax[d] - bbmin[d] ) / static_cast<double>( numCells[d] );

        // set higher dimensions to zero (in case they are contained in input file) --- *crucial*
        std::fill( numCells.begin()+dim, numCells.end(), 0 );
    }
        
public: // public on purpose!

    /// Spatial dimension (1, 2 or 3)
    unsigned       dim;

    /// Lower corner of domain
    ///
    /// Note: If you want to partition a Cartesian grid which is *not* mapped
    ///       to a cuboid in the physical domain, you can simply use the
    ///       zero indices here.
    VecD3           bbmin;

    /// Upper corner of domain
    ///
    /// Note: If you want to partition a Cartesian grid which is *not* mapped
    ///       to a cuboid in the physical domain, you can simply use #numCells
    ///       here.
    VecD3           bbmax;
        
    /// Number of cells in each space direction
    ///
    /// Note: They are assumed 0 for higher dimensions.
    VecI3           numCells;

    /// Size of cells in each space dimension
    ///
    /// Note: If you want to partition a Cartesian grid which is *not* mapped
    ///       to a cuboid in the physical domain, you can simply use (1,1,1) here.
    VecD3           cellSizes;

    /// Number of partitions
    int             numParts;
    /// Number of overlapping (aka ghost) cells in one direction
    int             overlapLayers;

    /// Write VTU file
    bool            writeVTU;
};

#endif
