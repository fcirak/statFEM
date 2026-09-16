//
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//
//                 Computational Structural Mechanics Lab
//                         University of Cambridge
//
// This software is copyrighted and all rights are retained by the CSMLab. All
// use, disclosure, and/or reproduction of any part not expressly authorized by
// F Cirak is prohibited. (C) 2010.
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//

#ifndef subdiv_surf_vertex_h
#define subdiv_surf_vertex_h

#include <map>
#include <algorithm>
#include <iostream>
#include <limits>

#include <Eigen/Core>
#include <corlib/eigenX.hpp>

#include <corlib/verify.hpp>

#include <subdiv/surf/MeshTags.hpp>

//------------------------------------------------------------------------------
namespace subdiv{
    namespace surf{

        class Vertex;     
        namespace eigenX = corlib::eigenX;
    }
}

//------------------------------------------------------------------------------
/// Vertex of refined subdivision mesh.
///
/// The vertex primarily holds a spatial coordinate. If it is refined than a
/// a coordinate for each refinement level is stored in a map, too.
/// 
/// To establish the subdivision shape functions using an picking/extraction
/// matrix approach, additional "coordinates" can be supplied.
class subdiv::surf::Vertex
{
public:
    /// Vertex tag type
    typedef enum subdiv::surf::vertexTag                            VertexTag;

    /// Dimension of embedding space
    static const unsigned dim = 3;

    /// Vector of first dim unknowns (coordinates)
    typedef eigenX::VectorSd<dim>     VecDim;
    /// Vector of all coordinates/variables #dim + #NUMADDCOORD
    typedef Eigen::VectorXd           Point;

    /// maps refinement-level (int) to coordinate (Point),
    /// sorted from largest to smallest level of refinement
    typedef std::map< int, Point, std::greater< int > >             LevPoint;
    typedef LevPoint::iterator                                      LevPointIter;

public:
    /// Constructor
    Vertex( int level = 0 ) :
        index_( std::numeric_limits< unsigned >::max( ) ),
        tag_( subdiv::surf::VERTEX_NOTAG )
    {
        if ( levelPoints_.find( level ) == levelPoints_.end( ) ) {
            levelPoints_[ level ] =
                Vertex::makePoint( std::numeric_limits< double >::max() );
        }
    }
    
    /// @name Methods for vector of generalised coordinates/variables
    //@{
    /// Get number of vertex variables in addition to dim=3 coordinates
    static
    unsigned getNumAddCoords( ) { return (pointSize_ - dim); }

    /// Set number of vertex variables in addition to dim=3 coordinates
    ///
    /// Warning: This method does not change the length of the #Point vectors
    ///          of these #Vectors_, but new ones will be created to the new length.
    static
    void setNumAddCoords( unsigned nac ) { pointSize_ = dim + nac; }

    /// Make and return a generalised coordinates/variables vector
    static
    Point makePoint( double val = 0.0 )
    {
        FTL_VERIFY( pointSize_ >= dim );
        return Eigen::VectorXd::Constant( pointSize_, val );
    }
    //@}

    /// @name iterators
    //@{
    /// Return iterator to beginning of Points map
    LevPointIter levBegin( ) { return levelPoints_.begin( ); }
    /// Return iterator to end of Points map
    LevPointIter levEnd( ) { return levelPoints_.end( ); }
    //@}

    /// @name IO
    //@{
    /// Read in the #dim=3 coordinates
    friend
    void operator>>( std::istream & is, subdiv::surf::Vertex * v )
    {
        Point & coor = v->levelPoints_[ 0 ];
        for ( unsigned d=0; d<dim; ++d ) is >> coor[ d ];
    }

    /// Write out the #dim=3 coordinates
    friend
    std::ostream & operator<<( std::ostream & os, subdiv::surf::Vertex * v )
    {
        const Point & coor = v->levelPoints_.begin()->second;
        for ( unsigned d=0; d<dim; ++d ) os << " " << coor[ d ];

        return os;
    }
    //@}
    
    /// Return Index
    unsigned index( ) const { return index_; } 
    /// Set index
    void setIndex( unsigned i ) { index_ = i; }
    

    /// Get the vertex tag
    const VertexTag & tag( ) const { return tag_; }
    /// Set the vertex tag
    void setTag( const VertexTag tg ) { tag_ = tg; }

    /// Add new refinement level
    ///
    /// \param[in]   newCoord   New coordinates/variables of the vertex
    /// \param[in]   level      Refinement level
    void addNewLevel( const Point & newCoord, const int level ) 
    {
        // add or overwrite vertex coordinate at level
        levelPoints_[ level ] = newCoord;
        return;
    }
    
    /// Return max active refinement level of the vertex
    int maxActiveLevel( ) { return levelPoints_.begin()->first; }

    /// Return vertex coordinates/variables at given level
    const Point & point( const int level )
    { 
        FTL_VERIFY( levelPoints_.find( level ) != levelPoints_.end( ) );
        return levelPoints_[ level ];
    }

    /// Return vertex coordinates (dim=3) at maximum refinement level
    VecDim giveCoordinates( ) const
    { 
        return (levelPoints_.begin()->second).head( dim );
    }

protected:
    /// Remove level
    void removeLevel( int level )
    {
        LevPointIter itr = levelPoints_.find( level );
        assert( itr != levelPoints_.end( ) );
        levelPoints_.erase( itr );
        return;
    }

    /// Remove maximum refinement level
    void removeMaxLevel( )
    {
        FTL_VERIFY( levelPoints_.size( ) );
        levelPoints_.erase( levelPoints_.begin() );
        return;
    }

    /// Set vertex coordinate at given refinement level
    void setPoint( const Point & p, int level = 0 )
    { 
        FTL_VERIFY( levelPoints_.find( level ) != levelPoints_.end( ) );
        levelPoints_[ level ] = p;
        return;
    }
    
protected:
    /// Number of generalised coordinates of the vertex containing the #dim=3
    /// coordinates followed by additional coefficients, which are subdivided
    /// like the physical coordinates.
    ///
    /// pointSize_ should be the same for all vertices in the mesh
    static unsigned pointSize_;

private:
    unsigned          index_;        ///< vertex index
    VertexTag         tag_;          ///< vertex tag
    LevPoint          levelPoints_;  ///< point at different levels (if refined)
};

//------------------------------------------------------------------------------
unsigned subdiv::surf::Vertex::pointSize_ = subdiv::surf::Vertex::dim;
    
#endif 
