// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//
//                 Computational Structural Mechanics Lab
//                         University of Cambridge
//
// This software is copyrighted and all rights are retained by the University of
// Cambridge. All use, disclosure, and/or reproduction of any part not expressly
// authorised by the University of Cambridge is prohibited. (C) 2010
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

//! @file NodalForces.hpp

#ifndef solid_fem_nodalforces_h
#define solid_fem_nodalforces_h
//------------------------------------------------------------------------------
#include <map>
#include <functional>
#include <boost/numeric/ublas/vector.hpp>
#include <corlib/PropertiesParser.hpp>
//------------------------------------------------------------------------------
namespace solid{
    namespace fem{
        template< typename NODE > class NodalForces;
    }
}


//------------------------------------------------------------------------------
/** \brief Application of nodal forces
 *
 *  \details
 *  Store a map containing node IDs and a nodal force vectors attached to it.
 *  The map can be build either by looking into a file or by a function.
 *  The file contains a list of node IDs and force vectors components.
 *  The function returns the nodal force considering the nodal co-ordinate.
 *
 *  The application of these forces is achieved
 *  by the implemented function call operator().
 *
 *  \tparam NODE  Type of node
 */
template< typename NODE >
class solid::fem::NodalForces
    : public std::binary_function< NODE *, const double, void >
{
public:
    //! DOF vector at node, i.e. here nodal load vector type
    typedef boost::numeric::ublas::bounded_vector< double, NODE::dof >  VecDof;

public:
    //! Empty constructor
    NodalForces( )
    { }

    //! Read nodal force data from an input stream
    NodalForces( std::istream & inp )
    {
        this->add( inp );
        return;
    }

    //! Apply nodal forces according to a functor
    template< typename MESH, typename FUNC >
    NodalForces( MESH & mesh, FUNC func )
    {
        this->add( mesh, func );
        return;
    }

    //! Cleanup
    ~NodalForces() { nodalForces_.clear(); }

    //! Add nodal forces contained in file
    //!
    //! Reads the data for nodal forces from an input stream of the form
    //! <pre>
    //!      number of nodalForces
    //!      nodeId direction value
    //!      ...
    //! </pre>
    //! and stores the data in the map
    //!
    //! \param[in]  inp   input stream
    void add( std::istream & inp )
    {
        corlib::skip_comment( inp );

        // read number of forces
        unsigned nForces = 0;
        inp >> nForces;

        // read all forces from stream
        for ( unsigned f = 0; f < nForces; f ++  ) {

            corlib::skip_comment( inp );

            unsigned id, dir;
            double val;
            inp >> id >> dir >> val;
            
            VecDof nf;
            nf.clear();
            nf(dir) = val;

            // insert new force or add upon existing force
            this->add( id, nf );
        }
    }

    //! Add nodal forces described with function
    //!
    //! \param[in] mesh   Mesh containing nodes
    //! \param[in] func   Nodal force function, i.e. provided a nodal co-ordinate
    //!                   it returns true if the co-ordinate is
    //!                   subject to nodal force. The nodal force is returned as well.
    //!
    //! \tparam    MESH   Mesh type
    //! \tparam    FUNC   Function type
    template< typename MESH, typename FUNC >
    void add( MESH & mesh, FUNC func )
    {
        typename MESH::NodeConstIterator first = mesh.nodesBegin();
        typename MESH::NodeConstIterator  last = mesh.nodesEnd( );
        for( ; first != last; ++ first ) {
            typename NODE::VecDim x = (*first) -> giveCoordinates();
            VecDof nodalForce; 
            nodalForce.clear();
            const bool subjectToForce = func( x, nodalForce );

            if ( subjectToForce ) {
                const unsigned nodeId = (*first) -> giveId();
                this->add( nodeId, nodalForce );
            }
        }
    }

    //! Add nodal force by node ID and load vector
    //!
    //! \param[in]  nodeId            The node ID
    //! \param[in]  nodalForce        Nodal load vector
    //! \param[in]  overWriteValues   Instead of adding contrib. at recurring node,
    //!                               last contrib. is taken
    void add( const unsigned nodeId, const VecDof nodalForce,
              const bool overWriteValues = false )
    {
        typename std::map<unsigned,VecDof>::iterator nfIter =
            nodalForces_.find( nodeId );
        if ( nfIter == nodalForces_.end() ) {
            nodalForces_[ nodeId ] = nodalForce;
        } else {
            if ( overWriteValues )
                (nfIter->second) = nodalForce;
            else 
                (nfIter->second) += nodalForce;
        }
        return;
    }

    //! Apply nodal forces to nodes (multiplied by factor)
    void operator()( NODE * np, const double factor ) const
    {
        const unsigned id = np -> giveId();

        typename std::map<unsigned,VecDof>::const_iterator nfIter = 
            nodalForces_.find( id );
        if ( nfIter != nodalForces_.end() ) {

            VecDof nf = factor * (nfIter -> second);

            np -> addToForce( nf );
        }

        return;
    }

    std::ostream & write( std::ostream & out ) const
    {
        typename std::map<unsigned,VecDof>::const_iterator iter = nodalForces_.begin();
        for ( ; iter != nodalForces_.end(); ++ iter ) {
            out << iter -> first << ": " << iter -> second << std::endl;
        }
        return out;
    }

private:
    //! list of nodal load vectors
    std::map< unsigned, VecDof >        nodalForces_;
};

#endif
