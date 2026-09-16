// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//
//                 Computational Structural Mechanics Lab
//                         University of Cambridge
//
// This software is copyrighted and all rights are retained by the University of
// Cambridge. All use, disclosure, and/or reproduction of any part not expressly
// authorised by the University of Cambridge is prohibited. (C) 2010
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

//! @file VTUwriter.hpp

#ifndef corlib_vtuwriter_h
#define corlib_vtuwriter_h
//------------------------------------------------------------------------------
#include <iomanip>
#include <sstream>
#include <string>
#include <map>

#include <boost/tuple/tuple.hpp>
#include <boost/function.hpp>

#include <corlib/Shape.hpp>
#include <corlib/ShapeToVtk.hpp>
#include <corlib/UniqueFilename.hpp>
#include <corlib/eigenX.hpp>

//------------------------------------------------------------------------------
namespace corlib {

    template<typename MESH>  class VTUwriter;
    class VTUanim;
    class VTUparaAnim;
    
    template<typename ELEMENT> class AssignConstantNumber;

    namespace detail_{
        //----------------------------------------------------------------------
        //! \brief Provide names for the data types 
        template<typename T> struct VtkDataType;
        //! \cond SKIPDOX
        template<> struct VtkDataType<double>  
        { 
            static const std::string name() { return std::string( "Float64" ); } 
        };

        template<> struct VtkDataType<unsigned>
        { 
            static const std::string name() { return std::string( "Int32" ); } 
        };

        template<> struct VtkDataType<int>
        { 
            static const std::string name() { return std::string( "Int32" ); } 
        };

        template<> struct VtkDataType<bool>
        {
            static const std::string name() { return std::string( "Int32" ); } 
        };
        //! \endcond

        //----------------------------------------------------------------------
        //! \brief change the number of output nodes for some element types
        template<shape SHAPE, unsigned numPoints> struct VtkOutputNodes
        { 
            static const unsigned vtkTypeNum = VtkType< SHAPE, numPoints>::vtk;
            
            static const unsigned num = VtkTypeNodeNum<vtkTypeNum>::np;	
        };

        //----------------------------------------------------------------------
        //! \brief Convenience functor for vtk-types of the elements
        template<typename ELEMENT> 
        class VtkTypeFunctor 
            : public boost::function<unsigned(const ELEMENT*)>
        {
        public:
            static const std::string name() { return std::string("types"); } 
            static const unsigned numComponents = 1;

            typedef unsigned  ReturnType;
            typedef unsigned   ValueType;
        
            ReturnType operator()( const ELEMENT * ep ) const
            {
                return VtkType<ELEMENT::myShape, ELEMENT::numNodes>::vtk;
            }
        };

        //--------------------------------------------------------------------------
        //! \brief Convenience functor for accessing the element's node indices
        template<typename ELEMENT>
        class ElementConnectivity 
            : public boost::function<std::vector<unsigned>(const ELEMENT*)>
        {
        public:
            ElementConnectivity( std::map<unsigned, unsigned> * global2LocalMap ) 
                : global2LocalMap_( global2LocalMap ) { }

    
            static const std::string name() { return std::string("connectivity"); } 
            static const unsigned numComponents = 1;
        
            static const unsigned numOutputNodes = 
                corlib::detail_::VtkOutputNodes<ELEMENT::myShape,ELEMENT::numNodes>::num;
        
            typedef std::vector<unsigned>                     ReturnType;
            typedef unsigned                                  ValueType;
        
            ReturnType operator() (const ELEMENT* ep) 
            {
                ReturnType aux;
                typename ELEMENT::NodeConstIterator first = ep -> nodesBegin();
                typename ELEMENT::NodeConstIterator  last = ep -> nodesEnd();
                for ( ; first != last; ++ first ) {
                    aux.push_back( (*first) -> giveId() );
                }
                // truncate for some element types
                ReturnType connec( aux.begin(), aux.begin() + numOutputNodes );
                // translate connec to local numbering
                ReturnType::iterator iterB = connec.begin();
                ReturnType::iterator iterE = connec.end();
                for ( ; iterB != iterE; ++iterB ) {
                    std::map<unsigned, unsigned>::iterator pos = global2LocalMap_ -> find( *iterB );
                    FTL_VERIFY_DESCRIPTIVE( pos != global2LocalMap_ -> end(),
                                            "Index of global node %d not found within local nodes. \n",  
                                            *iterB );
                    // store local index 
                    *iterB = pos -> second;
                }
                return connec;
            }
        private:
            std::map<unsigned, unsigned> * global2LocalMap_;
        };

        //--------------------------------------------------------------------------
        //! Convenience functor for the element offsets
        template<typename ELEMENT>
        class ElementOffset 
            : public boost::function<unsigned(const ELEMENT*)>
        {
        public:
            static const std::string name() { return std::string("offsets"); } 
            static const unsigned numComponents = 1;

            typedef unsigned int              ReturnType;
            typedef unsigned int              ValueType;

            ElementOffset() 
                : index(0), 
                  offset( VtkOutputNodes<ELEMENT::myShape,ELEMENT::numNodes>::num )
            { }

            ReturnType operator() (const ELEMENT* ep) 
            {
                index += offset;
                return index;
            }

        private:
            unsigned index;
            const unsigned offset;
        };

        //----------------------------------------------------------------------
        template <typename T>
        void padWithZeros( eigenX::VectorSd<3> & out, const T & in )
        {
            FTL_VERIFY( false );
        }

        template <>
        void padWithZeros( eigenX::VectorSd<3> & out,
                           const eigenX::VectorSd<3> & in )
        {
            out = in;
        }

        template <>
        void padWithZeros( eigenX::VectorSd<3> & out,
                           const eigenX::VectorSd<2> & in )
        {
            out.setZero();
            out.head( 2 ) = in;
        }

        template <>
        void padWithZeros( eigenX::VectorSd<3> & out,
                           const eigenX::VectorSd<1> & in )
        {
            out.setZero();
            out.head( 1 ) = in;
        }

        //----------------------------------------------------------------------
        //! @name Free functions to write primitive quantities
        //@{
        std::ostream & writeQuantity( std::ostream & out, const double & quantity );
        std::ostream & writeQuantity( std::ostream & out, const eigenX::VectorSd<1> & quantity );
        std::ostream & writeQuantity( std::ostream & out, const eigenX::VectorSd<2> & quantity );
        std::ostream & writeQuantity( std::ostream & out, const eigenX::VectorSd<3> & quantity );
        std::ostream & writeQuantity( std::ostream & out, const eigenX::MatrixSd<3> & quantity );
        std::ostream & writeQuantity( std::ostream & out, const std::vector<unsigned> & quantity );
        std::ostream & writeQuantity( std::ostream & out, const unsigned int & quantity );
        std::ostream & writeQuantity( std::ostream & out, const int & quantity );
        std::ostream & writeQuantity( std::ostream & out, const bool & quantity );
        //@}



    }// namespace detail_
}// namespace corlib

//------------------------------------------------------------------------------
//! class which manages data written as a VTU file (for Paraview, Mayavi2, etc.)
template< typename MESH > 
class corlib::VTUwriter
{
private:
    typedef MESH*                                                       MeshPtr_;
public:
    typedef boost::function< void( const typename MESH::Node * ) >      NodeFun;
    typedef NodeFun (MESH::*IterateOver)( NodeFun ) ;

public:
    VTUwriter( MeshPtr_  m, std::ostream & v  ) : meshPtr_( m ), vtu_( v ), 
                                                  dataFormat_( "ascii" ) { }
    virtual ~VTUwriter() { meshPtr_ = NULL; }

    // standard functions 
    virtual void writeMesh( );

    virtual void finalize(   );

    // here follow special functions for nodal and element quantities
    void  openPointData(  );
    void closePointData(  );
    template< typename NFUN > void writeNodalQuantity( NFUN nodeFunctor, const bool coord = false );

    void  openCellData(  );
    void closeCellData(  );
    template< typename EFUN > void writeElementQuantity( EFUN elementFunctor );

protected:
    virtual unsigned setNumberOfNodes();
    virtual unsigned setNumberOfElements();

protected:
    MeshPtr_            meshPtr_;
    std::ostream &      vtu_;

    std::string dataFormat_;
};

//------------------------------------------------------------------------------
/** \brief Convenience functor giving a constant number to each object
 *  \details This functor assigns a constant (!) to every object for the 
 *  output. Possible application is to pass the process number to the output.
 *  \tparam OBJ  Type of object which is given the number
 */
template<typename OBJ>
class corlib::AssignConstantNumber
    : public boost::function<unsigned(const OBJ*)>
{
public:
    static const unsigned numComponents = 1;
    
    typedef unsigned int              ReturnType;
    typedef unsigned int              ValueType;
    
    AssignConstantNumber( const unsigned & number, 
                          const std::string & name ) 
        : number_( number ), name_( name ) { }
    
    ReturnType operator() (const OBJ* ep) 
    {
        return number_;
    }
    
    std::string name() { return name_; }
    
private:
    const unsigned    number_; //!< Constant number to assign each object
    const std::string name_;   //!< Name of the field
};


//------------------------------------------------------------------------------
//! object to manage animations
class corlib::VTUanim
{
public:
    VTUanim( const std::string & b, const std::string & s )
        : basename_( b ) , suffix_( s )  { }

    ~VTUanim( ) { animation_.clear( ); }
    
    std::string snapshotName( const double time, const unsigned stepNum );

    std::ostream & writeAnimationFile( std::ostream & vta ) const;

private:
    std::string                               basename_;  //!< output file base
    std::string                               suffix_;    //!< output file end
    typedef std::pair< double, std::string >  Snapshot_;  
    std::vector< Snapshot_ >                  animation_; //!< animation file
};


//------------------------------------------------------------------------------
//! object to manage animations with parts (parallel comp.) per time step
class corlib::VTUparaAnim
{
public:
    VTUparaAnim( const std::string & b, 
                 const std::string & s, 
                 const unsigned numProcesses ) 
        : basename_( b ) , 
          suffix_( s ), 
          numProcesses_( numProcesses )
        { }

    ~VTUparaAnim( ) { paraAnimation_.clear( ); }
    
    void storeSnapshots( const double time, const unsigned stepNum );
    void storeSnapshots( std::string basename, const double time, const unsigned stepNum );
    std::ostream & writeAnimationFile( std::ostream & vta ) const;
    std::ostream & writeGroupFile( std::ostream & vtg ) const;
    std::string getBasename() const { return basename_; }

private:
    std::string                               basename_;  //!< output file base
    std::string                               suffix_;    //!< output file end

    typedef boost::tuple< double, unsigned, std::string > ParaSnapshot_;
    std::vector< ParaSnapshot_ >              paraAnimation_;
    const unsigned                            numProcesses_;
};


#include "VTUwriter.ipp"
//------------------------------------------------------------------------------
#endif
