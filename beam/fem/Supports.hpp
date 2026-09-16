// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//
//                 Computational Structural Mechanics Lab
//                         University of Cambridge
//
// This software is copyrighted and all rights are retained by the University of
// Cambridge. All use, disclosure, and/or reproduction of any part not expressly
// authorised by the University of Cambridge is prohibited. (C) 2010
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

//! @file Supports.hpp

#ifndef beam_fem_supports_h
#define beam_fem_supports_h

#include <algorithm>
#include <map>
#include <cctype>
#include <string>

#include <Eigen/Core>

#include <corlib/eigenX.hpp>
#include <corlib/Constraints.hpp>

//------------------------------------------------------------------------------
// forward declarations
namespace beam {
    namespace fem {

        class Support;
        class SupportsFromFile;

        template< typename LBASIC, unsigned DEGREE >
        class SupportsAtNode;

        template< typename LBASIC >
        class SupportsAtNode< LBASIC, 3 >;

        namespace eigenX = corlib::eigenX;

        //! Types of supports
        enum supportType {
            UNDEFINED,  //!< Nothing -> don't use it
            DISP,       //!< General displacement constraint
            ROT,        //!< General rotation constraint
            DISPX,      //!< x-axis displacement constraint
            DISPY,      //!< y-axis displacement constraint
            ROTZ,       //!< rotation around z-axis constraint
            PIN         //!< pin-point dof
        };

        namespace detail_{
            
            //! Convert string to supportType enum
            enum supportType stringToEnum( const std::string & input )
            {
                // make sure to have lower case only
                std::string name = input;
                std::transform (name.begin(), name.end(), name.begin(), ::tolower);

                if      ( name.find( "dispx" ) != std::string::npos ) return DISPX;
                else if ( name.find( "dispy" ) != std::string::npos ) return DISPY;
                else if ( name.find( "disp"  ) != std::string::npos ) return DISP;
                else if ( name.find( "rotz"  ) != std::string::npos ) return ROTZ;
                else if ( name.find( "rot"   ) != std::string::npos ) return ROT;
                else if ( name.find( "pin"   ) != std::string::npos ) return PIN;
                else FTL_VERIFY( false );
                return UNDEFINED;
            }

            //! convert supportType enum to string
            std::string enumToString( const supportType & type ) 
            {
                switch( type ) {
                case UNDEFINED: return "undefined"; break;
                case DISP:      return "disp";      break;
                case DISPX:     return "dispX";     break;
                case DISPY:     return "dispY";     break;
                case ROT:       return "rot";       break;
                case ROTZ:      return "rotZ";      break;
                case PIN:       return "pin";       break;
                default: return "undefined";
                }
                return "undefined";
            }

        }

    }
}

//==============================================================================

//==============================================================================
/** \brief Handle a displacement or rotation constraint along a given axis */
class beam::fem::Support
{
public:
    typedef eigenX::VectorSd<3>  Vec3;
    
    //! Constructor with data from stream
    Support( std::istream & inp ) { 

        // get type of support
        std::string aux;
        inp >> aux;
        type_ = detail_::stringToEnum( aux );
        
        // clear axis
        axis_.setZero();
        // check type
        if      ( type_ == DISPX ) axis_[0] = 1.;
        else if ( type_ == DISPY ) axis_[1] = 1.;
        else if ( type_ == ROTZ  ) axis_[2] = 1.; 
        else { // expect axis in inp
            for ( unsigned d = 0; d < 3; d ++ ) 
                inp >> axis_[d];
        }

        // normalise axis
        axis_ /= axis_.norm( );

        // read value
        inp >> value_;

    }

    Support( const supportType type, 
             const double value, 
             const Vec3 & axis )
        : type_(  type  ),
          value_( value ), 
          axis_(  axis  ) {
    }
    
    supportType giveType()  const { return type_; }
    double      giveValue() const { return value_; }
    Vec3        giveAxis()  const { return axis_;  }

    //! Print
    std::ostream & print( std::ostream & os ) const
    {
        os << "Type: " << detail_::enumToString( type_ )
           << ", value=" << value_
           << ", direction: " << axis_
           << std::endl;
        return os;
    }

private:
    enum supportType type_;
    double           value_;
    Vec3             axis_;
};

//==============================================================================
/** \brief Handle supports that are given by an input file
 *  \details This link handler is based on an input stream with the (old)
 *   format
 *   <pre>
 *   numberOfSupports
 *   nodeId1  flag1  value1
 *   nodeId2  flag2  value2
 *   ...      ...    ...
 *   </pre>
 *   (flag from {dispX,dispY,rotZ) with the number of supports (numberOfSupports) 
 *   given in the first line and each following line contains the node number 
 *   (nodeId), a flag and the value to apply. Alternatively, the more general 
 *   format can be used which looks like
 *   <pre>
 *   numberOfSupports
 *   nodeId1  flag1  axis1_1 axis1_2 axis1_3 value1
 *   nodeId2  flag2  axis2_1 axis2_2 axis2_3 value2
 *   ...      ...    ...
 *   </pre>
 *   using flag from {disp,rot} and the three coordinates of an axis. The axis
 *   describes in which direction the displacement (flag=disp) is prescribed
 *   or around which axis (flag=rot) the rotation is prescribed.
 */
class beam::fem::SupportsFromFile
{
private:
    //! Maps node ID to read in support type 
    typedef std::multimap< unsigned, Support >  MapNodeIdToSupport_;

public:
    typedef MapNodeIdToSupport_::const_iterator SupportConstIter;

public:
    //! Constructor with given file stream inp
    //!
    //! The constructor reads the provided input file stream and extracts
    //! the listed supports. The extracted supports are stored in a
    //! vector.
    SupportsFromFile( std::istream & inp );

    //! Begin of support container
    SupportConstIter begin() const { return supportContainer_.begin(); }

    //! End of support container
    SupportConstIter end() const { return supportContainer_.end(); }

    //! Print supports
    std::ostream & print( std::ostream& os ) const;

private:
    //! Container of supports
    MapNodeIdToSupport_ supportContainer_;

};

//------------------------------------------------------------------------------
/** \brief Apply Dirichlet boundary conditions, i.e. supports, on beam
 *         nodes by linearly relating the displacement DOFs of
 *         a patch of nodes
 *
 *  \tparam   LBASIC  link type
 *  \tparam   DEGREE  spline degree
 */
template< typename LBASIC, unsigned DEGREE >
class beam::fem::SupportsAtNode
{
public:
    typedef typename LBASIC::NodePtr         NodePtr;

private:
    typedef LBASIC                           LinkBasic_;

public:
    //! Main method to apply supports
    void set( std::vector< LinkBasic_ * > & links,
              const Support link,
              const std::array< NodePtr, DEGREE > & nodes )
    {
        FTL_VERIFY_DESCRIPTIVE( false, "Not implemented" );
    }
    
};

/** \brief Partial specialisation of beam support for cubic splines
 */
template< typename LBASIC >
class beam::fem::SupportsAtNode< LBASIC, 3 >
{
public:
    typedef typename LBASIC::Node            Node;
    typedef typename LBASIC::Node::VecDim    VecDim;

    static const unsigned dim = Node::dim;

private:
    typedef LBASIC                           LinkBasic_;

public:
    //! Main method to apply supports for cubic splines
    void set( std::vector< LinkBasic_ * > & links,
              const Support supp, //const Support supp,
              const std::array< const Node *, 3 > & nodes );


};

//------------------------------------------------------------------------------
// implementations
#include "Supports.ipp"

#endif
