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

//! @author Burkhard Bornemann, Kosala Bandara
//! @date   2012

#ifndef gshell_fem_nodeshell_h
#define gshell_fem_nodeshell_h

#include <iostream>
#include <cassert>
#include <array>

#include <Eigen/Core>
#include <corlib/eigenX.hpp>

//------------------------------------------------------------------------------
namespace gshell{
    namespace fem{

        class NodeShell;

        namespace eigenX = corlib::eigenX;
    }
}

//==============================================================================
/** \brief Base class for nodes
 *  \details Stores an ID-number and the coordinates
 */
class gshell::fem::NodeShell
{
public:
    //! The spatial dimension of this node
    static const unsigned dim = 3;

    //! Coordinate array
    typedef eigenX::VectorSd<dim>      VecDim;

    typedef std::map< NodeShell *, double >           MapNPtrDouble;

protected:
    //! Spatial dimension of shell
    static const unsigned localDim_ = 2;
    //! Container for subdivision weights to compute tangents
    typedef std::array< MapNPtrDouble, localDim_ >  MapNPtrDoubleArray_;

public:
    //--------------------------------------------------------------------------
    //! \brief Empty constructor
    //! \details Invalidates the stored data
    NodeShell()
        : id_( std::numeric_limits<unsigned>::max() ),
          coord_( Eigen::VectorXd::Ones( dim ) * std::numeric_limits<double>::max() )
    { }
    
    //--------------------------------------------------------------------------
    //! @name Coordinates
    //@{

    //! Set coordinates \details
    //!
    //! \param[in] c  Vector containing the coordinates  */
    void setCoordinates( const VecDim & c ) { coord_ = c; }

    //! Give the coordinates \details
    //!
    //! \return The coordinates of this node */
    VecDim giveCoordinates() const { return coord_; };

    //! Insert coordinates to given insert iterator
    //! deactivated for now 
    //template< typename IT >
    //IT passCoordinates( IT backIter ) const;
    //@}

    //--------------------------------------------------------------------------
    //! @name ID
    //@{
    
    //! Give node's ID \details 
    //!
    //! \return       The ID of this node */
    unsigned giveId() const { return id_; }

    //! Set the ID of this node \details
    //!
    //! \param[in] i  New ID of this node */
    void setId( const unsigned & i ) { id_ = i; return; }
    //@}
    
    //--------------------------------------------------------------------------
    //! @name I/O
    //@{

    //! Write id to a stream
    std::ostream & writeId( std::ostream & out ) const;
    //! Write coordinates to a stream
    std::ostream & writeCoordinates( std::ostream & out ) const;
    //! Read coordinates from a stream 
    std::istream & readSelf( std::istream & inp );
    //@}

    //--------------------------------------------------------------------------
    //! @name Actions on tangents
    //@{

    //! Clear coefficients of tangents
    void clearTangentCoefficients()
    {
        std::for_each( tangCoeffs_.begin(), tangCoeffs_.end(),
                       std::bind( &MapNPtrDouble::clear, std::placeholders::_1 ) );
        return;
    }

    //! Return number of limit coefficients
    //!
    //! \param[in]  dir     Tangential direction (0,1)
    //! \return             Number of limit coefficients
    unsigned getNumberTangentCoefficients( const unsigned dir ) const
    { 
        return tangCoeffs_[ dir ].size();
    }

    //! Set additional coefficients of tangent at #this w.r.t. #node
    //!
    //! \param[in]  node    Node of which the coefficient is added
    //! \param[in]  dir     Tangential direction (0,1)
    //! \param[in]  coeff   Coefficients of #dir-th tangent w.r.t. #node
    void addTangentCoefficient( NodeShell * node, const unsigned dir,
                                const double coeff );

    //! Get coefficient to tangent at #this w.r.t. #node
    //!
    //! \param[in]  node   Node to which the coefficient is searched
    //! \param[in]  dir    Tangential direction (0,1)
    //! \return            Coefficient, if node not found deliver zero
    double getTangentCoefficient( NodeShell * node, const unsigned dir ) const;

    //@}

    //--------------------------------------------------------------------------
    //! @name Actions for limit surface
    //@{

    //! Clear limit surface coefficients
    void clearLimitCoefficients() { limitCoeff_.clear(); }

    //! Return number of limit coefficients
    unsigned getNumberLimitCoefficients() const { return limitCoeff_.size(); }

    //! Set additional coefficient of limit surface at #this w.r.t. #node
    //!
    //! \param[in]  node    Node of which the coefficient is added
    //! \param[in]  coeff   The limit surface coefficient
    void addLimitCoefficient( NodeShell * node, const double coeff );

    //! Get coefficient of limit surface at #this w.r.t. #node
    //!
    //! \param[in]  node    Node of which the coefficient is searched
    //! \return     coeff   The limit surface coefficient
    double getLimitCoefficient( NodeShell * node ) const;

    //@}


protected:
    //--------------------------------------------------------------------------
    unsigned     id_;              //!< identifiying number
    VecDim       coord_;           //!< coordinates

    //! Subdivision coefficients to compute unique tangents
    MapNPtrDoubleArray_  tangCoeffs_;

    //! Subdivision coefficients to compute limit surface
    MapNPtrDouble        limitCoeff_;

};

#include "NodeShell.ipp"

#endif
