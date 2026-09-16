// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//
//                 Computational Structural Mechanics Lab
//                         University of Cambridge
//
// This software is copyrighted and all rights are retained by the University of
// Cambridge. All use, disclosure, and/or reproduction of any part not expressly
// authorised by the University of Cambridge is prohibited. (C) 2010
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

//! @file ElementDiffusion.hpp

#ifndef elementdiffusion_h
#define elementdiffusion_h
//------------------------------------------------------------------------------
#include <Eigen/Core>
#include <corlib/eigenX.hpp>

#include <del2/fem/ElementPotential.hpp>

//------------------------------------------------------------------------------
namespace del2{
    namespace fem{
        template<typename BELEMENT> class ElementDiffusion;
        namespace eigenX = corlib::eigenX;
    }
}

//------------------------------------------------------------------------------
/** \brief Adds dynamic features to #ElementPotential for parabolic equations
 *  \detail The computation of the mass matrix and the nodal forces for a 
 *  one-step time integrator are implemented in this element.
 *  \tparam BELEMENT Type of basis element
 */
template<typename BELEMENT>
class del2::fem::ElementDiffusion : 
    public del2::fem::ElementPotential<BELEMENT>
{
private:
    typedef BELEMENT                                           BasisElement;
    typedef typename del2::fem::ElementPotential<BasisElement> ElementPotential;

public:
    typedef typename ElementPotential::VecLDim        VecLDim;
    typedef typename ElementPotential::VecNN          VecNN;
    typedef typename ElementPotential::MatDofNN       MatDofNN;

    static const unsigned numNodes = ElementPotential::numNodes;

    //--------------------------------------------------------------------------
    /** Computation of the bilinear form corresponding to the scalar mass matrix
     *  \f[
     *      M[i,j] = \int  \rho \phi_i(x) \phi_j(x) dx
     *  \f]
     *  \param[in]  xi      Quadrature coordinate
     *  \param[in]  weight  Quadrature weight
     *  \param[out] result  Container for result
     */
    void massIntegrand( const VecLDim & xi, 
                        const double & weight,
                        Eigen::MatrixXd & result ) const
    {
        // debug check
        FTL_VERIFY( (result.rows() == numNodes) and
                    (result.cols() == numNodes) );

        // get Jacobian
        const double detJ = BasisElement::jacobian( xi );

        // evaluate shape functions
        VecNN phi; BasisElement::sfun( xi, phi );
        
        result += ( phi * phi.transpose( ) ) * ( weight * detJ * rho_ );

        return;
    }

    //--------------------------------------------------------------------------
    /** Computation of inertia forces for a one-step time integrator
     *  \f[
     *        F[i] = \frac{\rho}{\Delta t} \int \phi_i(x) u(x) dx
     *  \f]
     *  \param[in] xi      Quadrature coordinate
     *  \param[in] weight  Quadrature weight
     *  \param[in] factor  Multiplier
     */
    void computeInertialForces( const VecLDim & xi, 
                                const double & weight,
                                const double stepSize )
    {
        // get nodal solutions
        MatDofNN P; 
        typename BasisElement::NodeConstIterator first = BasisElement::nodesBegin();
        typename BasisElement::NodeConstIterator last  = BasisElement::nodesEnd();
        for ( unsigned n = 0; first != last; ++first, n++ ) {
            P.col( n ) = (*first) -> getPotential( );
        }
        // evaluate shape functions
        VecNN  phi; this -> sfun( xi, phi );
        // get Jacobian
        const double detJ = this -> jacobian( xi );
        // compute nodal force contributions
        const MatDofNN forces = ( rho_ * weight * detJ / stepSize ) * 
                ( P * ( phi * phi.transpose( ) ) );
        // pass to nodes
        first = BasisElement::nodesBegin();
        for ( unsigned v = 0; first != last; ++first, v++ ) {
            (*first) -> addToForce( forces.col( v ) );
        }
        return;
    }

    //--------------------------------------------------------------------------
    //! Store the mass density of the material, \param[in] rho Mass density
    void setDensity( const double rho ) { rho_ = rho; }

protected:
    double rho_; //!< Mass density
};


#endif
