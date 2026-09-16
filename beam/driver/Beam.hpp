// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//
//                 Computational Structural Mechanics Lab
//                         University of Cambridge
//
// This software is copyrighted and all rights are retained by the University of
// Cambridge. All use, disclosure, and/or reproduction of any part not expressly
// authorised by the University of Cambridge is prohibited. (C) 2010
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

//! @file Beam.hpp

#ifndef beam_driver_beam_h
#define beam_driver_beam_h

//------------------------------------------------------------------------------
// use for nice output
#include <boost/format.hpp>
#include <boost/type_traits.hpp>
#include <boost/noncopyable.hpp>

#include <Eigen/Core>

// corlib includes
#include <corlib/Shapefun.hpp>
#include <corlib/Quadrature.hpp>
#include <corlib/Integrator.hpp>
#include <corlib/Accessor.hpp>
#include <corlib/MatrixAssembler.hpp>
#include <corlib/VectorAssembler.hpp>
#include <corlib/Distributor.hpp>
#include <corlib/NewmarkAssembler.hpp>
#include <corlib/BackwardEulerAssembler.hpp>
#include <corlib/Constraints.hpp>
#include <corlib/BodyForce.hpp>
#include <corlib/DofNumberer.hpp>
#include <corlib/PropertiesParser.hpp>
#include <corlib/VTUwriter.hpp>
#include <corlib/TensorSpline.hpp>
#include <corlib/AccumulateQuantity.hpp>

// beam::fem includes
#include <beam/fem/NodeDynamic.hpp>
#include <beam/fem/ElementPlanarDynamic.hpp>
#include <beam/fem/Mesh.hpp>
#include <beam/fem/Supports.hpp>
#include <beam/fem/misc.hpp>

// solid::material includes 
#include <solid/material/MaterialContainer.hpp>
#include <solid/material/MaterialBase.hpp>

//------------------------------------------------------------------------------
namespace beam{
    namespace driver{

        template<typename ELEMENT, typename SOLVER, unsigned NGP>
        class Beam;
        
        //! Convenience traits for the user 
        template< unsigned EDIM, unsigned DEGREE, bool DYNAMIC >
        struct Traits;

        //!\cond SKIPDOX
        template< unsigned EDIM, unsigned DEGREE >
        struct Traits<EDIM,DEGREE,false>
        {
            typedef beam::fem::NodeStatic<EDIM,EDIM>                     Node;
            typedef solid::material::MaterialBase                        Material;
            typedef corlib::TensorSpline<DEGREE,1>                       SFun;
            typedef beam::fem::ElementPlanarStatic<Node,SFun,Material>   Element;
        };
        
        template< unsigned EDIM, unsigned DEGREE >
        struct Traits<EDIM,DEGREE,true>
        {
            typedef beam::fem::NodeDynamic<EDIM,EDIM>                    Node;
            typedef solid::material::MaterialBase                        Material;
            typedef corlib::TensorSpline<DEGREE,1>                       SFun;
            typedef beam::fem::ElementPlanarDynamic<Node,SFun,Material>  Element;
        };
        //!\endcond

        //----------------------------------------------------------------------
        namespace detail_{

            //! \cond SKIPDOX
            //! Carries a flag for distinction between dynamic and static cases
            template<typename NODE>
            struct IsDynamic
            {
                typedef beam::fem::NodeDynamic<NODE::dim,NODE::dof> DynamicNode;
                static const bool result = std::is_base_of<DynamicNode,NODE>::value;
            };
            
            //! Small writer allowing for dynamic data when available
            template<typename VTUWRITER, typename NODE, bool DYNAMIC>
            struct WriteDynamicSolution;

            template<typename VTUWRITER, typename NODE>
            struct WriteDynamicSolution<VTUWRITER,NODE,false>
            {
                void operator()( VTUWRITER & vtuwriter ) { ; }
            };

            template<typename VTUWRITER, typename NODE>
            struct WriteDynamicSolution<VTUWRITER,NODE,true>
            {
                void operator()( VTUWRITER & vtuwriter ) { 
                    vtuwriter.writeNodalQuantity( corlib::accessorFun( &NODE::giveVelocities   , "Velocities"    ) );
                    vtuwriter.writeNodalQuantity( corlib::accessorFun( &NODE::giveAccelerations, "Accelerations" ) );
                }
            };
            //! \endcond

        }
    }
}

//------------------------------------------------------------------------------
/** \brief Driver object for more concise applications
 *  \details This object encapsulates the typical data construction and method
 *  calls for applications using this beam module.
 *  \tparam ELEMENT  Type of beam element (for convenience, use DriverTraits)
 *  \tparam SOLVER   Type of system solver
 *  \tparam NGP      Number of Gauss points (defaults to 2)
 */
template<typename ELEMENT, typename SOLVER, unsigned NGP = 2 >
class beam::driver::Beam
    : boost::noncopyable
{
public:
    // element and mesh
    typedef ELEMENT                                                Element;
    typedef typename Element::Node                                 Node;
    typedef solid::material::MaterialBase                          Material;
    typedef beam::fem::Mesh<Element>                               Mesh;
    typedef beam::fem::LinkBasic<Node>                             Link;

    // system solver
    typedef SOLVER                                                 Solver;

    // quadrature definition
    static const unsigned numGaussPoints = NGP;
    typedef corlib::Quadrature< corlib::LINE, numGaussPoints >     Quadrature;
    
    // Integration types
    typedef boost::function< void(Element*, 
                                  const typename Element::VecLDim &, 
                                  const double & )>                Integrand;
    typedef corlib::Integrator< Quadrature, Integrand, Element >   Integrator;

    // Matrix return object
    typedef Eigen::MatrixXd                                        Matrix;
    typedef std::function< void( const Element*,
                                   Matrix & )>                     MatrixDonator;

    // small vectors
    typedef typename Node::VecDim                                  VecDim;
    typedef typename Node::VecDof                                  VecDof;

    // body force types
    typedef std::function< VecDof( VecDim ) >                      BodyForceFun;
    typedef corlib::BodyForce<Element, BodyForceFun >              BodyForce;

    typedef std::function<bool(const double)>                      DoubleToBool;

    /// Internal enum for switch time integrator type
    enum TimeIntegrator {
        TI_UNDEFINED,        ///< default value, time integrator is undefined
        TI_STATIC,           ///< static case (no dynamics, pseudo integrator)
        TI_BACKWARDEULER,    ///< Backward Euler
        TI_NEWMARK           ///< Newmark's method
    };

public:
    //! Constructor which reads the mesh file
    Beam( std::istream & smf ) 
        : mesh_( smf ),
          bodyForceFun_( NULL ), 
          externalTractionFun_( NULL ),
          solver_( NULL ),
          timeIntegrator_( TI_UNDEFINED ), 
          newmarkBeta_( 0.25), newmarkGamma_( 0.5 ),
          rayleighM_( 0. ),    rayleighK_( 0. )
    { }

    //! Destroy dynamic memory
    virtual ~Beam() { this -> deAllocateSolver(); }

    //! Set time integrator
    void setTimeInterator( const TimeIntegrator ti,
                           const double beta  = 0.25,
                           const double gamma = 0.5 )
    { 
        timeIntegrator_ = ti;
        newmarkBeta_    = beta;
        newmarkGamma_   = gamma;
    }

    //! Pass the material pointer to the elements
    void setMaterial( Material * material ) 
    { 
        mesh_.iterateOverElements( std::bind( &Element::cacheMaterial, std::placeholders::_1, material ) );
        return;
    }

    //! Set structural damping
    void setRayleighDamping( const double facM, const double facK )
    {
        rayleighM_ = facM;
        rayleighK_ = facK;
    }

    //! Set elements thickness data (rectangular cross-section assumed)
    void setThickness( const double thickness )
    {
        mesh_.iterateOverElements( std::bind( &Element::setThickness, std::placeholders::_1, thickness ) );

        return;
    }
    
    //! Set elements cross-section data (arbitrary cross-section)
    void setCrossSection( const double area, const double secondMoment )
    {
        mesh_.iterateOverElements( std::bind( &Element::setArea,       std::placeholders::_1, area ) );
        mesh_.iterateOverElements( std::bind( &Element::setSecMomArea, std::placeholders::_1, secondMoment ) );
                
        return;
    }

    //! Read supports from file and set the links
    void prepareLinks( std::istream & supp )
    {
        beam::fem::SupportsFromFile supportsReader( supp );
        mesh_.setLinksOfSupports( supportsReader );
    }

    //! Set a functor for the body force
    void setBodyForceFun( BodyForceFun bodyForceFun )
    {
        bodyForceFun_ = bodyForceFun;
    }

    //! Set a functioning functor for the external traction functor
    void setExternalTractionFun( Integrand externalTraction )
    {
        externalTractionFun_ = externalTraction;
    }

    //! Number the dofs (nodes and links)
    std::pair<unsigned,unsigned> numberDofs()
    {
        const unsigned numDofs = mesh_.iterateOverNodes( corlib::dofNumberFun( &Node::numberDOFs, 0 ) );
        const unsigned numLags = mesh_.iterateOverLinks( corlib::dofNumberFun( &Link::numberDOFs, numDofs ) ) 
            - numDofs;
        return std::make_pair( numDofs, numLags );
    }

    //! @name Matrix computation and assembly
    //@{
    void computeElementMass();
    void computeElementStiffness();
    void assembleStiffness();
    void assembleMatricesNewmark( const double stepSize );
    void assembleMatricesBE(      const double stepSize );
    void computeAndAssembleLHS( const double stepSize = 0. );
    //@}
    
    //! @name Force computation and assembly
    //@{
    void clearForces();
    void computeForces();
    virtual void computeExternalForces( );
    void computeInternalForces( );
    void computeForcesNewmark(  const double stepSize );
    void computeForcesBE(       const double stepSize );
    void assembleForces();
    void computeAndAssembleRHS( const double stepSize = 0. );
    //@}
    
    void assembleLinks();

    //--------------------------------------------------------------------------
    //! @name Solver related interfaces
    //@{
    void allocateSolver( const unsigned size ) {
        solver_ = new Solver( size );
    }

    void deAllocateSolver() {
        delete solver_;
        solver_ = NULL;
    }

    //! Solve the global linear system of equations
    void solveSystem() { solver_ -> solveSystem(); }

    //! Return pointer to solver
    Solver * getSolver() { return solver_; }
    //@}

    //--------------------------------------------------------------------------
    //! @name Pass back and post-process the solution
    //@{
    void distributeSolution();
    //! Return norm of max increment
    double maxDisplacementIncrement() const;
    double normSolution() { return solver_ -> normSol( ); }
    void updateNewmark( const double stepSize );
    void updateBE(      const double stepSize );
    void dynamicUpdate( const double stepSize );
    void updateSolution( ); 
    void updateSolution( const double stepSize )
    {
        this -> dynamicUpdate( stepSize );
        this -> updateSolution();
    }
    //@}

    //! Write a mesh file
    void writeMeshData( const std::string & vtuFile );

    //! Convenience monitoring of a specific node
    void monitorNode( const unsigned nodeId, const double time, std::ostream & out ) const
    {
        const Node * node = mesh_.getNodePointer( nodeId );
        VecDof u = node -> giveDisplacements();
        VecDof v = node -> giveVelocities();
        VecDof a = node -> giveAccelerations();
        out << boost::format( "%1% %|15t|%2% %|30t|%3% %|45t|%4% %|60t|%5% %|75t|%6% %|90t|%7% \n" )
            % time % u[0] %u[1] % v[0] % v[1] % a[0] % a[1];
        out << std::flush;
    }

    //! Return value of a link
    double linkValue( const unsigned linkId ) const 
    {
        const typename Mesh::Link * link = mesh_.getLinkPointer( linkId );
        return (link -> giveMultiplier())[0];
    }

    //! Allow to access the mesh
    Mesh & accessMesh() { return mesh_; }

    //! Short cut to apply functor to nodes
    template<typename OP>
    OP iterateOverNodes( OP op ) { return mesh_.iterateOverNodes( op ); }

    //! Short cut to apply functor to elements
    template<typename OP>
    OP iterateOverElements( OP op ) { return mesh_.iterateOverElements( op ); }

    //! do one time step (without update)
    bool iterate( const double timeStepSize,
                  const unsigned maxIterations,
                  DoubleToBool convCheckRhs,
                  DoubleToBool convCheckSolution );
    
private:
    Mesh         mesh_;                 //!< Mesh object
    BodyForceFun bodyForceFun_;         //!< Functor for the body force
    Integrand    externalTractionFun_;  //!< Functor for load by traction
    Solver *     solver_;               //!< Pointer to system solver
    Quadrature   quadrature_;           //!< Quadrature object

    //! @name Time integration parameters
    //@{ 
    TimeIntegrator  timeIntegrator_;  //!< The time integrator
    double          newmarkBeta_;     //!< Newmark's beta parameter
    double          newmarkGamma_;    //!< Newmark's gamma parameter
    //@}

    //! @name Structural damping
    //@{
    double  rayleighM_; //!< mass      matrix contribution to damping
    double  rayleighK_; //!< stiffness matrix contribution to damping
    //@}
};
//------------------------------------------------------------------------------
#include "Beam.ipp"

#endif

