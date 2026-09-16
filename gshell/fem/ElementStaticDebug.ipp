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

// @author Burkhard Bornemann, Kosala Bandara
// @date   2010

// //------------------------------------------------------------------------------
// //! Write shape functions and 1st derivatives to Gnuplot files
// //!
// //! Compute the shape functions \f$\phi^K\f$ at a series of points,
// //! i.e. \f$\phi^K(\xi^\alpha_\mathrm{centre}+\epsilon\eta^\alpha)\f$,
// //! and their local gradient \f$\phi_{,\alpha}\f$ at the centre of the interval.
// //! The results are written to data files and Gnuplot driver files
// //! for each shape function \f$K\f$.
// //!
// //! Directional tangent at centre
// //!\f[
// //!     t^K(\xi^\alpha_\mathrm{centre},\epsilon,\eta^\beta)
// //!     = \phi^K(\xi^\alpha_\mathrm{centre})
// //!     + \left.\frac{\partial \phi^K}{\partial \xi^\beta}\right|_{\xi^\alpha_\mathrm{centre}}
// //!       \epsilon\eta^\beta
// //!\f]
// //! in which \f$\epsilon\in[-h/2,h/2]\f$ and \f$\eta^\beta\f$ is a unit-vector.
// //!
// //! <h4>NOTE:</h4>
// //! This method is only intended for debug purposes.
// //!
// //! \param[in]  baseName         Folder name  to which data is written,
// //! \param[in]  centre           Centre point \f$\xi^\alpha_\mathrm{centre}\f$
// //! \param[in]  size             Interval size \f$h\f$
// //! \param[in]  direction        Direction unit vector \f$\eta^\beta\f$
// //!                                (It is scaled to become unit vector for convenience.)
// //! \param[in]  numSteps         Number of evaluation steps in which interval
// //!                                is divided
// //! \param[in]  scaleToOne       Scale shape function and tangent to be one at centre
// template< typename BELEMENT, typename MAT >
// void gshell::fem::ElementStatic<BELEMENT,MAT>::gnuPlotShapeFunctions(
//     const std::string & baseName,
//     const VecLDim & centre,
//     const double & size,
//     const VecLDim & direction,
//     const unsigned & numSteps,
//     const bool & scaleToOne
//     )
// {
//     shapeFun_->gnuPlotShapeFunctions( baseName, centre, size, direction,
//                                       numSteps, scaleToOne );
//     return;
// }

// //------------------------------------------------------------------------------
// //! Write 1st and 2nd derivatives of shape functions to Gnuplot files
// //!
// //! Compute the 1st derivatives of the shape functions \f$\phi^K_{,\xi^\gamma}(\xi^\alpha)\f$
// //! at a series of points and their local 2nd derivatives \f$\phi^K_{,\beta,\gamma}(\xi^\alpha)\f$
// //! at the centre of the interval.
// //! The results are written to data files and Gnuplot driver files
// //! for each shape function \f$K\f$.
// //!
// //! Directional tangent at centre
// //!\f[
// //!     t^K_\gamma(\xi^\alpha_\mathrm{centre},\epsilon,\eta^\beta)
// //!     = \phi^K_{,\gamma}(\xi^\alpha_\mathrm{centre})
// //!     + \left.\frac{\partial \phi^K_{,\gamma}}{\partial \xi^\beta}\right|_{\xi^\alpha_\mathrm{centre}}
// //!       \epsilon\eta^\beta
// //!\f]
// //! in which \f$\epsilon\in[-h/2,h/2]\f$ and \f$\eta^\beta\f$ is a unit-vector.
// //!
// //! <h4>NOTE:</h4>
// //! This method is only intended for debug purposes.
// //!
// //! \param[in]  baseName         Folder name  to which data is written,
// //! \param[in]  centre           Centre point \f$\xi^\alpha_\mathrm{centre}\f$
// //! \param[in]  size             Interval size \f$h\f$
// //! \param[in]  direction        Direction (unit) vector \f$\eta^\beta\f$
// //!                                (is scaled to become unit vector for convenience)
// //! \param[in]  numSteps         Number of evaluation steps in which interval
// //!                                is divided
// //! \param[in]  scaleToOne       Scale shape function and tangent to be one at centre
// template< typename BELEMENT, typename MAT >
// void gshell::fem::ElementStatic<BELEMENT,MAT>::gnuPlotShapeDeriv1(
//     const std::string & baseName,
//     const VecLDim & centre,
//     const double & size,
//     const VecLDim & direction,
//     const unsigned & numSteps,
//     const bool & scaleToOne
//     )
// {
//     shapeFun_->gnuPlotShapeDeriv1( baseName, centre, size, direction,
//                                    numSteps, scaleToOne );
//     return;
// }

//------------------------------------------------------------------------------
//! Write mid-surface vectors and derivatives to Gnuplot files
//!
//! Compute mid-surface vectors (like bases, director, etc) and their
//! 1st derivatives w.r.t. nodal DOFs. The vectors are computed for a range
//! of values and a tangent is computed at the middle of the range as well.
//! The vectors can be visualised with Gnuplot.
//!
//! <h4>NOTE:</h4>
//! This method is only intended for debug purposes.
//!
//! \param[in]  vectorName       Name of vector to investigate
//! \param[in]  a                Component index \f$A\f$
//! \param[in]  be               Local gradient index \f$\beta\f$
//!                                --- only for certain quantities
//! \param[in]  fileBaseName     Basename of file to which data is written
//! \param[in]  mid              Evaluation middle
//! \param[in]  width            Evaluation range
//! \param[in]  steps            Number of evaluation steps in range
//! \param[in]  dofLow           Lowest DOF index to test
//! \param[in]  dofUp            Highest DOF + 1 index to test
template< typename BELEMENT, typename MAT >
void gshell::fem::ElementStatic<BELEMENT,MAT>::gnuPlotVector(
    const std::string vectorName,
    const unsigned a,
    const unsigned be,
    const std::string fileBaseName,
    const double mid,
    const double width,
    const unsigned steps,
    const unsigned dofLow,
    const unsigned dofUp
    )
{
    // constants
    const unsigned numFunctions = this->numFunctions( );
    const double halfWidth = 0.5 * width;
    double delta[dim][dim];  la::identity<dim>( delta );
    double zeroDim[dim];     la::zero<dim>( zeroDim );

    // checks
    FTL_VERIFY( a < dim );
    //FTL_VERIFY( be < localDim );
    FTL_VERIFY( be < dim );
    FTL_VERIFY_DESCRIPTIVE( steps%2==0, "Only even number of steps allowed\n" );
    FTL_VERIFY( dofLow < dof );
    FTL_VERIFY( ( dofUp > 0 ) and ( dofUp <= dof ) );
    FTL_VERIFY( dofLow < dofUp );

    // open files
    std::stringstream gnuplotFileName;
    gnuplotFileName << fileBaseName << ".gplt";
    std::ofstream gnuplotFile( gnuplotFileName.str( ).c_str( ) );
    FTL_VERIFY( gnuplotFile.is_open( ) );

    // write Gnuplot driver
    gnuplotFile << "# Gnuplot script" << std::endl;
    gnuplotFile << "unset key" << std::endl;
    gnuplotFile << "set palette rgbformulae 22,13,-31" << std::endl;
    gnuplotFile << "plot \\" << std::endl;
    const unsigned maxLines = (dofUp-dofLow)*numFunctions;
    unsigned lines = 0;
    
    // loop Dofs
    for ( unsigned l=0; l<numFunctions; ++l ) {
        for ( unsigned dd=0; dd<dof; ++dd ) {
            const unsigned d = dd % dim;

            std::stringstream vectorFileName;
            vectorFileName << fileBaseName << "." << l << "." << dd << ".vector";
            std::ofstream vectorFile( vectorFileName.str( ).c_str( ) );
            FTL_VERIFY( vectorFile.is_open( ) );

            std::stringstream tangentFileName;
            tangentFileName << fileBaseName << "." << l << "." << dd << ".tangent";
            std::ofstream tangentFile( tangentFileName.str( ).c_str( ) );
            FTL_VERIFY( tangentFile.is_open( ) );

            const std::string ffn = vectorFileName.str( ).substr(
                vectorFileName.str( ).find("/") + 1,
                vectorFileName.str( ).length( ) );
            const std::string tfn = tangentFileName.str( ).substr(
                tangentFileName.str( ).find("/") + 1,
                tangentFileName.str( ).length( ) );
            const double frac = static_cast< double >( lines ) / 
                                static_cast< double >( maxLines );

            if ( ( dd >= dofLow ) and ( dd < dofUp ) ) {
                gnuplotFile << " \"" << tfn << "\""
                            << " u 1:" << 2
                            << " t \"tangent." << l << "." << dd << "\""
                            << " w l lw 2 lt 2"
                            << " lc palette frac " << frac;
                gnuplotFile << ",";
                gnuplotFile << " \"" << ffn << "\""
                            << " u 1:" << 2
                            << " t \"" << vectorName << "." << l << "." << dd << "\""
                            << " w lp lw 1"
                            << " lc palette frac " << frac;
                lines += 1;
                if ( lines != maxLines ) gnuplotFile << ", \\";
                gnuplotFile << std::endl;
            }

            // blank displacements
            for ( unsigned i = 0; i < numFunctions; ++i ) {
                this->BasisElement::giveSupportNodePtr( i )->clearDisplacements();
                this->BasisElement::giveSupportNodePtr( i )->clearIncrement();
            }

            double low = mid - halfWidth;
            for ( unsigned s=0; s<=steps; ++s ) {
            
                // prescribe displacement
                VecDof disp = eigenX::VectorSd<dof>::Zero( );
                disp[ dd ] = low;
                this->BasisElement::giveSupportNodePtr( l )->setDisplacements( disp );

                // Evaluation point
                const VecLDim xi = corlib::ShapeTraits< myShape >::centroid( );                

                // Collect shape functions and their derivatives
                VecNF phi;            this->BasisElement::sfun(  xi, phi );
                MatLDimNF dPhiDXi;    this->BasisElement::sfunGrad( xi, dPhiDXi );
                MatSDimNF ddPhiDDXi;  this->BasisElement::sfunHess( xi, ddPhiDDXi );

                // local shears
                MatDimNF uNF, wNF;   this->nodalDisplacements( uNF, wNF );
                double wLocal[dim];       SurfaceGeometer::interpolate( phi, wNF, wLocal );
                double zLocal[localDim][dim];  SurfaceGeometer::interpolateDeriv1( dPhiDXi, wNF, zLocal );

                // get reference mid-vectors
                double aRef[dim][dim];        double bRef[dim][localDim][dim];
                double wRef[dim];             double zRef[localDim][dim];
                double pRef[dim];             double pRefLGrad[localDim][dim];
                double qRef[dim];             double qRefLGrad[localDim][dim];
                double dRef[dim];             double cRef[localDim][dim];
                this->vectors_( REFERENCE, xi, phi, dPhiDXi, ddPhiDDXi,
                                aRef, bRef, wRef, zRef,
                                pRef, pRefLGrad, qRef, qRefLGrad, dRef, cRef );
                
                // get current mid-vectors
                double aCur[dim][dim];        double bCur[dim][localDim][dim];
                double wCur[dim];             double zCur[localDim][dim];
                double pCur[dim];             double pLGrad[localDim][dim];
                double qCur[dim];             double qLGrad[localDim][dim];
                double dCur[dim];             double cCur[localDim][dim];
                this->vectors_( CURRENT, xi, phi, dPhiDXi, ddPhiDDXi,
                                aCur, bCur, wCur, zCur,
                                pCur, pLGrad, qCur, qLGrad, dCur, cCur );

                // 1st derivatives of vectors
                double pDU_dl[dim];           double pLGradDU_dl[localDim][dim];
                double qDU_dl[dim];           double qLGradDU_dl[localDim][dim];
                double qDW_dl[dim];           double qLGradDW_dl[localDim][dim];
                double dDU_dl[dim];           double cDU_dl[localDim][dim];
                double dDW_dl[dim];           double cDW_dl[localDim][dim];
                this->vectorsGradCur_( d, l, xi, phi, dPhiDXi, ddPhiDDXi,
                                       aRef, aCur, bCur, wCur, zCur, 
                                       pCur, pLGrad, qCur, qLGrad, dCur, cCur,
                                       pDU_dl, pLGradDU_dl, qDU_dl, qDW_dl, qLGradDU_dl, qLGradDW_dl,
                                       dDU_dl, dDW_dl, cDU_dl, cDW_dl );

                // Green-Lagrange strain resultants
                double alphaCur[dim][dim], betaCur[dim][dim];
                double alphaDU_dl[dim][dim], alphaDW_dl[dim][dim];
                double betaDU_dl[dim][dim], betaDW_dl[dim][dim];
                if ( vectorName == "alpha" or vectorName == "beta" or
                     vectorName == "n" or vectorName == "m" ) {
                    this->strainRes_( aRef, dRef, cRef, aCur, dCur, cCur,
                                      alphaCur, betaCur );
                    this->strainResGradCur_( d, l, 
                                             aCur, dCur, dDU_dl, dDW_dl, cCur, cDU_dl, cDW_dl, dPhiDXi,
                                             alphaDU_dl, alphaDW_dl, betaDU_dl, betaDW_dl );
                }

                // 2PK stress
                double nCur[dim][dim], mCur[dim][dim];
                double hRef[dim][dim][dim][dim];
                if ( vectorName == "n" or vectorName == "m" ) {
                    this->stressResVK_( aRef, alphaCur, betaCur, nCur, mCur );
                    this->elasticityTensorVK_( aRef, hRef );
                }

                // initialise vectors;
                double v[dim];       la::zero<dim>( v );
                double vDU_dl[dim];  la::zero<dim>( vDU_dl );
                double vDW_dl[dim];  la::zero<dim>( vDW_dl );

                // set choice accordings to #vectorName
                if      ( vectorName == "a0" ) {
                    const unsigned al = 0;
                    la::assign<dim>( v, aCur[ al ] );
                    la::update<dim>( vDU_dl, dPhiDXi( al, l ), delta[ d ] );
                }
                else if ( vectorName == "a1" ) {
                    const unsigned al = 1;
                    la::assign<dim>( v, aCur[ al ] );
                    la::update<dim>( vDU_dl, dPhiDXi( al, l ), delta[ d ] );
                }
//                else if ( vectorName == "a2" ) {
//                    la::assign<dim>( v, aCur[ 2 ] );
//                    la::assign<dim>( vDU_dl, a2DU_dl );
//                }
                else if ( vectorName == "b0" ) {
                    const unsigned al = 0;
                    la::assign<dim>( v, bCur[ al ][ be ] );
                    la::update<dim>( vDU_dl, ddPhiDDXi( voigtForward[al][be], l ), delta[ d ] );
                }
                else if ( vectorName == "b1" ) {
                    const unsigned al = 1;
                    la::assign<dim>( v, bCur[ al ][ be ] );
                    la::update<dim>( vDU_dl, ddPhiDDXi( voigtForward[al][be], l ), delta[ d ] );
                }
//                else if ( vectorName == "b2" ) {
//                    la::assign<dim>( v, bCur[ 2 ][ be ] );
//                    la::assign<dim>( vDU_dl, b2DU_dl[ be ] );
//                }
                else if ( vectorName == "w" ) {
                    la::assign<dim>( v, wCur );
                    vDU_dl[ d ] = wLocal[ 0 ] * dPhiDXi( 0, l ) +
                        wLocal[ 1 ] * dPhiDXi( 1, l );
                    if ( d < localDim )
                        la::update<dim>( vDW_dl, phi( l ), aCur[ d ] );
                    else
                        la::assign<dim>( vDW_dl, zeroDim );
                }
//                else if ( vectorName == "z" ) {
//                    la::assign<dim>( v, zCur[ be ] );
//                    la::update<dim>( vDW_dl, dPhiDXi( be, l ), delta[ d ] );
//                }
                else if ( vectorName == "p" ) {
                    la::assign<dim>( v, pCur );
                    la::assign<dim>( vDU_dl, pDU_dl );
                    la::assign<dim>( vDW_dl, zeroDim );
                }
                else if ( vectorName == "pLGrad" ) {
                    la::assign<dim>( v, pLGrad[ be ] );
                    la::assign<dim>( vDU_dl, pLGradDU_dl[ be ] );
                    la::assign<dim>( vDW_dl, zeroDim );
                }
                else if ( vectorName == "q" ) {
                    la::assign<dim>( v, qCur );
                    la::assign<dim>( vDU_dl, qDU_dl );
                    la::assign<dim>( vDW_dl, qDW_dl );
                }
                else if ( vectorName == "qLGrad" ) {
                    la::assign<dim>( v, qLGrad[ be ] );
                    la::assign<dim>( vDU_dl, qLGradDU_dl[ be ] );
                    la::assign<dim>( vDW_dl, qLGradDW_dl[ be ] );
                }
                else if ( vectorName == "d" ) {
                    v  = dCur;
                    la::assign<dim>( vDU_dl, dDU_dl );
                    la::assign<dim>( vDW_dl, dDW_dl );
                }
                else if ( vectorName == "c" ) {
                    la::assign<dim>( v, cCur[ be ] );
                    la::assign<dim>( vDU_dl, cDU_dl[ be ] );
                    la::assign<dim>( vDW_dl, cDW_dl[ be ] );
                }
                else if ( vectorName == "alpha" ) {
                    la::assign<dim>( v, alphaCur[ be ] );
                    la::assign<dim>( vDU_dl, alphaDU_dl[ be ] );
                    la::assign<dim>( vDW_dl, alphaDW_dl[ be ] );
                }
                else if ( vectorName == "beta" ) {
                    la::assign<dim>( v, betaCur[ be ] );
                    la::assign<dim>( vDU_dl, betaDU_dl[ be ] );
                    la::assign<dim>( vDW_dl, betaDW_dl[ be ] );
                }
                else if ( vectorName == "n" ) {
                    const double t = thickness_;
                    la::assign<dim>( v, nCur[ be ] );
                    double nDU_dl[dim][dim];
                    for ( unsigned i=0; i<dim; ++i ) {
                        for ( unsigned j=0; j<dim; ++j ) {
                            la::innerProduct<dim,dim>( nDU_dl[ i ][ j ], hRef[ i ][ j ], alphaDU_dl );
                            nDU_dl[ i ][ j ] *= t;
                        }
                    }
                    la::assign<dim>( vDU_dl, nDU_dl[ be ] );
                    double nDW_dl[dim][dim];
                    for ( unsigned i=0; i<dim; ++i ) {
                        for ( unsigned j=0; j<dim; ++j ) {
                            la::innerProduct<dim,dim>( nDW_dl[ i ][ j ], hRef[ i ][ j ], alphaDW_dl );
                            nDW_dl[ i ][ j ] *= t;
                        }
                    }
                }
                else if ( vectorName == "m" ) {
                    const double t = thickness_;
                    const double s = t*t*t/12.;
                    la::assign<dim>( v, mCur[ be ] );
                    double mDU_dl[dim][dim];
                    for ( unsigned i=0; i<dim; ++i ) {
                        for ( unsigned j=0; j<dim; ++j ) {
                            la::innerProduct<dim,dim>( mDU_dl[ i ][ j ], hRef[ i ][ j ], betaDU_dl );
                            mDU_dl[ i ][ j ] *= s;
                        }
                    }
                    la::assign<dim>( vDU_dl, mDU_dl[ be ] );
                    double mDW_dl[dim][dim];
                    for ( unsigned i=0; i<dim; ++i ) {
                        for ( unsigned j=0; j<dim; ++j ) {
                            la::innerProduct<dim,dim>( mDW_dl[ i ][ j ], hRef[ i ][ j ], betaDW_dl );
                            mDW_dl[ i ][ j ] *= s;
                        }
                    }
                    la::assign<dim>( vDW_dl, mDW_dl[ be ] );
                }
                else {
                    FTL_VERIFY_DESCRIPTIVE( false, "Vector %s not found\n", vectorName.c_str( ) );
                }
         
                // write to vector to file
                vectorFile << low;
                vectorFile << " " << v[ a ];
                vectorFile << std::endl;
                
                // tangent
                if ( s == steps/2 ) {

                    // write tangent
                    tangentFile << mid-halfWidth;
                    if ( dd < dim )
                        tangentFile << " " <<  v[ a ] - vDU_dl[ a ] * halfWidth;
                    else
                        tangentFile << " " <<  v[ a ] - vDW_dl[ a ] * halfWidth;
                    tangentFile << std::endl;
                    tangentFile << mid
                                << " " <<  v[ a ]
                                << std::endl;
                    tangentFile << mid+halfWidth;
                    if ( dd < dim )
                        tangentFile << " " <<  v[ a ] + vDU_dl[ a ] * halfWidth;
                    else
                        tangentFile << " " <<  v[ a ] + vDW_dl[ a ] * halfWidth;
                    tangentFile << std::endl;
                }

                low += width / static_cast< double >( steps );

            }

            // close files
            vectorFile.close( );
            tangentFile.close( );
        }
    }
    gnuplotFile << "pause -1 \"" << gnuplotFileName.str( )
                << " --- Hit return to continue\"" << std::endl;
    gnuplotFile.close( );

    return;
}

//------------------------------------------------------------------------------
//! Write 1st and 2nd derivatives of mid-surface vectors to Gnuplot file
//!
//! Compute mid-surface vectors (like bases, director, etc) and their
//! 1st or 2nd derivatives w.r.t. nodal DOFs. The result of the derivatives
//! is compared to finite difference counterpart
//!
//! <h4>NOTE:</h4>
//! This method is only intended for debug purposes.
//!
//! \param[out]    out              Output stream
//! \param[in]     vectorName       Name of vector to investigate
//! \param[in]     a                Component index \f$A\f$
//! \param[in]     be               Local gradient index \f$\beta\f$
//!                                      --- only for certain quantities
//! \param[in]     mid              Evaluation middle
//! \param[in]     width            Evaluation range
//! \param[in]     secondDeriv      Type of 2nd derivative:
//!                                      0=do not test,
//!                                      1=w.r.t. nodal displacements,
//!                                      2=w.r.t. nodal shears
template< typename BELEMENT, typename MAT >
void gshell::fem::ElementStatic<BELEMENT,MAT>::finiteDifferenceVector(
    std::ostream & out,
    const std::string vectorName,
    const unsigned a,
    const unsigned be,
    const double mid,
    const double width,
    const unsigned secondDeriv
    )
{
    // constants
    const unsigned steps = 2;
    const unsigned numFunctions = this->numFunctions( );
    const double halfWidth = 0.5 * width;
    double delta[dim][dim];             la::identity<dim>( delta );
    double zeroDim[dim];                la::zero<dim>( zeroDim );
    double zeroLDimDim[localDim][dim];  la::zero<localDim,dim>( zeroLDimDim );

    // checks
    FTL_VERIFY( a < dim );
    FTL_VERIFY( be < dim );//FTL_VERIFY( be < localDim );
    FTL_VERIFY( secondDeriv < 3 );

    // loop Dofs
    const unsigned kMax = secondDeriv ? numFunctions : 1;
    const unsigned ccMax = secondDeriv ? dof : 1;

    // loop DOFs
    for ( unsigned k=0; k<kMax; ++k ) {
        for ( unsigned cc=0; cc<ccMax; ++cc ) {
            const unsigned c = cc % dim;

            // loop DOFs
            for ( unsigned l=0; l<numFunctions; ++l ) {
                for ( unsigned dd=0; dd<dof; ++dd ) {
                    const unsigned d = dd % dim;

                    // blank displacements
                    for ( unsigned i = 0; i < numFunctions; ++i ) {
                        this->BasisElement::giveSupportNodePtr( i )->clearDisplacements();
                        this->BasisElement::giveSupportNodePtr( i )->clearIncrement();
                    }

                    double tangent = 0.0;
                    double tangentFD = 0.0;

                    double low = mid;
                    double sig = -1.0;
                    for ( unsigned s=0; s<steps; ++s ) {
            
                        // prescribe displacement
                        VecDof disp = eigenX::VectorSd<dof>::Zero( );
                        disp[ dd ] = low;
                        this->BasisElement::giveSupportNodePtr( l )->setDisplacements( disp );

                        // Evaluation point
                        const VecLDim xi = corlib::ShapeTraits< myShape >::centroid( );

                        // Collect shape functions and their derivatives
                        VecNF phi;            this->BasisElement::sfun(  xi, phi );
                        MatLDimNF dPhiDXi;    this->BasisElement::sfunGrad( xi, dPhiDXi );
                        MatSDimNF ddPhiDDXi;  this->BasisElement::sfunHess( xi, ddPhiDDXi );

                        // local shears
                        MatDimNF uNF, wNF;             this->nodalDisplacements( uNF, wNF );
                        double wLocal[dim];            SurfaceGeometer::interpolate( phi, wNF, wLocal );
                        double zLocal[localDim][dim];  SurfaceGeometer::interpolateDeriv1( dPhiDXi, wNF, zLocal );

                        // get reference mid-vectors
                        double aRef[dim][dim];        double bRef[dim][localDim][dim];
                        double wRef[dim];             double zRef[localDim][dim];
                        double pRef[dim];             double pRefLGrad[localDim][dim];
                        double qRef[dim];             double qRefLGrad[localDim][dim];
                        double dRef[dim];             double cRef[localDim][dim];
                        this->vectors_( REFERENCE, xi, phi, dPhiDXi, ddPhiDDXi,
                                        aRef, bRef, wRef, zRef,
                                        pRef, pRefLGrad, qRef, qRefLGrad, dRef, cRef );

                        // get current mid-vectors
                        double aCur[dim][dim];        double bCur[dim][localDim][dim];
                        double wCur[dim];             double zCur[localDim][dim];
                        double pCur[dim];             double pLGrad[localDim][dim];
                        double qCur[dim];             double qLGrad[localDim][dim];
                        double dCur[dim];             double cCur[localDim][dim];
                        this->vectors_( CURRENT, xi, phi, dPhiDXi, ddPhiDDXi,
                                        aCur, bCur, wCur, zCur,
                                        pCur, pLGrad, qCur, qLGrad, dCur, cCur );
                        double pLen;  la::norm2<dim>( pLen, pCur );
                        double qLen;  la::norm2<dim>( qLen, qCur );
                        const double (&a0Cur)[dim] = aCur[ 0 ];
                        const double (&a1Cur)[dim] = aCur[ 1 ];
                        const double (&a2Cur)[dim] = aCur[ 2 ];
                        double b0Cur[localDim][dim], b1Cur[localDim][dim], b2Cur[localDim][dim];
                        for ( unsigned ga=0; ga<localDim; ++ga ) {
                            la::assign<dim>( b0Cur[ ga ], bCur[ 0 ][ ga ] );
                            la::assign<dim>( b1Cur[ ga ], bCur[ 1 ][ ga ] );
                            la::assign<dim>( b2Cur[ ga ], bCur[ 2 ][ ga ] );/*(
                                ublas::column( pLGrad, ga ) -
                                a2Cur * ublas::inner_prod( a2Cur, ublas::column( pLGrad, ga ) )
                                ) / pLen;*/
                        }
                        double sCur[localDim][dim];
                        la::assign<localDim,dim>( sCur, qLGrad );
                        la::divide<localDim,dim>( sCur, qLen );

                        double wDU_dl[dim];  la::zero<dim>( wDU_dl );
                        for ( unsigned al=0; al<localDim; ++al )
                            wDU_dl[ d ] += wLocal[ al ] * dPhiDXi( al, l );
                        double wDW_dl[dim];  la::zero<dim>( wDW_dl );
                        if ( d < localDim )  la::update<dim>( wDW_dl, phi( l ), aCur[ d ] );
                        double zDU_dl[localDim][dim];  la::zero<localDim,dim>( zDU_dl );
                        for ( unsigned ga=0; ga<localDim; ++ga )
                            for ( unsigned al=0; al<localDim; ++al )
                                zDU_dl[ ga ][ d ] += zLocal[ al ][ ga ] * dPhiDXi( al, l ) +
                                    wLocal[ al ] * ddPhiDDXi( voigtForward[al][ga], l );
                        double zDW_dl[localDim][dim];
                        for ( unsigned ga=0; ga<localDim; ++ga ) {
                            la::zero<dim>( zDW_dl[ ga ] );
                            if ( d < localDim ) {
                                la::update<dim>( zDW_dl[ ga ], dPhiDXi( ga, l ), aCur[ d ] );
                                la::update<dim>( zDW_dl[ ga ], phi( l ), bCur[ d ][ ga ] );
                            }
                        }

                        double pDU_dl[dim];                double pLGradDU_dl[localDim][dim];
                        double qDU_dl[dim];                double qLGradDU_dl[localDim][dim];
                        double qDW_dl[dim];                double qLGradDW_dl[localDim][dim];
                        double dDU_dl[dim];                double cDU_dl[localDim][dim];
                        double dDW_dl[dim];                double cDW_dl[localDim][dim];
                        this->vectorsGradCur_( d, l, xi, phi, dPhiDXi, ddPhiDDXi,
                                               aRef, aCur, bCur, wCur, zCur, 
                                               pCur, pLGrad, qCur, qLGrad, dCur, cCur,
                                               pDU_dl, pLGradDU_dl, qDU_dl, qDW_dl, qLGradDU_dl, qLGradDW_dl,
                                               dDU_dl, dDW_dl, cDU_dl, cDW_dl );
                        double a2DU_dl[dim];          double b2DU_dl[localDim][dim];
                        SurfaceGeometer::unitVectorGradDof( pLen, pDU_dl, a2Cur, a2DU_dl );
//                        SurfaceGeometer::unitVectorGradDof2( pCur, pDU_dl, a2DU_dl );
                        double sLGrad[localDim][dim];
                        la::assign<localDim,dim>( sLGrad, pLGrad );
                        la::divide<localDim,dim>( sLGrad, pLen );
                        SurfaceGeometer::unitVectorDeriv1GradDof( pLen, pDU_dl, pLGradDU_dl,
                                                                  a2Cur, a2DU_dl, sLGrad,
                                                                  b2DU_dl );

                        // Green-Lagrange strain resultants
                        double alphaCur[dim][dim], betaCur[dim][dim];
                        double alphaDU_dl[dim][dim], alphaDW_dl[dim][dim];
                        double betaDU_dl[dim][dim], betaDW_dl[dim][dim];
                        if ( vectorName == "alpha" or vectorName == "beta" or
                             vectorName == "n" or vectorName == "m" ) {
                            this->strainRes_( aRef, dRef, cRef, aCur, dCur, cCur,
                                              alphaCur, betaCur );
                            this->strainResGradCur_( d, l, 
                                                     aCur, dCur, dDU_dl, dDW_dl, cCur, cDU_dl, cDW_dl, dPhiDXi,
                                                     alphaDU_dl, alphaDW_dl, betaDU_dl, betaDW_dl );
                        }

                        // stuff for 2nd deriv
                        double a2DU_ck[dim];          double b2DU_ck[localDim][dim];
                        double pDU_ck[dim];           double pLGradDU_ck[localDim][dim];
                        double qDU_ck[dim];           double qLGradDU_ck[localDim][dim];
                        double qDW_ck[dim];           double qLGradDW_ck[localDim][dim];
                        double dDU_ck[dim];           double cDU_ck[localDim][dim];
                        double dDW_ck[dim];           double cDW_ck[localDim][dim];
                        double wDU_ck[dim];           double zDU_ck[localDim][dim];
                        double wDW_ck[dim];           double zDW_ck[localDim][dim];
                        //
                        double a2DUDU_dlck[dim];      double b2DUDU_dlck[localDim][dim];
                        double pDUDU_dlck[dim];       double pLGradDUDU_dlck[localDim][dim];
                        double qDUDU_dlck[dim];       double qLGradDUDU_dlck[localDim][dim];
                        double qDUDW_dlck[dim];       double qLGradDUDW_dlck[localDim][dim];
                        double qDWDU_dlck[dim];       double qLGradDWDU_dlck[localDim][dim];
                        double qDWDW_dlck[dim];       double qLGradDWDW_dlck[localDim][dim];
                        double dDUDU_dlck[dim];       double cDUDU_dlck[localDim][dim];
                        double dDUDW_dlck[dim];       double cDUDW_dlck[localDim][dim];
                        double dDWDU_dlck[dim];       double cDWDU_dlck[localDim][dim];
                        double dDWDW_dlck[dim];       double cDWDW_dlck[localDim][dim];
                        double wDWDU_dlck[dim];       double zDUDW_dlck[localDim][dim];
                        double wDUDW_dlck[dim];       double zDWDU_dlck[localDim][dim];
                        if ( secondDeriv ) {

                            la::zero<dim>( wDU_ck );
                            for ( unsigned ga=0; ga<localDim; ++ga )
                                wDU_ck[ c ] += wLocal[ ga ] * dPhiDXi( ga, k );
                            la::zero<dim>( wDW_ck );
                            if ( c < localDim ) la::update<dim>( wDW_ck, phi( k ), aCur[ c ] );
                            la::zero<localDim,dim>( zDU_ck ); 
                            for ( unsigned ga=0; ga<localDim; ++ga )
                                for ( unsigned al=0; al<localDim; ++al )
                                    zDU_ck[ ga ][ c ] += zLocal[ al ][ ga ] * dPhiDXi( al, k ) +
                                        wLocal[ al ] * ddPhiDDXi( voigtForward[al][ga], k );
                            la::zero<localDim,dim>( zDW_ck );
                            for ( unsigned ga=0; ga<localDim; ++ga ) {
                                la::zero<dim>( zDW_ck[ ga ] );
                                if ( c < localDim ) {
                                    la::update<dim>( zDW_ck[ ga ], dPhiDXi( ga, k ), aCur[ c ] );
                                    la::update<dim>( zDW_ck[ ga ], phi( k ), bCur[ c ][ ga ] );
                                }
                            }
                            
                            la::zero<dim>( wDUDW_dlck );
                            if ( c < localDim ) wDUDW_dlck[ d ] = phi( k ) * dPhiDXi( c, l );
                            la::zero<dim>( wDWDU_dlck );
                            if ( d < localDim ) wDWDU_dlck[ d ] = phi( l ) * dPhiDXi( d, k );
                            for ( unsigned ga=0; ga<localDim; ++ga ) {
                                la::zero<dim>( zDUDW_dlck[ ga ] );
                                if ( c < localDim ) {
                                    zDUDW_dlck[ ga ][ d ] = ( phi( k ) * ddPhiDDXi( voigtForward[c][ga], l ) +
                                                              dPhiDXi( ga, k ) * dPhiDXi( c, l ) );
                                }
                                la::zero<dim>( zDWDU_dlck[ ga ] );
                                if ( d < localDim )
                                    zDWDU_dlck[ ga ][ c ] = ( phi( l ) * ddPhiDDXi( voigtForward[d][ga], k ) +
                                                              dPhiDXi( ga, l ) * dPhiDXi( d, k ) );
                            }

                            this->vectorsGradCur_( c, k, xi, phi, dPhiDXi, ddPhiDDXi,
                                                   aRef, aCur, bCur, wCur, zCur, 
                                                   pCur, pLGrad, qCur, qLGrad, dCur, cCur,
                                                   pDU_ck, pLGradDU_ck, qDU_ck, qDW_ck, qLGradDU_ck, qLGradDW_ck,
                                                   dDU_ck, dDW_ck, cDU_ck, cDW_ck );

                            double a0DU_dl[dim]; la::zero<dim>( a0DU_dl ); a0DU_dl[ d ] = dPhiDXi( 0, l );
                            double a1DU_dl[dim]; la::zero<dim>( a1DU_dl ); a1DU_dl[ d ] = dPhiDXi( 1, l );
                            double a0DU_ck[dim]; la::zero<dim>( a0DU_ck ); a0DU_ck[ c ] = dPhiDXi( 0, k );
                            double a1DU_ck[dim]; la::zero<dim>( a1DU_ck ); a1DU_ck[ c ] = dPhiDXi( 1, k );
                            SurfaceGeometer::crossprodVectorHessDof( a0Cur, a0DU_dl, a0DU_ck, zeroDim,
                                                                     a1Cur, a1DU_dl, a1DU_ck, zeroDim,
                                                                     pDUDU_dlck );
                            SurfaceGeometer::unitVectorGradDof( pLen, pDU_ck, a2Cur, a2DU_ck );
                            double sLGrad[localDim][dim];
                            la::assign<localDim,dim>( sLGrad, pLGrad );
                            la::divide<localDim,dim>( sLGrad, pLen );
                            SurfaceGeometer::unitVectorDeriv1GradDof( pLen, pDU_ck, pLGradDU_ck,
                                                                      a2Cur, a2DU_ck, sLGrad,
                                                                      b2DU_ck );

                            double b0DU_dl[localDim][dim];  la::zero<localDim,dim>( b0DU_dl );
                            double b1DU_dl[localDim][dim];  la::zero<localDim,dim>( b1DU_dl );
                            double b0DU_ck[localDim][dim];  la::zero<localDim,dim>( b0DU_ck );
                            double b1DU_ck[localDim][dim];  la::zero<localDim,dim>( b1DU_ck );
                            for ( unsigned ga=0; ga<localDim; ++ga ) {
                                b0DU_dl[ ga ][ d ] = ddPhiDDXi( voigtForward[ 0 ][ ga ], l );
                                b1DU_dl[ ga ][ d ] = ddPhiDDXi( voigtForward[ 1 ][ ga ], l );
                                b0DU_ck[ ga ][ c ] = ddPhiDDXi( voigtForward[ 0 ][ ga ], k );
                                b1DU_ck[ ga ][ c ] = ddPhiDDXi( voigtForward[ 1 ][ ga ], k );
                            }
                            SurfaceGeometer::crossprodVectorDeriv1HessDof( a0Cur, a0DU_dl, a0DU_ck, zeroDim,
                                                                           b0Cur, b0DU_dl, b0DU_ck, zeroLDimDim,
                                                                           a1Cur, a1DU_dl, a1DU_ck, zeroDim,
                                                                           b1Cur, b1DU_dl, b1DU_ck, zeroLDimDim,
                                                                           pLGradDUDU_dlck );
                            SurfaceGeometer::unitVectorHessDof( pLen, pDU_dl, pDU_ck, pDUDU_dlck,
                                                                a2Cur, a2DU_dl, a2DU_ck, a2DUDU_dlck );
                            SurfaceGeometer::unitVectorDeriv1HessDof( pLen, pDU_dl, pDU_ck, pDUDU_dlck,
                                                                      pLGradDU_dl, pLGradDU_ck, pLGradDUDU_dlck,
                                                                      a2Cur, a2DU_dl, a2DU_ck, a2DUDU_dlck, sLGrad,
                                                                      b2DUDU_dlck );
                            la::assign<dim>( qDUDU_dlck, pDUDU_dlck );
                            la::zero<dim>( qDUDW_dlck );
                            if ( c < localDim )
                                qDUDW_dlck[ d ] = phi( k ) * dPhiDXi( c, l );
                            la::zero<dim>( qDWDU_dlck );
                            if ( d < localDim )
                                qDWDU_dlck[ c ] = phi( l ) * dPhiDXi( d, k );
                            la::assign<dim>( qDWDW_dlck, zeroDim );

                            SurfaceGeometer::unitVectorHessDof( qLen, qDU_dl, qDU_ck, qDUDU_dlck,
                                                                dCur, dDU_dl, dDU_ck, dDUDU_dlck );
                            SurfaceGeometer::unitVectorHessDof( qLen, qDU_dl, qDW_ck, qDUDW_dlck,
                                                                dCur, dDU_dl, dDW_ck, dDUDW_dlck );
                            SurfaceGeometer::unitVectorHessDof( qLen, qDW_dl, qDU_ck, qDWDU_dlck,
                                                                dCur, dDW_dl, dDU_ck, dDWDU_dlck );
                            SurfaceGeometer::unitVectorHessDof( qLen, qDW_dl, qDW_ck, qDWDW_dlck,
                                                                dCur, dDW_dl, dDW_ck, dDWDW_dlck );

                            la::assign<localDim,dim>( qLGradDUDU_dlck, pLGradDUDU_dlck );
                            la::assign<localDim,dim>( qLGradDUDW_dlck, zeroLDimDim );
                            if ( c < localDim )
                                for ( unsigned ga=0; ga<localDim; ++ga )
                                    qLGradDUDW_dlck[ ga ][ d ] = (
                                        dPhiDXi( ga, k ) * dPhiDXi( c, l ) +
                                        phi( k ) * ddPhiDDXi( voigtForward[c][ga], l )
                                        );
                            la::assign<localDim,dim>( qLGradDWDU_dlck, zeroLDimDim );
                            if ( d < localDim )
                                for ( unsigned ga=0; ga<localDim; ++ga )
                                    qLGradDWDU_dlck[ ga ][ c ] = (
                                        dPhiDXi( ga, l ) * dPhiDXi( d, k ) +
                                        phi( l ) * ddPhiDDXi( voigtForward[d][ga], k )
                                        );
                            la::assign<localDim,dim>( qLGradDWDW_dlck, zeroLDimDim );

                            
                            SurfaceGeometer::unitVectorDeriv1HessDof( qLen, qDU_dl, qDU_ck, qDUDU_dlck,
                                                                      qLGradDU_dl, qLGradDU_ck, qLGradDUDU_dlck,
                                                                      dCur, dDU_dl, dDU_ck, dDUDU_dlck,
                                                                      sCur, cDUDU_dlck );
                            SurfaceGeometer::unitVectorDeriv1HessDof( qLen, qDU_dl, qDW_ck, qDUDW_dlck,
                                                                      qLGradDU_dl, qLGradDW_ck, qLGradDUDW_dlck,
                                                                      dCur, dDU_dl, dDW_ck, dDUDW_dlck,
                                                                      sCur, cDUDW_dlck );
                            SurfaceGeometer::unitVectorDeriv1HessDof( qLen, qDW_dl, qDU_ck, qDWDU_dlck,
                                                                      qLGradDW_dl, qLGradDU_ck, qLGradDWDU_dlck,
                                                                      dCur, dDW_dl, dDU_ck, dDWDU_dlck,
                                                                      sCur, cDWDU_dlck );
                            SurfaceGeometer::unitVectorDeriv1HessDof( qLen, qDW_dl, qDW_ck, qDWDW_dlck,
                                                                      qLGradDW_dl, qLGradDW_ck, qLGradDWDW_dlck,
                                                                      dCur, dDW_dl, dDW_ck, dDWDW_dlck,
                                                                      sCur, cDWDW_dlck );
                        }

                        // initialise vectors;
                        double  v[dim];       la::zero<dim>( v );
                        double  vDU_dl[dim];  la::zero<dim>( vDU_dl );
                        double  vDW_dl[dim];  la::zero<dim>( vDW_dl );

                        // set choice accordings to #vectorName
                        if      ( vectorName == "a0" ) {
                            const unsigned al = 0;
                            if ( secondDeriv == 1 ) {
                                la::update<dim>( v, dPhiDXi( al, k ), delta[ c ] );  // a0DU_ck
                                la::assign<dim>( vDU_dl, zeroDim );  // a0DUDU_dlck
                                la::assign<dim>( vDW_dl, zeroDim );  // a0DWDU_dlck
                            }
                            else if ( secondDeriv == 2 ) {
                                la::assign<dim>( v, zeroDim );  // a0DW_ck
                                la::assign<dim>( vDU_dl, zeroDim );  // a0DUDW_dlck
                                la::assign<dim>( vDW_dl, zeroDim );  // a0DWDW_dlck
                            }
                            else {
                                la::assign<dim>( v, aCur[ al ] );  // a0
                                la::update<dim>( vDU_dl, dPhiDXi( al, l ), delta[ d ] );  // a0DU_dl
                                la::assign<dim>( vDW_dl, zeroDim );  // a0DW_dl
                            }
                        }
                        else if ( vectorName == "a1" ) {
                            const unsigned al = 1;
                            if ( secondDeriv == 1 ) {
                                la::update<dim>( v, dPhiDXi( al, k ), delta[ c ] );  // a1DU_ck
                                la::assign<dim>( vDU_dl, zeroDim );  // a1DUDU_dlck
                                la::assign<dim>( vDW_dl, zeroDim );  // a1DWDU_dlck
                            }
                            else if ( secondDeriv == 2 ) {
                                la::assign<dim>( v, zeroDim );  // a1DW_ck
                                la::assign<dim>( vDU_dl, zeroDim );  // a1DUDW_dlck
                                la::assign<dim>( vDW_dl, zeroDim );  // a1DWDW_dlck
                            }
                            else {
                                la::assign<dim>( v, aCur[ al ] );  // a1
                                la::update<dim>( vDU_dl, dPhiDXi( al, l ), delta[ d ] );  // a1DU_dl
                                la::assign<dim>( vDW_dl, zeroDim );  // a1DW_dl
                            }
                        }
                        else if ( vectorName == "a2" ) {
                            if ( secondDeriv == 1 ) {
                                la::assign<dim>( v, a2DU_ck );
                                la::assign<dim>( vDU_dl, a2DUDU_dlck );  // a2DUDU_dlck
                                //vDU_dl = (
                                //    pDUDU_dlck -
                                //    ublas::inner_prod( a2Cur, pDU_ck ) * a2DU_dl -
                                //    a2DU_ck * ublas::inner_prod( a2Cur, pDU_dl ) -
                                //    a2Cur * (
                                //        ublas::inner_prod( a2DU_ck, pDU_dl ) +
                                //        ublas::inner_prod( a2Cur, pDUDU_dlck )
                                //        )
                                //    ) / pLen;
                                la::assign<dim>( vDW_dl, zeroDim );  // a2DWDU_dlck
                            }
                            else if ( secondDeriv == 2 ) {
                                la::assign<dim>( v, zeroDim );
                                la::assign<dim>( vDU_dl, zeroDim );
                                la::assign<dim>( vDW_dl, zeroDim );
                            }
                            else {
                                la::assign<dim>( v, aCur[ 2 ] );  // a2
                                la::assign<dim>( vDU_dl, a2DU_dl );  // a2DU_dl
                                la::assign<dim>( vDW_dl, zeroDim );  // a2DW_dl
                            }
                        }
                        else if ( vectorName == "b0" ) {
                            const unsigned al = 0;
                            if ( secondDeriv == 1 ) {
                                la::update<dim>( v, ddPhiDDXi( voigtForward[al][be], k ), delta[ c ] );  // b0DU_ck
                                la::assign<dim>( vDU_dl, zeroDim );  // b0DUDU_dlck
                                la::assign<dim>( vDW_dl, zeroDim );  // b0DWDU_dlck
                            }
                            else if ( secondDeriv == 2 ) {
                                la::assign<dim>( v, zeroDim );  // b0DW_ck
                                la::assign<dim>( vDU_dl, zeroDim );  // b0DUDW_dlck
                                la::assign<dim>( vDW_dl, zeroDim );  // b0DWDW_dlck
                            }
                            else {
                                la::assign<dim>( v, bCur[ al ][ be ] );
                                la::update<dim>( vDU_dl, ddPhiDDXi( voigtForward[al][be], l ), delta[ d ] );
                            }
                        }
                        else if ( vectorName == "b1" ) {
                            const unsigned al = 1;
                            if ( secondDeriv == 1 ) {
                                la::update<dim>( v, ddPhiDDXi( voigtForward[al][be], k ), delta[ c ] );  // b1DU_ck
                                la::assign<dim>( vDU_dl, zeroDim );  // b1DUDU_dlck
                                la::assign<dim>( vDW_dl, zeroDim );  // b1DWDU_dlck
                            }
                            else if ( secondDeriv == 2 ) {
                                la::assign<dim>( v, zeroDim );  // b1DW_ck
                                la::assign<dim>( vDU_dl, zeroDim );  // b1DUDW_dlck
                                la::assign<dim>( vDW_dl, zeroDim );  // b1DWDW_dlck
                            }
                            else {
                                la::assign<dim>( v, bCur[ al ][ be ] );
                                la::update<dim>( vDU_dl, ddPhiDDXi( voigtForward[al][be], l ), delta[ d ] );
                            }
                        }
                        else if ( vectorName == "b2" ) {
                            if ( secondDeriv == 1 ) {
                                la::assign<dim>( v, b2DU_ck[ be ] );  // b2DU_ck
                                la::assign<dim>( vDU_dl, b2DUDU_dlck[ be ] );  // b2DUDU_dlck
                                la::assign<dim>( vDW_dl, zeroDim );  // b2DWDU_dlck
                            }
                            else if ( secondDeriv == 2 ) {
                                la::assign<dim>( v, zeroDim );  // b2DW_ck
                                la::assign<dim>( vDU_dl, zeroDim );  // b2DUDW_dlck
                                la::assign<dim>( vDW_dl, zeroDim );  // b2DWDW_dlck
                            }
                            else {
                                la::assign<dim>( v, b2Cur[ be ] );
                                la::assign<dim>( vDU_dl, b2DU_dl[ be ] );
                                la::assign<dim>( vDW_dl, zeroDim );
                            }
                        }
                        else if ( vectorName == "w" ) {
                            if ( secondDeriv == 1 ) {
                                la::assign<dim>( v, wDU_ck );
                                la::assign<dim>( vDU_dl, zeroDim );  // wDUDU_dlck
                                la::assign<dim>( vDW_dl, wDWDU_dlck );
                            }
                            else if ( secondDeriv == 2 ) {
                                la::assign<dim>( v, wDW_ck );
                                la::assign<dim>( vDU_dl, wDUDW_dlck );  // wDUDW_dlck
                                la::assign<dim>( vDW_dl, zeroDim );  // wDWDW_dlck
                            }
                            else {
                                la::assign<dim>( v, wCur );
                                la::assign<dim>( vDU_dl, wDU_dl );
                                la::assign<dim>( vDW_dl, wDW_dl );
                            }
                        }
                        else if ( vectorName == "z" ) {
                            if ( secondDeriv == 1 ) {
                                la::assign<dim>( v, zDU_ck[ be ] );
                                la::assign<dim>( vDU_dl, zeroDim );
                                la::assign<dim>( vDW_dl, zDWDU_dlck[ be ] );
                            }
                            else if ( secondDeriv == 2 ) {
                                la::assign<dim>( v, zDW_ck[ be ] );
                                la::assign<dim>( vDU_dl, zDUDW_dlck[ be ] );
                                la::assign<dim>( vDW_dl, zeroDim );
                            }
                            else {
                                la::assign<dim>( v, zCur[ be ] );
                                la::assign<dim>( vDU_dl, zDU_dl[ be ] );
                                la::assign<dim>( vDW_dl, zDW_dl[ be ] );
                            }
                        }
                        else if ( vectorName == "p" ) {
                            if ( secondDeriv == 1 ) {
                                la::assign<dim>( v, pDU_ck );
                                la::assign<dim>( vDU_dl, pDUDU_dlck );
                            }
                            else if ( secondDeriv == 2 ) {
                                la::assign<dim>( v, zeroDim );
                                la::assign<dim>( vDU_dl, zeroDim );
                            }
                            else {
                                la::assign<dim>( v, pCur );
                                la::assign<dim>( vDU_dl, pDU_dl );
                            }
                        }
                        else if ( vectorName == "pLGrad" ) {
                            if ( secondDeriv == 1) {
                                la::assign<dim>( v, pLGradDU_ck[ be ] );
                                la::assign<dim>( vDU_dl, pLGradDUDU_dlck[ be ] );
                                la::assign<dim>( vDW_dl, zeroDim );
                            }
                            else if ( secondDeriv == 2 ) {
                                la::assign<dim>( v, zeroDim );
                                la::assign<dim>( vDU_dl, zeroDim );
                                la::assign<dim>( vDW_dl, zeroDim );
                            }
                            else {
                                la::assign<dim>( v, pLGrad[ be ] );
                                la::assign<dim>( vDU_dl, pLGradDU_dl[ be ] );
                            }
                        }
                        else if ( vectorName == "q" ) {
                            if ( secondDeriv == 1) {
                                la::assign<dim>( v, qDU_ck );
                                la::assign<dim>( vDU_dl, qDUDU_dlck );
                                la::assign<dim>( vDW_dl, qDWDU_dlck );
                            }
                            else if ( secondDeriv == 2 ) {
                                la::assign<dim>( v, qDW_ck );
                                la::assign<dim>( vDU_dl, qDUDW_dlck );
                                la::assign<dim>( vDW_dl, qDWDW_dlck );
                            }
                            else {
                                la::assign<dim>( v, qCur );
                                la::assign<dim>( vDU_dl, qDU_dl );
                                la::assign<dim>( vDW_dl, qDW_dl );
                            }
                        }
                        else if ( vectorName == "qLGrad" ) {
                            if ( secondDeriv == 1) {
                                la::assign<dim>( v, qLGradDU_ck[ be ] );
                                la::assign<dim>( vDU_dl, qLGradDUDU_dlck[ be ] );
                                la::assign<dim>( vDW_dl, qLGradDWDU_dlck[ be ] );
                            }
                            else if ( secondDeriv == 2 ) {
                                la::assign<dim>( v, qLGradDW_ck[ be ] );
                                la::assign<dim>( vDU_dl, qLGradDUDW_dlck[ be ] );
                                la::assign<dim>( vDW_dl, qLGradDWDW_dlck[ be ] );
                            }
                            else {
                                la::assign<dim>( v, qLGrad[ be ] );
                                la::assign<dim>( vDU_dl, qLGradDU_dl[ be ] );
                                la::assign<dim>( vDW_dl, qLGradDW_dl[ be ] );
                            }
                        }
                        else if ( vectorName == "d" ) {
                            if ( secondDeriv == 1 ) {
                                la::assign<dim>( v, dDU_ck );
                                la::assign<dim>( vDU_dl, dDUDU_dlck );
                                la::assign<dim>( vDW_dl, dDWDU_dlck );
                            }
                            else if ( secondDeriv == 2 ) {
                                la::assign<dim>( v, dDW_ck );
                                la::assign<dim>( vDU_dl, dDUDW_dlck );
                                la::assign<dim>( vDW_dl, dDWDW_dlck );
                            }
                            else {
                                la::assign<dim>( v, dCur );
                                la::assign<dim>( vDU_dl, dDU_dl );
                                la::assign<dim>( vDW_dl, dDW_dl );
                            }
                        }
                        else if ( vectorName == "c" ) {
                            if ( secondDeriv == 1 ) {
                                la::assign<dim>( v, cDU_ck[ be ] );
                                la::assign<dim>( vDU_dl, cDUDU_dlck[ be ] );
                                la::assign<dim>( vDW_dl, cDWDU_dlck[ be ] );
                            }
                            else if ( secondDeriv == 2 ) {
                                la::assign<dim>( v, cDW_ck[ be ] );
                                la::assign<dim>( vDU_dl, cDUDW_dlck[ be ] );
                                la::assign<dim>( vDW_dl, cDWDW_dlck[ be ] );
                            }
                            else {
                                la::assign<dim>( v, cCur[ be ] );
                                la::assign<dim>( vDU_dl, cDU_dl[ be ] );
                                la::assign<dim>( vDW_dl, cDW_dl[ be ] );
                            }
                        }
                        else if ( vectorName == "alpha" ) {
                            if ( secondDeriv ) {
                                FTL_VERIFY( false );
                            }
                            else {
                                la::assign<dim>( v, alphaCur[ be ] );
                                la::assign<dim>( vDU_dl, alphaDU_dl[ be ] );
                                la::assign<dim>( vDW_dl, alphaDW_dl[ be ] );
                            }
                        }
                        else if ( vectorName == "beta" ) {
                            if ( secondDeriv ) {
                                FTL_VERIFY( false );
                            }
                            else {
                                la::assign<dim>( v, betaCur[ be ] );
                                la::assign<dim>( vDU_dl, betaDU_dl[ be ] );
                                la::assign<dim>( vDW_dl, betaDW_dl[ be ] );
                            }
                        }
                        else {
                            FTL_VERIFY_DESCRIPTIVE( false, "Vector %s not found\n",
                                                    vectorName.c_str( ) );
                        }

                        // write to vector to file
                        tangentFD += sig * v[ a ];
                        sig *= -1.0;
                
                        // tangent
                        if ( s == 0 ) {
                            if ( dd < dim )
                                tangent = vDU_dl[ a ];
                            else
                                tangent = vDW_dl[ a ];
                        }

                        low += width;

                    }

                    // finite differencing
                    tangentFD /= width;

                    // relative tolerance
                    const double tol = std::fabs( tangentFD ) > 1.e-10 ?
                        1000.0 * width * std::pow( 10.0, std::floor( std::log10( std::fabs( tangentFD ) ) ) ) :
                        1000.0 * width;
                    // sign
                    const double signFD = std::fabs( tangentFD ) > 1.e-10 ?
                        copysign( 1.0, tangentFD ) :
                        copysign( 1.0, tangent );
                    // check
                    const bool isEqual = ( ( std::fabs( tangent - tangentFD ) < tol ) and 
                                           ( copysign( 1.0, tangent ) == signFD ) );

                    // inform user
                    out << "vector=" << vectorName
                        << " a=" << a
                        << " be=" << be
                        << " l=" << std::setw(2) << l
                        << " dd=" << dd;
                    if ( secondDeriv != 0 )
                        out << " sd=" << secondDeriv
                            << " k=" << std::setw(2) << k
                            << " cc=" << cc;

                    if ( not isEqual )
                        out << " : ERROR";
                    else  
                        out << " : OK";

                    out << " : tang=" << std::setw(13) << std::scientific << std::setprecision(6) << tangent
                        << ", tangFD=" << std::setw(13) << std::scientific << std::setprecision(6) << tangentFD
                        << ", err=" << std::setw(13) << std::scientific << std::setprecision(6) << std::fabs(tangent-tangentFD)
                        << ", tol=" << std::setw(13) << std::scientific << std::setprecision(6) << tol;
                    out << std::endl;

                }
            }
        }
    }

    return;
}

//------------------------------------------------------------------------------
//! Write nodal tangents
template< typename BELEMENT, typename MAT >
void gshell::fem::ElementStatic<BELEMENT,MAT>::writeNodalTangents(
    std::ostream & out
    )
{
    // header
    out << "-------------" << std::endl;
    out << "element="
        //<< this->giveId()
        << std::endl;
    
    // loop vertices
    for ( unsigned k = 0; k < numVertices; ++k ) {
        Node * node = this->giveNodePtr( k );
        out << "vertexId=" << node->giveId()
            << ", eVertexId=" << k
            << ", coord=" << node->giveCoordinates()
            << std::endl;
    }

    // loop support nodes
    for ( unsigned k = 0; k < this->numFunctions(); ++k ) {
        Node * node_k = this->BasisElement::giveSupportNodePtr( k );
        std::array< VecDim, localDim > tang_k;
//        this->computeTangentsAtNode( ElementStatic_::CURRENT, k,
//                                     tang_k[0], tang_k[1], NULL, NULL );
        tang_k[ 0 ] = node_k->giveTangent( 0 );
        tang_k[ 1 ] = node_k->giveTangent( 1 );
        out << "nodeId=" << node_k->giveId()
            << ", eNodeId=" << k
            << ", tang0=" << tang_k[0]
            << ", tang1=" << tang_k[1]
            << std::endl;
    }
}
