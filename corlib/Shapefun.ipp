// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//
//                 Computational Structural Mechanics Lab
//                         University of Cambridge
//
// This software is copyrighted and all rights are retained by the University of
// Cambridge. All use, disclosure, and/or reproduction of any part not expressly
// authorised by the University of Cambridge is prohibited. (C) 2010
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

//! @file Shapefun.ipp

namespace corlib {
    //------------------------------------------------------------------------------
    /** \brief Pseudo shape function on a point
     */
    template<> 
    void Shapefun<corlib::POINT, 1>::evaluate( const VecDim & xi, 
                                              VecNFun & phi ) const
    {
        return;
    }
    
    template<> 
    void Shapefun<corlib::POINT, 1>::giveInterpolationPoints( MatDimNFun & ipoints )
    {
        return;
    }

    //------------------------------------------------------------------------------
    /** \brief Linear shape function for one dimension
     *  \details <pre>
     *
     *   0----------1 --> x 
     *
     *  </pre>
     */   
    template<> 
    void Shapefun<corlib::LINE, 2>::evaluate( const VecDim & xi,
                                              VecNFun & phi ) const
    {
        phi( 0 ) = 1. - xi( 0 );
        phi( 1 ) =      xi( 0 );
        return;
    }

    template<> 
    void Shapefun<corlib::LINE, 2>::evaluateGradient( const VecDim & xi,
                                                      MatDimNFun & dphi ) const
    {
        dphi( 0, 0 ) = -1.;
        dphi( 0, 1 ) =  1.;
        return;
    }

    template<> 
    void Shapefun<corlib::LINE, 2>::giveInterpolationPoints( MatDimNFun & ipoints )
    {
        ipoints( 0, 0 ) = 0.;
        ipoints( 0, 1 ) = 1.;
        return;
    }

    //------------------------------------------------------------------------------
    /** \brief Quadratic shape function for one dimension (third node is middle node!)
     *  \details <pre>
     *
     *   0----2----1 --> x  
     *
     *  </pre>
     */   
    template<> 
    void Shapefun<corlib::LINE, 3>::evaluate( const VecDim & xi,
                                              VecNFun & phi ) const
    {
        phi( 0 ) = (1. - xi(0)) * (1. - 2.*xi(0));
        phi( 1 ) = xi(0) * (2. * xi(0) - 1.);
        phi( 2 ) = 4. * xi(0) * (1. - xi(0));
        return;
    }

    template<> 
    void Shapefun<corlib::LINE, 3>::evaluateGradient( const VecDim & xi,
                                                      MatDimNFun & dphi ) const
    {
        dphi( 0, 0 ) = 4. * xi(0) - 3.;
        dphi( 0, 1 ) = 4. * xi(0) - 1.;
        dphi( 0, 2 ) = 4. - 8. * xi(0);
        return;
    }

    template<> 
    void Shapefun<corlib::LINE, 3>::evaluateHessian( const VecDim & xi,
                                                     MatVecNFunDimDim & ddphi ) const
    {
        ddphi(0,0)( 0 ) =  4.;
        ddphi(0,0)( 1 ) =  4.;
        ddphi(0,0)( 2 ) = -8.;
    }

    template<>
    void Shapefun<corlib::LINE, 3>::giveInterpolationPoints( MatDimNFun & ipoints )
    {
        ipoints( 0, 0 ) = 0.; 
        ipoints( 0, 1 ) = 1.; 
        ipoints( 0, 2 ) = 0.5;
    }

    //------------------------------------------------------------------------------
    /** \brief Cubic shape function for a line element (0,1) 
     *  \details <pre>
     *
     *   0----2----3----1 --> x 
     *
     *  </pre>
     */       
    template<> 
    void Shapefun<corlib::LINE, 4>::evaluate( const VecDim & xi,
                                              VecNFun & phi ) const
    {
        const double zeta0 = 1. - xi(0);
        const double zeta1 =      xi(0);
        phi( 0 ) = 0.5 * zeta0 * (3. * zeta0 - 1.) * (3. * zeta0 - 2.);
        phi( 1 ) = 0.5 * zeta1 * (3. * zeta1 - 1.) * (3. * zeta1 - 2.);
        phi( 2 ) = 4.5 * zeta0 * (3. * zeta0 - 1.) * zeta1;
        phi( 3 ) = 4.5 * zeta1 * (3. * zeta1 - 1.) * zeta0;
        return;
    }

    template<> 
    void Shapefun<corlib::LINE, 4>::evaluateGradient( const VecDim & xi,
                                                      MatDimNFun & dphi ) const
    {
        const double zeta0 = 1. - xi(0);
        const double zeta1 =      xi(0);
        dphi( 0, 0 ) = -0.5 * ( (3.*zeta0 - 1.)*(3.*zeta0 - 2.) + 3.*zeta0 * (6.*zeta0 - 3.) );
        dphi( 0, 1 ) =  0.5 * ( (3.*zeta1 - 1.)*(3.*zeta1 - 2.) + 3.*zeta1 * (6.*zeta1 - 3.) );
        dphi( 0, 2 ) = -4.5 * ( zeta1 * (6.*zeta0 - 1.) - zeta0 * (3.*zeta0 - 1.) );
        dphi( 0, 3 ) =  4.5 * ( zeta0 * (6.*zeta1 - 1.) - zeta1 * (3.*zeta1 - 1.) );
        return;
    }

    template<> 
    void Shapefun<corlib::LINE, 4>::evaluateHessian( const VecDim & xi,
                                                     MatVecNFunDimDim & ddphi ) const
    {
        ddphi(0,0)( 0 ) = -9. * (3. * xi(0) - 2.);
        ddphi(0,0)( 1 ) =  9. * (3. * xi(0) - 1.);
        ddphi(0,0)( 2 ) =  9. * (9. * xi(0) - 5.);
        ddphi(0,0)( 3 ) = -9. * (9. * xi(0) - 4.);
    }

    template<>
    void Shapefun<corlib::LINE, 4>::giveInterpolationPoints( MatDimNFun & ipoints )
    {
        ipoints( 0, 0 ) = 0.; 
        ipoints( 0, 1 ) = 1.; 
        ipoints( 0, 2 ) = 1./3.;
        ipoints( 0, 3 ) = 2./3.;
    }


    //------------------------------------------------------------------------------
    /** \brief Linear shape functions on a triangle
     *  \details <pre>
     *
     *    y
     *    ^                                            
     *    |                                            
     *    2                    
     *    |`\                  
     *    |  `\                
     *    |    `\              
     *    |      `\            
     *    |        `\          
     *    0----------1 --> x   
     *
     *  </pre>
     */       
    template<> 
    void Shapefun<corlib::TRIANGLE, 3>::evaluate( const VecDim & xi,
                                                  VecNFun & phi ) const
    {
        phi( 0 ) = 1. - ( xi( 0 ) + xi( 1 ) );
        phi( 1 ) = xi( 0 );
        phi( 2 ) = xi( 1 );
        return;
    }

    template<> 
    void Shapefun<corlib::TRIANGLE, 3>::evaluateGradient( const VecDim & xi,
                                                          MatDimNFun & dphi ) const
    {
        dphi( 0, 0 ) = -1.; dphi( 0, 1 ) =  1.; dphi( 0, 2 ) = 0.;
        dphi( 1, 0 ) = -1.; dphi( 1, 1 ) =  0.; dphi( 1, 2 ) = 1.;
        return;
    }

    template<>
    void Shapefun<corlib::TRIANGLE, 3>::giveInterpolationPoints( MatDimNFun & ipoints )
    {
        ipoints( 0, 0 ) = 0.; ipoints( 0, 1 ) = 1.; ipoints( 0, 2 ) = 0.;
        ipoints( 1, 0 ) = 0.; ipoints( 1, 1 ) = 0.; ipoints( 1, 2 ) = 1.;
    }

    //------------------------------------------------------------------------------
    /** \brief Quadratic shape functions on a triangle
     *  \details <pre>
     *
     *   y
     *   ^                                            
     *   |                                            
     *   2                     
     *   |`\                  
     *   |  `\                
     *   5    `4              
     *   |      `\            
     *   |        `\          
     *   0-----3----1 --> x         
     *
     *  </pre>
     */       
    template<> 
    void Shapefun<corlib::TRIANGLE, 6>::evaluate( const VecDim & xi,
                                                  VecNFun & phi ) const
    {
        phi( 0 ) = (1. - xi(0) - xi(1)) * (1. - 2.*xi(0) - 2.*xi(1) );
        phi( 1 ) = xi(0) * ( 2. * xi(0) - 1.);
        phi( 2 ) = xi(1) * ( 2. * xi(1) - 1.);

        phi( 3 ) = 4. * xi(0) * (1. - xi(0) - xi(1) );
        phi( 4 ) = 4. * xi(0) * xi(1);
        phi( 5 ) = 4. * xi(1) * (1. - xi(0) - xi(1) );
        return;
    }

    template<> 
    void Shapefun<corlib::TRIANGLE, 6>::evaluateGradient( const VecDim & xi,
                                                          MatDimNFun & dphi ) const
    {
        dphi( 0, 0 )  = 4. * xi(0) + 4. * xi(1) - 3.;
        dphi( 1, 0 )  = dphi( 0, 0 );
        dphi( 0, 1 )  = 4. * xi( 0 ) - 1.;
        dphi( 1, 1 )  = 0.;
        dphi( 0, 2 )  = 0.;
        dphi( 1, 2 )  = 4. * xi( 1 ) - 1.;

        dphi( 0, 3 )  = 4. - 8. * xi(0) - 4. * xi(1);
        dphi( 1, 3 )  = - 4. * xi(0);
        dphi( 0, 4 )  = 4. * xi(1);
        dphi( 1, 4 )  = 4. * xi(0);
        dphi( 0, 5 )  = - 4. * xi(1);
        dphi( 1, 5 )  = 4. - 4. * xi(0) - 8. * xi(1);
        return;
    }

    template<> 
    void Shapefun<corlib::TRIANGLE,6>::evaluateHessian( const VecDim & xi,
                                                        MatVecNFunDimDim & ddphi ) const
    {
        ddphi(0, 0)(0) = 4.; ddphi(0, 1)(0) = ddphi(1, 0)(0) = 4.; ddphi(1, 1)(0) = 4.;
        ddphi(0, 0)(1) = 4.; ddphi(0, 1)(1) = ddphi(1, 0)(1) = 0.; ddphi(1, 1)(1) = 0.;
        ddphi(0, 0)(2) = 0.; ddphi(0, 1)(2) = ddphi(1, 0)(2) = 0.; ddphi(1, 1)(2) = 4.;

        ddphi(0, 0)(3) =-8.; ddphi(0, 1)(3) = ddphi(1, 0)(3) =-4.; ddphi(1, 1)(3) = 0.;
        ddphi(0, 0)(4) = 0.; ddphi(0, 1)(4) = ddphi(1, 0)(4) = 4.; ddphi(1, 1)(4) = 0.;
        ddphi(0, 0)(5) = 0.; ddphi(0, 1)(5) = ddphi(1, 0)(5) =-4.; ddphi(1, 1)(5) =-8.;
        return;
    }

    template<>
    void Shapefun<corlib::TRIANGLE, 6>::giveInterpolationPoints( MatDimNFun & ipoints )
    {
        // corner points
        ipoints( 0, 0 ) = 0.;  ipoints( 0, 1 ) = 1.;  ipoints( 0, 2 ) = 0.;
        ipoints( 1, 0 ) = 0.;  ipoints( 1, 1 ) = 0.;  ipoints( 1, 2 ) = 1.;
        // edge mid-points
        ipoints( 0, 3 ) = 0.5; ipoints( 0, 4 ) = 0.5; ipoints( 0, 5 ) = 0.;
        ipoints( 1, 3 ) = 0.;  ipoints( 1, 4 ) = 0.5; ipoints( 1, 5 ) = 0.5;
    }

    //------------------------------------------------------------------------------
    /** \brief  Cubic shape functions on a triangle
     *  \details <pre>
     *
     *   y
     *   ^                                            
     *   |                                            
     *   2                     
     *   |`\                  
     *   7  `6                
     *   |    `\              
     *   8   9  `5            
     *   |        `\          
     *   0---3---4---1 --> x         
     *
     *  </pre>
     */       
    template<> 
    void Shapefun< corlib::TRIANGLE,10 >::evaluate( const VecDim & xi,
                                                    VecNFun & phi ) const
    {
        const double zeta0 = 1. - xi(0) - xi(1);
        const double zeta1 =      xi(0);
        const double zeta2 =              xi(1);

        // corner functions (0,1,2)
        phi( 0 ) = 0.5 * zeta0 * (3. * zeta0 - 1.) * (3. * zeta0 - 2.);
        phi( 1 ) = 0.5 * zeta1 * (3. * zeta1 - 1.) * (3. * zeta1 - 2.);
        phi( 2 ) = 0.5 * zeta2 * (3. * zeta2 - 1.) * (3. * zeta2 - 2.);

        // edge functions (0-1, 1-2, 2,0)
        phi( 3 ) = 4.5 * zeta0 * (3. * zeta0 - 1.) * zeta1;
        phi( 4 ) = 4.5 * zeta1 * (3. * zeta1 - 1.) * zeta0;

        phi( 5 ) = 4.5 * zeta1 * (3. * zeta1 - 1.) * zeta2;
        phi( 6 ) = 4.5 * zeta2 * (3. * zeta2 - 1.) * zeta1;

        phi( 7 ) = 4.5 * zeta2 * (3. * zeta2 - 1.) * zeta0;
        phi( 8 ) = 4.5 * zeta0 * (3. * zeta0 - 1.) * zeta2;
        
        // bubble function
        phi( 9 ) = 27. * zeta0 * zeta1 * zeta2;
        return;
    }

    template<> 
    void Shapefun< corlib::TRIANGLE, 10 >::evaluateGradient( const VecDim & xi,
                                                             MatDimNFun & dphi ) const
    {
        const double zeta0 = 1. - xi(0) - xi(1);
        const double zeta1 =      xi(0);
        const double zeta2 =              xi(1);
        
        // corner functions
        dphi( 0, 0 )  = -0.5 * ( (3. * zeta0 - 1.) * (3. * zeta0 - 2.) + 3. * zeta0 * (6. * zeta0 - 3.) );
        dphi( 1, 0 )  = -0.5 * ( (3. * zeta0 - 1.) * (3. * zeta0 - 2.) + 3. * zeta0 * (6. * zeta0 - 3.) );
        dphi( 0, 1 )  =  0.5 * ( (3. * zeta1 - 1.) * (3. * zeta1 - 2.) + 3. * zeta1 * (6. * zeta1 - 3.) );
        dphi( 1, 1 )  =  0.;
        dphi( 0, 2 )  =  0.;
        dphi( 1, 2 )  =  0.5 * ( (3.*zeta2 - 1.)*(3.*zeta2 - 2.) + 3.*zeta2 * (6.*zeta2 - 3.) );

        // edge functions (0-1, 1-2, 2,0)
        dphi( 0, 3 )  = -4.5 * ( zeta1 * (6. * zeta0 - 1.) - zeta0 * (3. * zeta0 - 1.) );
        dphi( 1, 3 )  = -4.5 * ( zeta1 * (6. * zeta0 - 1.) );
        dphi( 0, 4 )  =  4.5 * ( zeta0 * (6. * zeta1 - 1.) - zeta1 * (3. * zeta1 - 1.) );
        dphi( 1, 4 )  = -4.5 * ( zeta1 * (3. * zeta1 - 1.) );

        dphi( 0, 5 )  =  4.5 * (6. * zeta1 - 1.) * zeta2;
        dphi( 1, 5 )  =  4.5 * (3. * zeta1 - 1.) * zeta1;
        dphi( 0, 6 )  =  4.5 * (3. * zeta2 - 1.) * zeta2;
        dphi( 1, 6 )  =  4.5 * (6. * zeta2 - 1.) * zeta1;

        dphi( 0, 7 )  = -4.5 * ( zeta2 * (3. * zeta2 - 1.) );
        dphi( 1, 7 )  =  4.5 * ( zeta0 * (6. * zeta2 - 1.) - zeta2 * (3. * zeta2 - 1.) );
        dphi( 0, 8 )  = -4.5 * ( zeta2 * (6. * zeta0 - 1.) ); 
        dphi( 1, 8 )  = -4.5 * ( zeta2 * (6. * zeta0 - 1.) - zeta0 * (3. * zeta0 - 1.) );

        // bubble
        dphi( 0, 9 )  = 27. * zeta2 * (zeta0 - zeta1);
        dphi( 1, 9 )  = 27. * zeta1 * (zeta0 - zeta2);

        return;
    }

    template<> 
    void Shapefun<corlib::TRIANGLE,10>::evaluateHessian( const VecDim & xi,
                                                         MatVecNFunDimDim & ddphi ) const
    {
        const double xi1 = xi(0);
        const double xi2 = xi(1);

        // corner functions
        ddphi( 0, 0 )( 0 ) = -9. * (3. * xi2 + 3. * xi1 - 2.);
        ddphi( 1, 1 )( 0 ) = -9. * (3. * xi2 + 3. * xi1 - 2.);
        ddphi( 1, 0 )( 0 ) = -9. * (3. * xi2 + 3. * xi1 - 2.);

        ddphi( 0, 0 )( 1 ) = 9. * (3. * xi1 - 1.);
        ddphi( 1, 1 )( 1 ) = 0.;
        ddphi( 1, 0 )( 1 ) = 0.;

        ddphi( 0, 0 )( 2 ) = 0.;
        ddphi( 1, 1 )( 2 ) = 9. * (3. * xi2 - 1.);
        ddphi( 1, 0 )( 2 ) = 0.;

        // edge functions
        ddphi( 0, 0 )( 3 ) =  9.  * (6. * xi2 +  9. * xi1 - 5.);
        ddphi( 1, 1 )( 3 ) = 27.  * xi1;
        ddphi( 1, 0 )( 3 ) =  4.5 * (6. * xi2 + 12. * xi1 - 5.);
        
        ddphi( 0, 0 )( 4 ) = -9.  * (3. * xi2 + 9. * xi1 - 4.);
        ddphi( 1, 1 )( 4 ) =  0.;
        ddphi( 1, 0 )( 4 ) = -4.5 * (6. * xi1 - 1.);

        ddphi( 0, 0 )( 5 ) = 27.  * xi2;
        ddphi( 1, 1 )( 5 ) =  0.;
        ddphi( 1, 0 )( 5 ) =  4.5 * (6. * xi1 - 1.);
        
        ddphi( 0, 0 )( 6 ) =  0.;
        ddphi( 1, 1 )( 6 ) = 27.  * xi1;
        ddphi( 1, 0 )( 6 ) =  4.5 * (6. * xi2 - 1.);

        ddphi( 0, 0 )( 7 ) =  0.;
        ddphi( 1, 1 )( 7 ) = -9.  * (9. * xi2 + 3. * xi1 - 4.);
        ddphi( 1, 0 )( 7 ) = -4.5 * (6. * xi2 - 1.);

        ddphi( 0, 0 )( 8 ) = 27.  * xi2;
        ddphi( 1, 1 )( 8 ) =  9.  * ( 9. * xi2 + 6. * xi1 - 5.);
        ddphi( 1, 0 )( 8 ) =  4.5 * (12. * xi2 + 6. * xi1 - 5.);

        // bubble function
        ddphi( 0, 0 )( 9 ) = -54. * xi2;
        ddphi( 1, 1 )( 9 ) = -54. * xi1;
        ddphi( 1, 0 )( 9 ) = -27. * (2. * xi2 + 2. * xi1 - 1.);

        // symmetry of mixed partial derivatives: d^2/(d xi_1 d xi_2) = d^2/(d xi_2 d xi_1)
        for ( unsigned s = 0; s < 10; s ++ ) ddphi( 0, 1 )( s ) = ddphi( 1, 0 )( s );
    }


    template<>
    void Shapefun< corlib::TRIANGLE, 10 >::giveInterpolationPoints( MatDimNFun & ipoints )
    {
        // corners
        ipoints( 0, 0 ) = 0.;  ipoints( 0, 1 ) = 1.;  ipoints( 0, 2 ) = 0.;
        ipoints( 1, 0 ) = 0.;  ipoints( 1, 1 ) = 0.;  ipoints( 1, 2 ) = 1.;

        // edges (0-1, 1-2, 2,0)
        ipoints( 0, 3 ) = 1./3.; ipoints( 0, 4 ) = 2./3.; 
        ipoints( 1, 3 ) = 0.;    ipoints( 1, 4 ) = 0.;    
        ipoints( 0, 5 ) = 2./3.; ipoints( 0, 6 ) = 1./3.; 
        ipoints( 1, 5 ) = 1./3.; ipoints( 1, 6 ) = 2./3.;    
        ipoints( 0, 7 ) = 0.;    ipoints( 0, 8 ) = 0.; 
        ipoints( 1, 7 ) = 2./3.; ipoints( 1, 8 ) = 1./3.;    

        // bubble
        ipoints( 0, 9 ) = 1./3.;
        ipoints( 1, 9 ) = 1./3.;
    }


    //------------------------------------------------------------------------------
    /** \brief Bilinear shape functions for a quadrilateral
     *  \details <pre>
     *
     *   y                                                                 
     *   ^                                                                 
     *   |                                                                 
     *   3-----------2      
     *   |           |      
     *   |           |      
     *   |           |      
     *   |           |       
     *   |           |        
     *   0-----------1 --> x       
     *
     *  </pre>
     */       
    template<> 
    void Shapefun< corlib::QUADRILATERAL, 4>::evaluate( const VecDim & xi,
                                                        VecNFun & phi ) const
    {
        phi( 0 ) = (1. - xi(0) ) * (1. - xi(1) );
        phi( 1 ) =       xi(0)   * (1. - xi(1) );
        phi( 2 ) =       xi(0)   *       xi(1)  ;
        phi( 3 ) = (1. - xi(0) ) *       xi(1)  ;
    }

    template<> 
    void Shapefun< corlib::QUADRILATERAL, 4>::evaluateGradient( const VecDim & xi,
                                                                MatDimNFun & dphi ) const
    {
        dphi( 0, 0 ) = - (1. - xi(1)); dphi( 1, 0 ) = - (1. - xi(0) );
        dphi( 0, 1 ) =   (1. - xi(1)); dphi( 1, 1 ) =       - xi(0)  ;
        dphi( 0, 2 ) =         xi(1) ; dphi( 1, 2 ) =         xi(0)  ;
        dphi( 0, 3 ) =       - xi(1) ; dphi( 1, 3 ) =   (1. - xi(0) );
    }

    template<> 
    void Shapefun<corlib::QUADRILATERAL,4>::evaluateHessian( const VecDim & xi,
                                                             MatVecNFunDimDim & ddphi ) const
    {
        ddphi( 0, 0 )( 0 ) =  0.;
        ddphi( 1, 1 )( 0 ) =  0.;
        ddphi( 1, 0 )( 0 ) =  1.;
        ddphi( 0, 1 )( 0 ) =  1.;

        ddphi( 0, 0 )( 1 ) =  0.;
        ddphi( 1, 1 )( 1 ) =  0.;
        ddphi( 1, 0 )( 1 ) = -1.;
        ddphi( 0, 1 )( 1 ) = -1.;

        ddphi( 0, 0 )( 2 ) =  0.;
        ddphi( 1, 1 )( 2 ) =  0.;
        ddphi( 1, 0 )( 2 ) =  1.;
        ddphi( 0, 1 )( 2 ) =  1.;

        ddphi( 0, 0 )( 3 ) =  0.;
        ddphi( 1, 1 )( 3 ) =  0.;
        ddphi( 1, 0 )( 3 ) = -1.;
        ddphi( 0, 1 )( 3 ) = -1.;
    }

    template<>
    void Shapefun< corlib::QUADRILATERAL, 4>::giveInterpolationPoints( MatDimNFun & ipoints)
    {
        ipoints( 0, 0 ) = 0.; ipoints( 1, 0 ) = 0.;
        ipoints( 0, 1 ) = 1.; ipoints( 1, 1 ) = 0.;
        ipoints( 0, 2 ) = 1.; ipoints( 1, 2 ) = 1.;
        ipoints( 0, 3 ) = 0.; ipoints( 1, 3 ) = 1.;
    }


    //------------------------------------------------------------------------------
    /** \brief Biquadratic shape functions for a quadrilateral (serendipity)
     *  \details <pre>
     *
     *   y                                                                 
     *   ^                                                                 
     *   |   
     *   3-----6-----2      
     *   |           |      
     *   |           |      
     *   7           5      
     *   |           |      
     *   |           |      
     *   0-----4-----1 --> x        
     *
     *  </pre>
     */       
    template<> 
    void Shapefun< corlib::QUADRILATERAL, 8>::evaluate( const VecDim & xi,
                                                        VecNFun & phi ) const
    {
        const double x = 2. * xi(0) - 1.;
        const double y = 2. * xi(1) - 1.;

        // mid-edge functions
        phi( 4 ) = (1. - x*x) * (1. - y) / 2.;
        phi( 5 ) = (1. - y*y) * (1. + x) / 2.;
        phi( 6 ) = (1. - x*x) * (1. + y) / 2.;
        phi( 7 ) = (1. - y*y) * (1. - x) / 2.;

        // corner shape functions
        phi( 0 ) = (1. - xi(0) ) * (1. - xi(1) ) - ( phi(7) + phi(4) )/2.;
        phi( 1 ) =       xi(0)   * (1. - xi(1) ) - ( phi(4) + phi(5) )/2.;
        phi( 2 ) =       xi(0)   *       xi(1)   - ( phi(5) + phi(6) )/2.;
        phi( 3 ) = (1. - xi(0) ) *       xi(1)   - ( phi(6) + phi(7) )/2.;
    }

    template<> 
    void Shapefun< corlib::QUADRILATERAL, 8>::evaluateGradient( const VecDim & xi,
                                                                MatDimNFun & dphi ) const
    {
        const double x = 2. * xi(0) - 1.;  // x = 2 xi  - 1
        const double y = 2. * xi(1) - 1.;  // y = 2 eta - 1
        // const double j = 2.; // d x / dxi = d y / d eta = 2

        // mid-edge functions
        dphi( 0, 4 ) = - 2. * x * (1. - y); dphi( 1, 4 ) = - (1. - x*x);
        dphi( 0, 5 ) =   (1. - y*y);        dphi( 1, 5 ) = - 2. * (1. + x) * y;
        dphi( 0, 6 ) = - 2. * x * (1. + y); dphi( 1, 6 ) =   (1. - x*x);
        dphi( 0, 7 ) = - (1. - y*y);        dphi( 1, 7 ) = - 2. * (1. - x) * y;

        // corner functions
        dphi( 0, 0 ) = - (1. - xi(1)) - ( dphi( 0, 7 ) + dphi( 0, 4 ) ) / 2.;
        dphi( 0, 1 ) =   (1. - xi(1)) - ( dphi( 0, 4 ) + dphi( 0, 5 ) ) / 2.;
        dphi( 0, 2 ) =         xi(1)  - ( dphi( 0, 5 ) + dphi( 0, 6 ) ) / 2.;
        dphi( 0, 3 ) =       - xi(1)  - ( dphi( 0, 6 ) + dphi( 0, 7 ) ) / 2.;

        dphi( 1, 0 ) = - (1. - xi(0) )- ( dphi( 1, 7 ) + dphi( 1, 4 ) ) / 2.;
        dphi( 1, 1 ) =       - xi(0)  - ( dphi( 1, 4 ) + dphi( 1, 5 ) ) / 2.;
        dphi( 1, 2 ) =         xi(0)  - ( dphi( 1, 5 ) + dphi( 1, 6 ) ) / 2.;
        dphi( 1, 3 ) =   (1. - xi(0) )- ( dphi( 1, 6 ) + dphi( 1, 7 ) ) / 2.;
    }

    template<>
    void Shapefun< corlib::QUADRILATERAL, 8>::giveInterpolationPoints( MatDimNFun & ipoints)
    {
        // corners
        ipoints( 0, 0 ) = 0.; ipoints( 1, 0 ) = 0.;
        ipoints( 0, 1 ) = 1.; ipoints( 1, 1 ) = 0.;
        ipoints( 0, 2 ) = 1.; ipoints( 1, 2 ) = 1.;
        ipoints( 0, 3 ) = 0.; ipoints( 1, 3 ) = 1.;

        // edge mid-points
        ipoints( 0, 4 ) = 0.5; ipoints( 1, 4 ) = 0.;
        ipoints( 0, 5 ) = 1.;  ipoints( 1, 5 ) = 0.5;
        ipoints( 0, 6 ) = 0.5; ipoints( 1, 6 ) = 1.;
        ipoints( 0, 7 ) = 0.;  ipoints( 1, 7 ) = 0.5;

    }

    //------------------------------------------------------------------------------
    /** \brief  Biquadratic shape functions for a quadrilateral
     *  \details <pre>
     *
     *   y                                                                 
     *   ^                                                                 
     *   |   
     *   3-----6-----2     
     *   |           |     
     *   |           |     
     *   7     8     5     
     *   |           |     
     *   |           |     
     *   0-----4-----1 --> x            
     *
     *  </pre>
     */       
    template<> 
    void Shapefun< corlib::QUADRILATERAL, 9>::evaluate( const VecDim & xi,
                                                        VecNFun & phi ) const
    {
        // corner functions
        phi( 0 ) = ( 2.*xi(0)*xi(0) - 3.*xi(0) + 1.) * ( 2.*xi(1)*xi(1) - 3.*xi(1) + 1.);
        phi( 1 ) = ( 2.*xi(0)*xi(0) -    xi(0)     ) * ( 2.*xi(1)*xi(1) - 3.*xi(1) + 1.);
        phi( 2 ) = ( 2.*xi(0)*xi(0) -    xi(0)     ) * ( 2.*xi(1)*xi(1) -    xi(1)     );
        phi( 3 ) = ( 2.*xi(0)*xi(0) - 3.*xi(0) + 1.) * ( 2.*xi(1)*xi(1) -    xi(1)     );

        // mid-edge functions
        phi( 4 ) = (-4.*xi(0)*xi(0) + 4.*xi(0)     ) * ( 2.*xi(1)*xi(1) - 3.*xi(1) + 1.);
        phi( 5 ) = ( 2.*xi(0)*xi(0) -    xi(0)     ) * (-4.*xi(1)*xi(1) + 4.*xi(1)     );
        phi( 6 ) = (-4.*xi(0)*xi(0) + 4.*xi(0)     ) * ( 2.*xi(1)*xi(1) -    xi(1)     );
        phi( 7 ) = ( 2.*xi(0)*xi(0) - 3.*xi(0) + 1.) * (-4.*xi(1)*xi(1) + 4.*xi(1)     );

        // bubble
        phi( 8 ) = (-4.*xi(0)*xi(0) + 4.*xi(0)     ) * (-4.*xi(1)*xi(1) + 4.*xi(1)     );
    }

    template<> 
    void Shapefun< corlib::QUADRILATERAL, 9>::evaluateGradient( const VecDim & xi,
                                                                MatDimNFun & dphi ) const
    {
        // d() / d xi_1
        dphi(0, 0) = ( 4.*xi(0) - 3.) * ( 2.*xi(1)*xi(1) - 3.*xi(1) + 1.);
        dphi(0, 1) = ( 4.*xi(0) - 1.) * ( 2.*xi(1)*xi(1) - 3.*xi(1) + 1.);
        dphi(0, 2) = ( 4.*xi(0) - 1.) * ( 2.*xi(1)*xi(1) -    xi(1)     );
        dphi(0, 3) = ( 4.*xi(0) - 3.) * ( 2.*xi(1)*xi(1) -    xi(1)     );
        dphi(0, 4) = (-8.*xi(0) + 4.) * ( 2.*xi(1)*xi(1) - 3.*xi(1) + 1.);
        dphi(0, 5) = ( 4.*xi(0) - 1.) * (-4.*xi(1)*xi(1) + 4.*xi(1)     );
        dphi(0, 6) = (-8.*xi(0) + 4.) * ( 2.*xi(1)*xi(1) -    xi(1)     );
        dphi(0, 7) = ( 4.*xi(0) - 3.) * (-4.*xi(1)*xi(1) + 4.*xi(1)     );
        dphi(0, 8) = (-8.*xi(0) + 4.) * (-4.*xi(1)*xi(1) + 4.*xi(1)     );

        // d() / d xi_2
        dphi(1, 0) = ( 2.*xi(0)*xi(0) - 3.*xi(0) + 1.) * ( 4.*xi(1) - 3.);
        dphi(1, 1) = ( 2.*xi(0)*xi(0) -    xi(0)     ) * ( 4.*xi(1) - 3.);
        dphi(1, 2) = ( 2.*xi(0)*xi(0) -    xi(0)     ) * ( 4.*xi(1) - 1.);
        dphi(1, 3) = ( 2.*xi(0)*xi(0) - 3.*xi(0) + 1.) * ( 4.*xi(1) - 1.);
        dphi(1, 4) = (-4.*xi(0)*xi(0) + 4.*xi(0)     ) * ( 4.*xi(1) - 3.);
        dphi(1, 5) = ( 2.*xi(0)*xi(0) -    xi(0)     ) * (-8.*xi(1) + 4.);
        dphi(1, 6) = (-4.*xi(0)*xi(0) + 4.*xi(0)     ) * ( 4.*xi(1) - 1.);
        dphi(1, 7) = ( 2.*xi(0)*xi(0) - 3.*xi(0) + 1.) * (-8.*xi(1) + 4.);
        dphi(1, 8) = (-4.*xi(0)*xi(0) + 4.*xi(0)     ) * (-8.*xi(1) + 4.);

    }

    template<> 
    void Shapefun<corlib::QUADRILATERAL,9>::evaluateHessian( const VecDim & xi,
                                                             MatVecNFunDimDim & ddphi ) const
    {
        const double xi1 = xi(0);
        const double xi2 = xi(1);

        // corner functions
        ddphi( 0, 0 )( 0 ) = 4. * (xi2 - 1.) * (2. * xi2 - 1.);
        ddphi( 1, 1 )( 0 ) = 4. * (xi1 - 1.) * (2. * xi1 - 1.);
        ddphi( 1, 0 )( 0 ) = (4. * xi1 - 3.) * (4. * xi2 - 3.);

        ddphi( 0, 0 )( 1 ) = 4. * (xi2 - 1.) * (2. * xi2 - 1.);
        ddphi( 1, 1 )( 1 ) = 4. * xi1 * (2. * xi1 - 1.);
        ddphi( 1, 0 )( 1 ) = (4. * xi1 - 1.) * (4. * xi2 - 3.);

        ddphi( 0, 0 )( 2 ) = 4. * xi2 * (2. * xi2 - 1.);
        ddphi( 1, 1 )( 2 ) = 4. * xi1 * (2. * xi1 - 1.);
        ddphi( 1, 0 )( 2 ) = (4. * xi1 - 1.) * (4. * xi2 - 1.);

        ddphi( 0, 0 )( 3 ) = 4. * xi2 * (2. * xi2 - 1.);
        ddphi( 1, 1 )( 3 ) = 4. * (xi1 - 1.) * (2. * xi1 - 1.);
        ddphi( 1, 0 )( 3 ) = (4. * xi1 - 3.) * (4. * xi2 - 1.);

        // edge functions
        ddphi( 0, 0 )( 4 ) =  -8. * (xi2 - 1.) * (2. * xi2 - 1.);
        ddphi( 1, 1 )( 4 ) = -16. * (xi1 - 1.) * xi1; 
        ddphi( 1, 0 )( 4 ) =  -4. * (2. * xi1 - 1.) * (4. * xi2 - 3.);

        ddphi( 0, 0 )( 5 ) = -16. * (xi2 - 1.) * xi2;
        ddphi( 1, 1 )( 5 ) =  -8. * xi1 * (2. * xi1 - 1.);
        ddphi( 1, 0 )( 5 ) =  -4. * (4. * xi1 - 1.) * (2. * xi2 - 1.);

        ddphi( 0, 0 )( 6 ) =  -8. * xi2 * (2. * xi2 - 1.);
        ddphi( 1, 1 )( 6 ) = -16. * (xi1 - 1.) * xi1;
        ddphi( 1, 0 )( 6 ) =  -4. * (2. * xi1 - 1.) * (4. * xi2 - 1.);

        ddphi( 0, 0 )( 7 ) = -16. * (xi2 - 1.) * xi2;
        ddphi( 1, 1 )( 7 ) =  -8. * (xi1 - 1.) * (2. * xi1 - 1.);
        ddphi( 1, 0 )( 7 ) =  -4. * (4. * xi1 - 3.) * (2. * xi2 - 1.);

        // bubble functions
        ddphi( 0, 0 )( 8 ) = 32. * (xi2 - 1.) * xi2;
        ddphi( 1, 1 )( 8 ) = 32. * (xi1 - 1.) * xi1;
        ddphi( 1, 0 )( 8 ) = 16. * (2. * xi1 - 1.) * (2. * xi2 - 1.);

        // symmetry of mixed partial derivatives: d^2/(d xi_1 d xi_2) = d^2/(d xi_2 d xi_1)
        for ( unsigned s = 0; s < 9; s ++ ) ddphi( 0, 1 )( s ) = ddphi( 1, 0 )( s );
    }

    template<>
    void Shapefun< corlib::QUADRILATERAL, 9>::giveInterpolationPoints( MatDimNFun & ipoints)
    {
        // corners
        ipoints( 0, 0 ) = 0.;   ipoints( 1, 0 ) = 0.;
        ipoints( 0, 1 ) = 1.;   ipoints( 1, 1 ) = 0.;
        ipoints( 0, 2 ) = 1.;   ipoints( 1, 2 ) = 1.;
        ipoints( 0, 3 ) = 0.;   ipoints( 1, 3 ) = 1.;

        // edge mid-points
        ipoints( 0, 4 ) = 0.5;  ipoints( 1, 4 ) = 0.;
        ipoints( 0, 5 ) = 1.;   ipoints( 1, 5 ) = 0.5;
        ipoints( 0, 6 ) = 0.5;  ipoints( 1, 6 ) = 1.;
        ipoints( 0, 7 ) = 0.;   ipoints( 1, 7 ) = 0.5;

        // center
        ipoints( 0, 8 ) = 0.5;  ipoints( 1, 8 ) = 0.5;
    }

    //------------------------------------------------------------------------------
    /** \brief Linear shape functions for a tetrahedron
     *  \details <pre>
     *                  z
     *                ,/`
     *               /
     *             3                      
     *           ,/|`\                    
     *         ,/  |  `\                  
     *       ,/    '.   `\                
     *     ,/       |     `\              
     *   ,/         |       `\            
     *  0-----------'.--------2 --> y     
     *   `\.         |      ,/            
     *      `\.      |    ,/             
     *         `\.   '. ,/               
     *            `\. |/                 
     *              `1                  
     *                `\.
     *                  ` x
     *
     *  </pre>
     */       
    template<> 
    void Shapefun< corlib::TETRAHEDRON, 4 > :: evaluate( const VecDim & xi,
                                                         VecNFun & phi ) const
    {
        phi( 0 ) = 1. - xi( 0 ) - xi( 1 ) - xi( 2 );
        phi( 1 ) = xi( 0 );
        phi( 2 ) = xi( 1 );
        phi( 3 ) = xi( 2 );
    }

    template<> 
    void Shapefun< corlib::TETRAHEDRON, 4 > :: evaluateGradient( const VecDim & xi,
                                                                 MatDimNFun & dphi ) const
    {
        dphi( 0, 0 ) = -1.; dphi( 0, 1 ) = 1.; dphi( 0, 2 ) = 0.; dphi( 0, 3 ) = 0.;
        dphi( 1, 0 ) = -1.; dphi( 1, 1 ) = 0.; dphi( 1, 2 ) = 1.; dphi( 1, 3 ) = 0.;
        dphi( 2, 0 ) = -1.; dphi( 2, 1 ) = 0.; dphi( 2, 2 ) = 0.; dphi( 2, 3 ) = 1.;
    }

    template<>
    void Shapefun< corlib::TETRAHEDRON, 4 >::giveInterpolationPoints( MatDimNFun & ipoints )
    {
        ipoints( 0, 0 ) = 0.; ipoints( 1, 0 ) = 0.; ipoints( 2, 0 ) = 0.;
        ipoints( 0, 1 ) = 1.; ipoints( 1, 1 ) = 0.; ipoints( 2, 1 ) = 0.;
        ipoints( 0, 2 ) = 0.; ipoints( 1, 2 ) = 1.; ipoints( 2, 2 ) = 0.;
        ipoints( 0, 3 ) = 0.; ipoints( 1, 3 ) = 0.; ipoints( 2, 3 ) = 1.;
    }

    //------------------------------------------------------------------------------
    /** \brief Quadratic shape functions for a tetrahedron 
     *  \details <pre>
     *                   z
     *                ,/`
     *               /
     *             3                      
     *           ,/|`\                    
     *         ,/  |  `\                  
     *       ,7    '.   `9                
     *     ,/       8     `\              
     *   ,/         |       `\            
     *  0--------6--'.--------2 --> y     
     *   `\.         |      ,/            
     *      `\.      |    ,5             
     *         `4.   '. ,/               
     *            `\. |/                 
     *              `1                  
     *                `\.
     *                  ` x
     *
     *  </pre>
     */       
    template<> 
    void Shapefun< corlib::TETRAHEDRON, 10> :: evaluate( const VecDim & xi,
                                                         VecNFun & phi ) const
    {
        // Baricentric coordinates
        const double zeta1 = xi[0];
        const double zeta2 = xi[1];
        const double zeta3 = xi[2];
        const double zeta0 = 1. - zeta1 - zeta2 - zeta3;

        phi[0] = zeta0 * (2. * zeta0 - 1.);
        phi[1] = zeta1 * (2. * zeta1 - 1.);
        phi[2] = zeta2 * (2. * zeta2 - 1.);
        phi[3] = zeta3 * (2. * zeta3 - 1.);
        phi[4] = 4. * zeta0 * zeta1;
        phi[5] = 4. * zeta1 * zeta2;
        phi[6] = 4. * zeta2 * zeta0;
        phi[7] = 4. * zeta0 * zeta3;
        phi[8] = 4. * zeta1 * zeta3;
        phi[9] = 4. * zeta2 * zeta3;

        return;
    }

    template<> 
    void Shapefun< corlib::TETRAHEDRON, 10> :: evaluateGradient( const VecDim & xi,
                                                                 MatDimNFun & dphi ) const
    {
        // Barycentric coordinates
        const double zeta1 = xi[0];
        const double zeta2 = xi[1];
        const double zeta3 = xi[2];
        const double zeta0 = 1.-zeta1-zeta2-zeta3;

        // Derivatives of Barycentric with respect to Cartesian coordinates
        const double dzeta0dxi = -1.; const double dzeta0deta = -1.; const double dzeta0dzeta = -1.;
        const double dzeta1dxi =  1.; const double dzeta1deta =  0.; const double dzeta1dzeta =  0.;
        const double dzeta2dxi =  0.; const double dzeta2deta =  1.; const double dzeta2dzeta =  0.;
        const double dzeta3dxi =  0.; const double dzeta3deta =  0.; const double dzeta3dzeta =  1.;

        // d()/dxi
        dphi(0, 0) = (4. * zeta0 - 1.) * dzeta0dxi;
        dphi(0, 1) = (4. * zeta1 - 1.) * dzeta1dxi;
        dphi(0, 2) = (4. * zeta2 - 1.) * dzeta2dxi;
        dphi(0, 3) = (4. * zeta3 - 1.) * dzeta3dxi;
        dphi(0, 4) = 4. * (zeta0 * dzeta1dxi + dzeta0dxi * zeta1);
        dphi(0, 5) = 4. * (zeta1 * dzeta2dxi + dzeta1dxi * zeta2);
        dphi(0, 6) = 4. * (zeta0 * dzeta2dxi + dzeta0dxi * zeta2);
        dphi(0, 7) = 4. * (zeta0 * dzeta3dxi + dzeta0dxi * zeta3);
        dphi(0, 8) = 4. * (zeta1 * dzeta3dxi + dzeta1dxi * zeta3);
        dphi(0, 9) = 4. * (zeta2 * dzeta3dxi + dzeta2dxi * zeta3);

        // d()/deta
        dphi(1, 0) = (4. * zeta0 - 1.) * dzeta0deta;
        dphi(1, 1) = (4. * zeta1 - 1.) * dzeta1deta;
        dphi(1, 2) = (4. * zeta2 - 1.) * dzeta2deta;
        dphi(1, 3) = (4. * zeta3 - 1.) * dzeta3deta;
        dphi(1, 4) = 4. * (zeta0 * dzeta1deta + dzeta0deta * zeta1);
        dphi(1, 5) = 4. * (zeta1 * dzeta2deta + dzeta1deta * zeta2);
        dphi(1, 6) = 4. * (zeta0 * dzeta2deta + dzeta0deta * zeta2);
        dphi(1, 7) = 4. * (zeta0 * dzeta3deta + dzeta0deta * zeta3);
        dphi(1, 8) = 4. * (zeta1 * dzeta3deta + dzeta1deta * zeta3);
        dphi(1, 9) = 4. * (zeta2 * dzeta3deta + dzeta2deta * zeta3);


        // d()/dzeta
        dphi(2, 0) = (4. * zeta0 - 1.) * dzeta0dzeta;
        dphi(2, 1) = (4. * zeta1 - 1.) * dzeta1dzeta;
        dphi(2, 2) = (4. * zeta2 - 1.) * dzeta2dzeta;
        dphi(2, 3) = (4. * zeta3 - 1.) * dzeta3dzeta;
        dphi(2, 4) = 4. * (zeta0 * dzeta1dzeta + dzeta0dzeta * zeta1);
        dphi(2, 5) = 4. * (zeta1 * dzeta2dzeta + dzeta1dzeta * zeta2);
        dphi(2, 6) = 4. * (zeta0 * dzeta2dzeta + dzeta0dzeta * zeta2);
        dphi(2, 7) = 4. * (zeta0 * dzeta3dzeta + dzeta0dzeta * zeta3);
        dphi(2, 8) = 4. * (zeta1 * dzeta3dzeta + dzeta1dzeta * zeta3);
        dphi(2, 9) = 4. * (zeta2 * dzeta3dzeta + dzeta2dzeta * zeta3);

        return;
    }

    template<>
    void Shapefun< corlib::TETRAHEDRON, 10>::giveInterpolationPoints( MatDimNFun & ipoints )
    {
        ipoints( 0, 0 ) = 0.;  ipoints( 1, 0 ) = 0.;  ipoints( 2, 0 ) = 0.;
        ipoints( 0, 1 ) = 1.;  ipoints( 1, 1 ) = 0.;  ipoints( 2, 1 ) = 0.;
        ipoints( 0, 2 ) = 0.;  ipoints( 1, 2 ) = 1.;  ipoints( 2, 2 ) = 0.;
        ipoints( 0, 3 ) = 0.;  ipoints( 1, 3 ) = 0.;  ipoints( 2, 3 ) = 1.;

        ipoints( 0, 4 ) = 0.5; ipoints( 1, 4 ) = 0.;  ipoints( 2, 4 ) = 0.;
        ipoints( 0, 5 ) = 0.5; ipoints( 1, 5 ) = 0.5; ipoints( 2, 5 ) = 0.;
        ipoints( 0, 6 ) = 0.;  ipoints( 1, 6 ) = 0.5; ipoints( 2, 6 ) = 0.;
        ipoints( 0, 7 ) = 0.;  ipoints( 1, 7 ) = 0.;  ipoints( 2, 7 ) = 0.5;
        ipoints( 0, 8 ) = 0.5; ipoints( 1, 8 ) = 0.;  ipoints( 2, 8 ) = 0.5;
        ipoints( 0, 9 ) = 0.;  ipoints( 1, 9 ) = 0.5; ipoints( 2, 9 ) = 0.5;
        return;
    }

    //------------------------------------------------------------------------------
    /** \brief Trilinear shape functions for a hexahedron
     *  \details <pre>
     *
     *     z
     *   ^
     *   |                                                                        
     *   |                                                                        
     *   4----------7            
     *   |\         |\           
     *   | \        | \          
     *   |  \       |  \         
     *   |   5----------6        
     *   |   |      |   |        
     *   0---|------3---|--> y   
     *    \  |       \  |        
     *     \ |        \ |        
     *      \|         \|        
     *       1----------2        
     *        \
     *         \
     *         `'
     *           x
     *
     *  </pre>
     */       
    template<> 
    void Shapefun< corlib::HEXAHEDRON, 8 >::evaluate( const VecDim & xi,
                                                      VecNFun & phi ) const
    {
        phi( 0 ) = (1. - xi(0)) * (1. - xi(1)) * (1. - xi(2));
        phi( 1 ) =       xi(0)  * (1. - xi(1)) * (1. - xi(2));
        phi( 2 ) =       xi(0)  *       xi(1)  * (1. - xi(2));
        phi( 3 ) = (1. - xi(0)) *       xi(1)  * (1. - xi(2));

        phi( 4 ) = (1. - xi(0)) * (1. - xi(1)) * xi(2);
        phi( 5 ) =       xi(0)  * (1. - xi(1)) * xi(2);
        phi( 6 ) =       xi(0)  *       xi(1)  * xi(2);
        phi( 7 ) = (1. - xi(0)) *       xi(1)  * xi(2);
    }

    template<> 
    void Shapefun< corlib::HEXAHEDRON, 8 >::evaluateGradient( const VecDim & xi,
                                                              MatDimNFun & dphi ) const
    {
        // d() / d xi_1
        dphi( 0, 0 ) = - (1. - xi(1)) * (1. - xi(2));
        dphi( 0, 1 ) =   (1. - xi(1)) * (1. - xi(2));
        dphi( 0, 2 ) =         xi(1)  * (1. - xi(2));
        dphi( 0, 3 ) =       - xi(1)  * (1. - xi(2));
        dphi( 0, 4 ) = - (1. - xi(1)) * xi(2);
        dphi( 0, 5 ) =   (1. - xi(1)) * xi(2);
        dphi( 0, 6 ) =         xi(1)  * xi(2);
        dphi( 0, 7 ) =       - xi(1)  * xi(2);

        // d() / d xi_2
        dphi( 1, 0 ) = - (1. - xi(0)) * (1. - xi(2));
        dphi( 1, 1 ) =       - xi(0)  * (1. - xi(2));
        dphi( 1, 2 ) =         xi(0)  * (1. - xi(2));
        dphi( 1, 3 ) =   (1. - xi(0)) * (1. - xi(2));
        dphi( 1, 4 ) = - (1. - xi(0)) * xi(2);
        dphi( 1, 5 ) =       - xi(0)  * xi(2);
        dphi( 1, 6 ) =         xi(0)  * xi(2);
        dphi( 1, 7 ) =   (1. - xi(0)) * xi(2);

        // d() / d xi_3
        dphi( 2, 0 ) = - (1. - xi(0)) * (1. - xi(1));
        dphi( 2, 1 ) =       - xi(0)  * (1. - xi(1));
        dphi( 2, 2 ) =       - xi(0)  *       xi(1) ;
        dphi( 2, 3 ) = - (1. - xi(0)) *       xi(1) ;
        dphi( 2, 4 ) =   (1. - xi(0)) * (1. - xi(1));
        dphi( 2, 5 ) =         xi(0)  * (1. - xi(1));
        dphi( 2, 6 ) =         xi(0)  *       xi(1) ;
        dphi( 2, 7 ) =   (1. - xi(0)) *       xi(1) ;

    }

    //------------------------------------------------------------------------------
    /** \brief  Tri-quadratic serendipity shape functions on hexahedron \f$[0,1]^3\f$ 
     *  \details VTK_QUADRATIC_HEXAHEDRON (=25) ordering   
     *  <pre>
     *
     *     z
     *    ^
     *    |                                                                        
     *    |                                                                        
     *    4----15----7          
     *    |\         |\         
     *    |12        | 14       
     *   16  \       19 \        
     *    |   5----13----6      
     *    |   |      |   |      
     *    0---|-11---3 --|-----> y
     *     \ 17       \  18     
     *      8 |        10|      
     *       \|         \|      
     *        1-----9----2      
     *         \
     *          \
     *          `'
     *            x
     *
     *  </pre>
     */       
    template<> 
    void Shapefun< corlib::HEXAHEDRON, 20 >::evaluate( const VecDim & xi,
                                                       VecNFun & phi ) const
    {
        // map from [0,1] to [-1,1]
        const double r = 2.0 * xi( 0 ) - 1.0;
        const double s = 2.0 * xi( 1 ) - 1.0;
        const double t = 2.0 * xi( 2 ) - 1.0;

        // auxiliary variables
        const double rm  = 1.0 - r;
        const double rp  = 1.0 + r;
        const double sm  = 1.0 - s;
        const double sp  = 1.0 + s;
        const double tm  = 1.0 - t;
        const double tp  = 1.0 + t;
        const double rrm = 1.0 - r*r;
        const double ssm = 1.0 - s*s;
        const double ttm = 1.0 - t*t;

        // shape functions associated to vertex nodes k=1,...,8
        // N^k = 1/8 (1 + r^k r) (1 + s^k s) (1 + t^k k)
        //           (r^k r + s^k s + t^k t - 2)
        // with r^k,s^k,t^k = -1,+1
        // [Zienkiewicz, Methode der Finiten Elemente, Hanser, 1975]
        // However, here the slightly different notation is used
        // N^k = 1/8 (1 + r^k r) (1 + s^k s) (1 + t^k k)
        //           ( (1 + r^k r) + (1 + s^k s) + (1 + t^k t) - 2 - 3)
        phi(0) = 0.125*rm*sm*tm*(rm+sm+tm-5.0);
        phi(1) = 0.125*rp*sm*tm*(rp+sm+tm-5.0);
        phi(2) = 0.125*rp*sp*tm*(rp+sp+tm-5.0);
        phi(3) = 0.125*rm*sp*tm*(rm+sp+tm-5.0);
        phi(4) = 0.125*rm*sm*tp*(rm+sm+tp-5.0);
        phi(5) = 0.125*rp*sm*tp*(rp+sm+tp-5.0);
        phi(6) = 0.125*rp*sp*tp*(rp+sp+tp-5.0);
        phi(7) = 0.125*rm*sp*tp*(rm+sp+tp-5.0);

        // shape functions associated to middle nodes on edges k=8,...,19
        // N^k = 1/4 (1 - r r) (1 + s^k s) (1 + t^k t)
        // with r^k=0, s^k,t^k = -1,+1
        // analogously for s^k,t^k=0
        // [Zienkiewicz, Methode der Finiten Elemente, Hanser, 1975]
        phi(8) = 0.25*rrm*sm*tm;
        phi(9) = 0.25*rp*ssm*tm;
        phi(10) = 0.25*rrm*sp*tm;
        phi(11) = 0.25*rm*ssm*tm;

        phi(12) = 0.25*rrm*sm*tp;
        phi(13) = 0.25*rp*ssm*tp;
        phi(14) = 0.25*rrm*sp*tp;
        phi(15) = 0.25*rm*ssm*tp;

        phi(16) = 0.25*rm*sm*ttm;
        phi(17) = 0.25*rp*sm*ttm;
        phi(18) = 0.25*rp*sp*ttm;
        phi(19) = 0.25*rm*sp*ttm;
    }

    template<> 
    void Shapefun< corlib::HEXAHEDRON, 20 >::evaluateGradient( const VecDim & xi,
                                                               MatDimNFun & dphi ) const
    {
        // map from [0,1] to [-1,1]
        const double r = 2.0 * xi( 0 ) - 1.0;
        const double s = 2.0 * xi( 1 ) - 1.0;
        const double t = 2.0 * xi( 2 ) - 1.0;

        // auxiliary variables
        const double rm  = 1.0 - r;
        const double rp  = 1.0 + r;
        const double sm  = 1.0 - s;
        const double sp  = 1.0 + s;
        const double tm  = 1.0 - t;
        const double tp  = 1.0 + t;
        const double rrm = 1.0 - r*r;
        const double ssm = 1.0 - s*s;
        const double ttm = 1.0 - t*t;

        // corners
        dphi(0,0) = -0.125*sm*tm*(2.0*rm+sm+tm-5.0);  // d/dr
        dphi(1,0) = -0.125*rm*tm*(rm+2.0*sm+tm-5.0);  // d/ds
        dphi(2,0) = -0.125*rm*sm*(rm+sm+2.0*tm-5.0);  // d/dt
        dphi(0,1) = 0.125*sm*tm*(2.0*rp+sm+tm-5.0);
        dphi(1,1) = -0.125*rp*tm*(rp+2.0*sm+tm-5.0);
        dphi(2,1) = -0.125*rp*sm*(rp+sm+2.0*tm-5.0);
        dphi(0,2) = 0.125*sp*tm*(2.0*rp+sp+tm-5.0);
        dphi(1,2) = 0.125*rp*tm*(rp+2.0*sp+tm-5.0);
        dphi(2,2) = -0.125*rp*sp*(rp+sp+2.0*tm-5.0);
        dphi(0,3) = -0.125*sp*tm*(2.0*rm+sp+tm-5.0);
        dphi(1,3) = 0.125*rm*tm*(rm+2.0*sp+tm-5.0);
        dphi(2,3) = -0.125*rm*sp*(rm+sp+2.0*tm-5.0);
        dphi(0,4) = -0.125*sm*tp*(2.0*rm+sm+tp-5.0);
        dphi(1,4) = -0.125*rm*tp*(rm+2.0*sm+tp-5.0);
        dphi(2,4) = 0.125*rm*sm*(rm+sm+2.0*tp-5.0);
        dphi(0,5) = 0.125*sm*tp*(2.0*rp+sm+tp-5.0);
        dphi(1,5) = -0.125*rp*tp*(rp+2.0*sm+tp-5.0);
        dphi(2,5) = 0.125*rp*sm*(rp+sm+2.0*tp-5.0);
        dphi(0,6) = 0.125*sp*tp*(2.0*rp+sp+tp-5.0);
        dphi(1,6) = 0.125*rp*tp*(rp+2.0*sp+tp-5.0);
        dphi(2,6) = 0.125*rp*sp*(rp+sp+2.0*tp-5.0);
        dphi(0,7) = -0.125*sp*tp*(2.0*rm+sp+tp-5.0);
        dphi(1,7) = 0.125*rm*tp*(rm+2.0*sp+tp-5.0);
        dphi(2,7) = 0.125*rm*sp*(rm+sp+2.0*tp-5.0);
        // centres in (t=-1)-plane
        dphi(0,8) = 0.5*(rm-1.0)*sm*tm;
        dphi(1,8) = -0.25*rrm*tm;
        dphi(2,8) = -0.25*rrm*sm;
        dphi(0,9) = 0.25*ssm*tm;
        dphi(1,9) = 0.5*rp*(sm-1.0)*tm;
        dphi(2,9) = -0.25*rp*ssm;
        dphi(0,10) = 0.5*(rm-1.0)*sp*tm;
        dphi(1,10) = 0.25*rrm*tm;
        dphi(2,10) = -0.25*rrm*sp;
        dphi(0,11) = -0.25*ssm*tm;
        dphi(1,11) = 0.5*rm*(sm-1.0)*tm;
        dphi(2,11) = -0.25*rm*ssm;
        // centres in (t=1)-plane
        dphi(0,12) = 0.5*(rm-1.0)*sm*tp;
        dphi(1,12) = -0.25*rrm*tp;
        dphi(2,12) = 0.25*rrm*sm;
        dphi(0,13) = 0.25*ssm*tp;
        dphi(1,13) = 0.5*rp*(sm-1.0)*tp;
        dphi(2,13) = 0.25*rp*ssm;
        dphi(0,14) = 0.5*(rm-1.0)*sp*tp;
        dphi(1,14) = 0.25*rrm*tp;
        dphi(2,14) = 0.25*rrm*sp;
        dphi(0,15) = -0.25*ssm*tp;
        dphi(1,15) = 0.5*rm*(sm-1.0)*tp;
        dphi(2,15) = 0.25*rm*ssm;
        // centres in (t=0)-plane
        dphi(0,16) = -0.25*sm*ttm;
        dphi(1,16) = -0.25*rm*ttm;
        dphi(2,16) = 0.5*rm*sm*(tm-1.0);
        dphi(0,17) = 0.25*sm*ttm;
        dphi(1,17) = -0.25*rp*ttm;
        dphi(2,17) = 0.5*rp*sm*(tm-1.0);
        dphi(0,18) = 0.25*sp*ttm;
        dphi(1,18) = 0.25*rp*ttm;
        dphi(2,18) = 0.5*rp*sp*(tm-1.0);
        dphi(0,19) = -0.25*sp*ttm;
        dphi(1,19) = 0.25*rm*ttm;
        dphi(2,19) = 0.5*rm*sp*(tm-1.0);

        // convert derivatives from [-1,1] to [0,1]
        dphi *= 2.;

    }

    template<>
    void Shapefun< corlib::HEXAHEDRON, 20>::giveInterpolationPoints( MatDimNFun & ipoints )
    {
        ipoints( 0, 0 )  = 0.;   ipoints( 1, 0 )  = 0.;   ipoints( 2, 0 )  = 0.;
        ipoints( 0, 1 )  = 1.;   ipoints( 1, 1 )  = 0.;   ipoints( 2, 1 )  = 0.;
        ipoints( 0, 2 )  = 1.;   ipoints( 1, 2 )  = 1.;   ipoints( 2, 2 )  = 0.;
        ipoints( 0, 3 )  = 0.;   ipoints( 1, 3 )  = 1.;   ipoints( 2, 3 )  = 0.;

        ipoints( 0, 4 )  = 0.;   ipoints( 1, 4 )  = 0.;   ipoints( 2, 4 )  = 1.;
        ipoints( 0, 5 )  = 1.;   ipoints( 1, 5 )  = 0.;   ipoints( 2, 5 )  = 1.;
        ipoints( 0, 6 )  = 1.;   ipoints( 1, 6 )  = 1.;   ipoints( 2, 6 )  = 1.;
        ipoints( 0, 7 )  = 0.;   ipoints( 1, 7 )  = 1.;   ipoints( 2, 7 )  = 1.;

        ipoints( 0, 8 )  = 0.5;  ipoints( 1, 8 )  = 0.;   ipoints( 2, 8 )  = 0.;
        ipoints( 0, 9 )  = 1.;   ipoints( 1, 9 )  = 0.5;  ipoints( 2, 9 )  = 0.;
        ipoints( 0, 10 ) = 0.5;  ipoints( 1, 10 ) = 1.;   ipoints( 2, 10 ) = 0.;
        ipoints( 0, 11 ) = 0.;   ipoints( 1, 11 ) = 0.5;  ipoints( 2, 11 ) = 0.;

        ipoints( 0, 12 ) = 0.5;  ipoints( 1, 12 ) = 0.;   ipoints( 2, 12 ) = 1.;
        ipoints( 0, 13 ) = 1.;   ipoints( 1, 13 ) = 0.5;  ipoints( 2, 13 ) = 1.;
        ipoints( 0, 14 ) = 0.5;  ipoints( 1, 14 ) = 1.;   ipoints( 2, 14 ) = 1.;
        ipoints( 0, 15 ) = 0.;   ipoints( 1, 15 ) = 0.5;  ipoints( 2, 15 ) = 1.;

        ipoints( 0, 16 ) = 0.;   ipoints( 1, 16 ) = 0.;   ipoints( 2, 16 ) = 0.5;
        ipoints( 0, 17 ) = 1.;   ipoints( 1, 17 ) = 0.;   ipoints( 2, 17 ) = 0.5;
        ipoints( 0, 18 ) = 1.;   ipoints( 1, 18 ) = 1.;   ipoints( 2, 18 ) = 0.5;
        ipoints( 0, 19 ) = 0.;   ipoints( 1, 19 ) = 1.;   ipoints( 2, 19 ) = 0.5;

        return;
    }

    //------------------------------------------------------------------------------
    /** \brief Tri-quadratic Lagrangian shape functions on hexahedron \f$[0,1]^3\f$
     *  \details functions  0..19 follow VTK_QUADRATIC_HEXAHEDRON (=25) ordering
     *           functions 20..26 similarly 
     *  <pre>
     *
     *     z
     *    ^
     *    |                                                                        
     *    |                                                                        
     *    4----15----7   
     *    |\         |\    
     *    |12    21  | 14  
     *   16  \ 25    19 \  
     *    |   5----13----6 
     *    |22 |  26  | 24| 
     *    0---|-11---3 --|----> y
     *     \ 17    23 \  18
     *      8 |  20    10| 
     *       \|         \| 
     *        1-----9----2 
     *         \
     *          \
     *          `'
     *            x
     *
     *  </pre>
     */   
    template<> 
    void Shapefun< corlib::HEXAHEDRON, 27 >::evaluate( const VecDim & xi,
                                                       VecNFun & phi ) const
    {
        // map from [0,1] to [-1,1]
        const double r = 2.0 * xi( 0 ) - 1.0;
        const double s = 2.0 * xi( 1 ) - 1.0;
        const double t = 2.0 * xi( 2 ) - 1.0;

        // auxiliary variables
        const double rm  = 1.0 - r;
        const double rp  = 1.0 + r;
        const double sm  = 1.0 - s;
        const double sp  = 1.0 + s;
        const double tm  = 1.0 - t;
        const double tp  = 1.0 + t;
        //const double rrm = 1.0 - 2.0*r;
        //const double rrp = 1.0 + 2.0*r;
        //const double ssm = 1.0 - 2.0*s;
        //const double ssp = 1.0 + 2.0*s;
        //const double ttm = 1.0 - 2.0*t;
        //const double ttp = 1.0 + 2.0*t;

        // corner nodes
        phi( 0) = -0.125* r*rm * s*sm * t*tm;
        phi( 1) = 0.125 * r*rp * s*sm * t*tm;
        phi( 2) = -0.125 * r*rp * s*sp * t*tm;
        phi( 3) = 0.125 * r*rm * s*sp * t*tm;
        phi( 4) = 0.125 * r*rm * s*sm * t*tp;
        phi( 5) = -0.125 * r*rp * s*sm * t*tp;
        phi( 6) = 0.125 * r*rp * s*sp * t*tp;
        phi( 7) = -0.125 * r*rm * s*sp * t*tp;
        // edge centre nodes
        phi( 8) = 0.25 * rm*rp * s*sm * t*tm;
        phi( 9) = -0.25 * r*rp * sm*sp * t*tm;
        phi(10) = -0.25 * rm*rp * s*sp * t*tm;
        phi(11) = 0.25 * r*rm * sm*sp * t*tm;
        phi(12) = -0.25 * rm*rp * s*sm * t*tp;
        phi(13) = 0.25 * r*rp * sm*sp * t*tp;
        phi(14) = 0.25 * rm*rp * s*sp * t*tp;
        phi(15) = -0.25 * r*rm * sm*sp * t*tp;
        phi(16) = 0.25 * r*rm * s*sm * tm*tp;
        phi(17) = -0.25 * r*rp * s*sm * tm*tp;
        phi(18) = 0.25 * r*rp * s*sp * tm*tp;
        phi(19) = -0.25 * r*rm * s*sp * tm*tp;
        // side centre nodes
        phi(20) = -0.5 * rm*rp * sm*sp * t*tm;
        phi(21) = 0.5 * rm*rp * sm*sp * t*tp;
        phi(22) = -0.5 * rm*rp * s*sm * tm*tp;
        phi(23) = 0.5 * r*rp * sm*sp * tm*tp;
        phi(24) = 0.5 * rm*rp * s*sp * tm*tp;
        phi(25) = -0.5 * r*rm * sm*sp * tm*tp;
        // volume centre node
        phi(26) = rm*rp * sm*sp * tm*tp;

    }

    template<> 
    void Shapefun< corlib::HEXAHEDRON, 27 >::evaluateGradient( const VecDim & xi,
                                                               MatDimNFun & dphi ) const
    {
        // map from [0,1] to [-1,1]
        const double r = 2.0 * xi( 0 ) - 1.0;
        const double s = 2.0 * xi( 1 ) - 1.0;
        const double t = 2.0 * xi( 2 ) - 1.0;

        // auxiliary variables
        const double rm  = 1.0 - r;
        const double rp  = 1.0 + r;
        const double sm  = 1.0 - s;
        const double sp  = 1.0 + s;
        const double tm  = 1.0 - t;
        const double tp  = 1.0 + t;
        const double rrm = 1.0 - 2.0*r;
        const double rrp = 1.0 + 2.0*r;
        const double ssm = 1.0 - 2.0*s;
        const double ssp = 1.0 + 2.0*s;
        const double ttm = 1.0 - 2.0*t;
        const double ttp = 1.0 + 2.0*t;

        dphi(0,0) = -0.125 * rrm * s*sm * t*tm;  // d/dr
        dphi(1,0) = -0.125 * r*rm * ssm * t*tm;  // d/ds
        dphi(2,0) = -0.125 * r*rm * s*sm * ttm;  // d/dt

        dphi(0,1) = 0.125 * rrp * s*sm * t*tm;
        dphi(1,1) = 0.125 * r*rp * ssm * t*tm;
        dphi(2,1) = 0.125 * r*rp * s*sm * ttm;

        dphi(0,2) = -0.125 * rrp * s*sp * t*tm;
        dphi(1,2) = -0.125 * r*rp * ssp * t*tm;
        dphi(2,2) = -0.125 * r*rp * s*sp * ttm;

        dphi(0,3) = 0.125 * rrm * s*sp * t*tm;
        dphi(1,3) = 0.125 * r*rm * ssp * t*tm;
        dphi(2,3) = 0.125 * r*rm * s*sp * ttm;

        dphi(0,4) = 0.125 * rrm * s*sm * t*tp;
        dphi(1,4) = 0.125 * r*rm * ssm * t*tp;
        dphi(2,4) = 0.125 * r*rm * s*sm * ttp;

        dphi(0,5) = -0.125 * rrp * s*sm * t*tp;
        dphi(1,5) = -0.125 * r*rp * ssm * t*tp;
        dphi(2,5) = -0.125 * r*rp * s*sm * ttp;

        dphi(0,6) = 0.125 * rrp * s*sp * t*tp;
        dphi(1,6) = 0.125 * r*rp * ssp * t*tp;
        dphi(2,6) = 0.125 * r*rp * s*sp * ttp;

        dphi(0,7) = -0.125 * rrm * s*sp * t*tp;
        dphi(1,7) = -0.125 * r*rm * ssp * t*tp;
        dphi(2,7) = -0.125 * r*rm * s*sp * ttp;

        dphi(0,8) = -0.5 * r * s*sm * t*tm;
        dphi(1,8) = 0.25 * rm*rp * ssm * t*tm;
        dphi(2,8) = 0.25 * rm*rp * s*sm * ttm;

        dphi(0,9) = -0.25 * rrp * sm*sp * t*tm;
        dphi(1,9) = 0.5 * r*rp * s * t*tm;
        dphi(2,9) = -0.25 * r*rp * sm*sp * ttm;

        dphi(0,10) = 0.5 * r * s*sp * t*tm;
        dphi(1,10) = -0.25 * rm*rp * ssp * t*tm;
        dphi(2,10) = -0.25 * rm*rp * s*sp * ttm;

        dphi(0,11) = 0.25 * rrm * sm*sp * t*tm;
        dphi(1,11) = -0.5 * r*rm * s * t*tm;
        dphi(2,11) = 0.25 * r*rm * sm*sp * ttm;

        dphi(0,16) = 0.25 * rrm * s*sm * tm*tp;
        dphi(1,16) = 0.25 * r*rm * ssm * tm*tp;
        dphi(2,16) = -0.5 * r*rm * s*sm * t;

        dphi(0,17) = -0.25 * rrp * s*sm * tm*tp;
        dphi(1,17) = -0.25 * r*rp * ssm * tm*tp;
        dphi(2,17) = 0.5 * r*rp * s*sm * t;

        dphi(0,18) = 0.25 * rrp * s*sp * tm*tp;
        dphi(1,18) = 0.25 * r*rp * ssp * tm*tp;
        dphi(2,18) = -0.5 * r*rp * s*sp * t;

        dphi(0,19) = -0.25 * rrm * s*sp * tm*tp;
        dphi(1,19) = -0.25 * r*rm * ssp * tm*tp;
        dphi(2,19) = 0.5 * r*rm * s*sp * t;

        dphi(0,12) = 0.5 * r * s*sm * t*tp;
        dphi(1,12) = -0.25 * rm*rp * ssm * t*tp;
        dphi(2,12) = -0.25 * rm*rp * s*sm * ttp;

        dphi(0,13) = 0.25 * rrp * sm*sp * t*tp;
        dphi(1,13) = -0.5 * r*rp * s * t*tp;
        dphi(2,13) = 0.25 * r*rp * sm*sp * ttp;

        dphi(0,14) = -0.5 * r * s*sp * t*tp;
        dphi(1,14) = 0.25 * rm*rp * ssp * t*tp;
        dphi(2,14) = 0.25 * rm*rp * s*sp * ttp;

        dphi(0,15) = -0.25 * rrm * sm*sp * t*tp;
        dphi(1,15) = 0.5 * r*rm * s * t*tp;
        dphi(2,15) = -0.25 * r*rm * sm*sp * ttp;

        dphi(0,20) = r * sm*sp * t*tm;
        dphi(1,20) = rm*rp * s * t*tm;
        dphi(2,20) = -0.5 * rm*rp * sm*sp * ttm;

        dphi(0,21) = -r * sm*sp * t*tp;
        dphi(1,21) = -rm*rp * s * t*tp;
        dphi(2,21) = 0.5 * rm*rp * sm*sp * ttp;

        dphi(0,22) = r * s*sm * tm*tp;
        dphi(1,22) = -0.5 * rm*rp * ssm * tm*tp;
        dphi(2,22) = rm*rp * s*sm * t;

        dphi(0,23) = 0.5 * rrp * sm*sp * tm*tp;
        dphi(1,23) = -r*rp * s * tm*tp;
        dphi(2,23) = -r*rp * sm*sp * t;

        dphi(0,24) = -r * s*sp * tm*tp;
        dphi(1,24) = 0.5 * rm*rp * ssp * tm*tp;
        dphi(2,24) = -rm*rp * s*sp * t;

        dphi(0,25) = -0.5 * rrm * sm*sp * tm*tp;
        dphi(1,25) = r*rm * s * tm*tp;
        dphi(2,25) = r*rm * sm*sp * t;

        dphi(0,26) = -2.0 * r * sm*sp * tm*tp;
        dphi(1,26) = -2.0 * rm*rp * s * tm*tp;
        dphi(2,26) = -2.0 * rm*rp * sm*sp * t;

        // convert derivatives from [-1,1] to [0,1]
        dphi *= 2.;

    }

    template<>
    void Shapefun< corlib::HEXAHEDRON, 27>::giveInterpolationPoints( MatDimNFun & ipoints )
    {
        ipoints( 0, 0 )  = 0.;   ipoints( 1, 0 )  = 0.;   ipoints( 2, 0 )  = 0.;
        ipoints( 0, 1 )  = 1.;   ipoints( 1, 1 )  = 0.;   ipoints( 2, 1 )  = 0.;
        ipoints( 0, 2 )  = 1.;   ipoints( 1, 2 )  = 1.;   ipoints( 2, 2 )  = 0.;
        ipoints( 0, 3 )  = 0.;   ipoints( 1, 3 )  = 1.;   ipoints( 2, 3 )  = 0.;

        ipoints( 0, 4 )  = 0.;   ipoints( 1, 4 )  = 0.;   ipoints( 2, 4 )  = 1.;
        ipoints( 0, 5 )  = 1.;   ipoints( 1, 5 )  = 0.;   ipoints( 2, 5 )  = 1.;
        ipoints( 0, 6 )  = 1.;   ipoints( 1, 6 )  = 1.;   ipoints( 2, 6 )  = 1.;
        ipoints( 0, 7 )  = 0.;   ipoints( 1, 7 )  = 1.;   ipoints( 2, 7 )  = 1.;

        ipoints( 0, 8 )  = 0.5;  ipoints( 1, 8 )  = 0.;   ipoints( 2, 8 )  = 0.;
        ipoints( 0, 9 )  = 1.;   ipoints( 1, 9 )  = 0.5;  ipoints( 2, 9 )  = 0.;
        ipoints( 0, 10 ) = 0.5;  ipoints( 1, 10 ) = 1.;   ipoints( 2, 10 ) = 0.;
        ipoints( 0, 11 ) = 0.;   ipoints( 1, 11 ) = 0.5;  ipoints( 2, 11 ) = 0.;

        ipoints( 0, 12 ) = 0.5;  ipoints( 1, 12 ) = 0.;   ipoints( 2, 12 ) = 1.;
        ipoints( 0, 13 ) = 1.;   ipoints( 1, 13 ) = 0.5;  ipoints( 2, 13 ) = 1.;
        ipoints( 0, 14 ) = 0.5;  ipoints( 1, 14 ) = 1.;   ipoints( 2, 14 ) = 1.;
        ipoints( 0, 15 ) = 0.;   ipoints( 1, 15 ) = 0.5;  ipoints( 2, 15 ) = 1.;

        ipoints( 0, 16 ) = 0.;   ipoints( 1, 16 ) = 0.;   ipoints( 2, 16 ) = 0.5;
        ipoints( 0, 17 ) = 1.;   ipoints( 1, 17 ) = 0.;   ipoints( 2, 17 ) = 0.5;
        ipoints( 0, 18 ) = 1.;   ipoints( 1, 18 ) = 1.;   ipoints( 2, 18 ) = 0.5;
        ipoints( 0, 19 ) = 0.;   ipoints( 1, 19 ) = 1.;   ipoints( 2, 19 ) = 0.5;

        ipoints( 0, 20 ) = 0.5;  ipoints( 1, 20 ) = 0.5;  ipoints( 2, 20 ) = 0.;
        ipoints( 0, 21 ) = 0.5;  ipoints( 1, 21 ) = 0.5;  ipoints( 2, 21 ) = 1.;

        ipoints( 0, 22 ) = 0.5;  ipoints( 1, 22 ) = 0.;   ipoints( 2, 22 ) = 0.5;
        ipoints( 0, 23 ) = 1.;   ipoints( 1, 23 ) = 0.5;  ipoints( 2, 23 ) = 0.5;
        ipoints( 0, 24 ) = 0.5;  ipoints( 1, 24 ) = 1.;   ipoints( 2, 24 ) = 0.5;
        ipoints( 0, 25 ) = 0.;   ipoints( 1, 25 ) = 0.5;  ipoints( 2, 25 ) = 0.5;

        ipoints( 0, 26 ) = 0.5;  ipoints( 1, 26 ) = 0.5;  ipoints( 2, 26 ) = 0.5;
        return;
    }

}// end namespace corlib
