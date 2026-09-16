//
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//
//                 Computational Structural Mechanics Lab
//                         University of Cambridge
//
// This software is copyrighted and all rights are retained by the CSMLab. All
// use, disclosure, and/or reproduction of any part not expressly authorized by
// F Cirak is prohibited. (C) 2011.
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//

//! @author Kosala Bandara, Burkhard Bornemann
//! @date   2011


//------------------------------------------------------------------------------
template< typename FACET, typename SUBDIV >
subdiv::twoD::ShapeFunSubdivision< FACET, SUBDIV >::ShapeFunSubdivision(
    Edge * f, Subdiv * subdiv
    ) : facet_( f ), subdiv_( subdiv )
{


    return;
}

//------------------------------------------------------------------------------
template< typename FACET, typename SUBDIV >
void subdiv::twoD::ShapeFunSubdivision< FACET, SUBDIV >::evaluate(
    const VecLDim & xi,
    Vec & phi)
{
    // determine subdivison matrix
    Mat subMat( numFunctions, patchVertices_.size( ) );
    VecLDim subXi;
    MatLDimLDim subJac;
    this->computeSubdivisionMatrix_( xi, subMat, subXi, subJac );

    // determine shape functions over patch
    VecNF subPhi;
    spline_.evaluate( subXi, subPhi );
    phi = subPhi * subMat;

    return;
}

//------------------------------------------------------------------------------
template< typename FACET, typename SUBDIV >
void subdiv::twoD::ShapeFunSubdivision< FACET, SUBDIV >::evaluateGradient(
    const VecLDim & xi,
    Mat & dPhiDXi
    )
{
    // determine subdivison matrix
    Mat subMat;
    VecLDim subXi;
    MatLDimLDim subJac;
    this->computeSubdivisionMatrix_( xi, subMat, subXi, subJac );
    
    // determine 1st derivative shape functions over patch
    MatLDimNF subDPhiDXi;
    spline_.evaluateGradient( subXi, subDPhiDXi );
    dPhiDXi = ( subJac.transpose( ) * subDPhiDXi ) * subMat;

    return;
}

//------------------------------------------------------------------------------
template< typename FACET, typename SUBDIV >
void subdiv::twoD::ShapeFunSubdivision< FACET, SUBDIV >::evaluateHessian(
    const VecLDim & xi,
    Mat & ddPhiDDXi)
{
    // determine subdivison matrix
    Mat subMat;
    VecLDim subXi;
    MatLDimLDim subJac;
    this->computeSubdivisionMatrix_( xi, subMat, subXi, subJac );
    
    // determine 2nd derivative of shape functions over patch
    MatVecNFLDimLDim subDDPhiDDXi;
    spline_.evaluateHessian( subXi, subDDPhiDDXi );
    MatSDimNF temp2;
    for ( unsigned f = 0; f < numFunctions; ++f ) {
        MatLDimLDim ddPhiDDXi_f;
        for ( unsigned i = 0; i < localDim; ++i )
            for ( unsigned j = 0; j < localDim; ++j )
                ddPhiDDXi_f( i, j ) = subDDPhiDDXi( i, j )( f );
        ddPhiDDXi_f = subJac.transpose( ) * ( ddPhiDDXi_f * subJac );
        // store
        temp2( 0, f ) = ddPhiDDXi_f( 0, 0 );
        temp2( 1, f ) = ddPhiDDXi_f( 0, 1 );
        temp2( 2, f ) = ddPhiDDXi_f( 1, 1 );
    }
    ddPhiDDXi = temp2 * subMat;

    return;
}

//------------------------------------------------------------------------------
template< typename FACET, typename SUBDIV >
void subdiv::twoD::ShapeFunSubdivision< FACET, SUBDIV >::evaluateGradHess(
    const VecLDim & xi,
    Vec & phi,
    Mat & dPhiDXi,
    Mat & ddPhiDDXi
    )
{
    // determine subdivison matrix
    Mat subMat;
    VecLDim subXi;
    MatLDimLDim subJac;
    this->computeSubdivisionMatrix_( xi, subMat, subXi, subJac );

    // determine shape functions over patch
    VecNF subPhi;
    spline_.evaluate( subXi, subPhi );
    phi = subPhi * subMat;

    // determine 1st derivative of shape functions over patch
    MatLDimNF subDPhiDXi;
    spline_.evaluateGradient( subXi, subDPhiDXi );
    dPhiDXi = ( subJac.transpose( ) * subDPhiDXi ) * subMat;

    // determine 2nd derivative of shape functions over patch
    MatVecNFLDimLDim subDDPhiDDXi;
    spline_.evaluateHessian( subXi, subDDPhiDDXi );
    MatSDimNF temp2;
    for ( unsigned f = 0; f < numFunctions; ++f ) {
        MatLDimLDim ddPhiDDXi_f;
        for ( unsigned i = 0; i < localDim; ++i )
            for ( unsigned j = 0; j < localDim; ++j )
                ddPhiDDXi_f( i, j ) = subDDPhiDDXi( i, j )( f );
        ddPhiDDXi_f = subJac.transpose( ) * ( ddPhiDDXi_f * subJac );
        // store
        temp2( 0, f ) = ddPhiDDXi_f( 0, 0 );
        temp2( 1, f ) = ddPhiDDXi_f( 0, 1 );
        temp2( 2, f ) = ddPhiDDXi_f( 1, 1 );
    }
    ddPhiDDXi = temp2 * subMat;

    return;
}


//------------------------------------------------------------------------------
template< typename FACET, typename SUBDIV >
void subdiv::twoD::ShapeFunSubdivision< FACET, SUBDIV >::formPatchMesh_(
    Mesh_ & patch
    )
{
    // generate patch streams
    std::stringstream smf, tg;
    this->writePatch_( smf, tg );  // this contains #patchEdges_ and #extraEdges_
    smf.seekg( 0 );
    tg.seekg( 0 );

    // generate patch mesh
    patch.readSmf( smf );
    patch.readTags( tg );
    patch.buildEdgeTopology( );
    //patch.iterateOverEdges( boost::bind( &Edge::write, _1, boost::ref( std::cout ) ) );
    patch.implantEdgeTrees( );
    //patch.iterateOverEdgeTrees( boost::bind( &EdgeTree_::write, _1, boost::ref( std::cout ) ) );

    return;
}

//------------------------------------------------------------------------------
template< typename FACET, typename SUBDIV >
void subdiv::twoD::ShapeFunSubdivision< FACET, SUBDIV >::computeSubdivisionMatrix_(
    const VecLDim & xi,
    Mat & subMat,
    VecLDim & subXi,
    MatLDimLDim & subJac )
{
    double xiC = detail_::criticalDimension<myShape>(xi);    
    const unsigned n = static_cast<unsigned>(ceil(log2(1./xiC)));

    //find kx, ky and compute new local coordinates
    const double jfac = pow(2.0, n);
    double kxD;

    xiLocal[0]= modf(xi[0] * jfac, &kxD);

    // cast integer parts to unsigned ints
    unsigned kx = static_cast<unsigned>(kxD);

    //determine the jacobians of the linear  mapping
    MatLDimLDim_ j1_deriv;         j1_deriv.clear(); 
    MatDeriv2Deriv2_ j1_hessian;   j1_hessian.clear();    

    for(unsigned i = 0; i < localDim; i++)
        j1_deriv(i,i)=1/jfac;     
    
    for(unsigned i = 0; i < numDeriv2; i++)
        j1_hessian(i,i)=1/(jfac*jfac);   

    corlib::inverse(j1_deriv);
    corlib::inverse(j1_hessian);

    j_deriv = j1_deriv;
    j_hessian = j1_hessian;

    // compute the subdivision matrix
    unsigned K = subMat.size2();
    unsigned numAdVariables = Vertex::numVariables - Vertex::dim;
    //iterate K times and calculate each row of the subdivision matrix
    for (unsigned i_K = 0; i_K < K; i_K += numAdVariables){
	// reset the mesh each time 
	this->resetMesh();
	
        //set v[i_K] to one and others to zero
	for (unsigned i = 0; i < vertices_.size(); i++){

            if (i_K <= i and  i < i_K + numAdVariables){
                //find which additional variable to set to 1
                unsigned i_K2 = i - i_K;
                //check if all rows have already been assigned
                VecP vec = vertices_[i]->point();

                if (i_K + i_K2 < K)
                    vec[dim + i_K2] = 1.;
                else
                    vec[dim + i_K2] = 0.;

                vertices_[i]->setPoint(vec);
                //vertices_[i]->point() = vec;
                // std::cout<<vec<<" |"<<vertices_[i]->point()<<std::endl;//TEMP
            }
        }
	
	//obtain the reference vertex and reference element
        ETree* refETree = edgeTrees_[0];
        Vertex* refVertex = refETree->vertex(0);

        //check for extrodinary vertex
        unsigned locXvtx = 0;
        for (unsigned k = 0; k < numNeighbors; k++){
            if (refETree->neighbor(k)==NULL)
                locXvtx = k;
        }

        //adjust kx for special case
        if (locXvtx == 1 and
            corlib::fuzzyEqual(xiLocal[0], 0.) and kx >0){
            kx--;
            xiLocal[0]=1.0;
        }


        unsigned id = refVertex->index();//TEMP
        unsigned id1 = refETree->index();//TEMP
        //std::cout<<refVertex->point()<<std::endl;

	//iterate n times and  subdivide
        for (unsigned i = 0; i < edgeTrees_.size(); i++) {
            //edgeTrees_[i]->subdivide(n, true); //needed if writeMesh is used
            edgeTrees_[i]->subdivide(n); 
	}
        
        
        //update vertex labels
        

        //std::cout<<refETree->isLeaf()<<std::endl;
        // std::ofstream temp ("output/patch.smf");//TEMP
        // this->writeMesh(temp);

	// //build support vertices
        // for (unsigned i = 0; i < edgeTrees_.size(); i++) {
        //     edgeTrees_[i]->buildSupportVertices(); 
	// }

        for (unsigned i = 0 ; i < n ; i++)
            refETree = refETree->child(0);

        id = refETree->index();//TEMP

	//extract correct regular patch
        for (unsigned i = 0; i < kx; i++){
            int j = localVtxNumber(refETree, refVertex);
            refVertex = refETree->vertex(next_[j]);
            refETree = refETree->neighbor(next_[j]);
        }

        id = refETree->index();//TEMP

	//add data to Subdivision Matrix
        refETree->buildSupportVertices();
        VertexVec vReg = refETree->regularPatch();
	for (unsigned i = 0; i < vReg.size(); i++){
            VecP vec = vReg[i]->point(n);
            //std::cout<<vec<<std::endl;//TEMP
            for (unsigned j = 0; j < numAdVariables; j++)
                if (i_K+j < K)
                    subMat(i, i_K+j)= vec(dim+j);
	}
    }

    return;
}
