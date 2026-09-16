// 
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// 
//                              KosalaBandara                                   
//                 Computational Structural Mechanics Lab                       
//                         University of Cambridge                              
// 
// This software is copyrighted and all rights are retained by the CSMLab. All 
// use, disclosure, and/or reproduction of any part not expressly authorized by 
// F Cirak is prohibited. (C) 2011. 
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// 
using namespace subdiv::curve;

//! Local index of next vertex
template <typename V, typename E, typename ET>
const std::array<int,ShapeDim<ET::myShape>::numVertices>
Mesh<V, E, ET>::next_ = ShapeProp<ET::myShape>::next;

//! Local index of previous vertex
template <typename V, typename E, typename ET>
const std::array<int,ShapeDim<ET::myShape>::numVertices>
Mesh<V, E, ET>::previous_ = ShapeProp<ET::myShape>::previous;

template <typename V, typename E, typename ET>
Mesh<V, E, ET>::Mesh(std::istream & smf):
    refinementLevel_(0)
{    
    FTL_VERIFY(Edge::myShape == ETree::myShape);
    
    // validate input file
    corlib::SmfHead smfHead;
    smfHead.readValidated( smf, myShape, numVertices );
    
    unsigned vsize = 0, esize = 0;
    smf >> vsize >> esize;
    
    vertices_.reserve(vsize);
    edges_.reserve(esize);
    
    //add vertices
    for(unsigned i = 0; i < vsize; ++i ){
	Vertex* v = new Vertex();
	smf >> v;
	this->addVertex_(v);
        v->setIndex(i);//TEMP
    }
    
    // add elements
    for(unsigned i = 0; i < esize; ++i ){
        VecVPtrNV vec;
        //collect vertices
        for (unsigned j = 0; j < numVertices; j++){
            unsigned id;
            smf >> id;
            vec[j]= vertices_[id];
        }

	Edge* e = new Edge(vec);
	this->addEdge_(e);
        e->setIndex(i);//TEMP
    }  
}

//! Destructor
template <typename V, typename E, typename ET>
Mesh<V, E, ET>::~Mesh()
{    
    //delete edge trees
    this->iterateOverEdgeTrees( corlib::deleteFunctor( ) );
    edgeTrees_.clear( );

    //delete edges
    this -> iterateOverEdges( corlib::deleteFunctor( ) );

    //delete vertices
    this->iterateOverVertices( corlib::deleteFunctor( ) );
    vertices_.clear( );
}


//! Iterate over vertices
template <typename V, typename E, typename ET>
template< typename OP >
OP Mesh<V, E, ET>::iterateOverVertices( OP op )
{
    VertexIterator begin = vertices_.begin( );
    VertexIterator end = vertices_.end( );
    return std::for_each( begin, end, op );
}

//! Iterate over edges
template <typename V, typename E, typename ET>
template< typename OP >
OP Mesh<V, E, ET>::iterateOverEdges( OP op )
{
    EdgeIterator begin = edges_.begin( );
    EdgeIterator end = edges_.end( );
    return std::for_each( begin, end, op );
}

//! Iterate over edgeTrees
template <typename V, typename E, typename ET>
template< typename OP >
OP Mesh<V, E, ET>::iterateOverEdgeTrees( OP op )
{
    EdgeIterator  begin=edges_.begin(), end=edges_.end();
    for ( ; begin!=end; ++begin)
        op(edgeTrees_[*begin]);
    return op;
}


//! Build the tree structure
template <typename V, typename E, typename ET>
void Mesh<V, E, ET>::buildTopology(void ) 
{
    //add all edges to a map
    typedef std::map<Vertex*, Edge* >                       VEMap;
    VEMap vertexEdgeMap;
    std::pair<typename VEMap::iterator, bool> ret;

    EdgeIterator  begin=edges_.begin(), end=edges_.end();
    Edge *e;

    for ( ; begin!=end; ++begin){
	e = *begin;
        for (unsigned i =0; i < numVertices; i++){
            ret = vertexEdgeMap.insert(typename VEMap::value_type(e->vertex(i), e) );

            // if vertex appear in the map
            if (!ret.second){
                // Given two vertices that are the same, set up the edge adjacency
                Vertex* vo = ret.first->first;
                Edge* eo = ret.first->second;
	
                FTL_VERIFY(vo != NULL);
                FTL_VERIFY(eo != NULL);
	
                // set neighbor information
                FTL_VERIFY(e!=eo);            
                eo->opposite(vo) = e;
                e->opposite(vo) = eo;
            }
        }
    }

    // instantiate Edge Trees
    begin=edges_.begin();
    for ( ; begin!=end; ++begin){
        VecVPtrNV vec;
        for (unsigned i=0; i<numVertices; i++)
            vec[i]=(*begin)->vertex(i);

        ETree* etree = new ETree(vec, 0, NULL);
        std::pair<typename ETreeMap::iterator, bool> 
            result = edgeTrees_.insert(typename ETreeMap::value_type (*begin, etree));
        assert(result.second);
    }
    
    // set neighbor information
    begin=edges_.begin();
    for ( ; begin!=end; ++begin){
        ETree *e = edgeTrees_[*begin];
        for (unsigned i = 0; i < numNeighbors; i++)
            e->neighbor(i) = edgeTrees_[(*begin)->neighbor(i)];
    }

    // mark border edge vertices as nosub-vertex
    begin=edges_.begin();
    for ( ; begin!=end; ++begin){
        for (unsigned i = 0; i < numNeighbors; ++i) {
            if ((*begin)->neighbor(i)==NULL){
                (*begin)->vertex(i)->setTag(subdiv::curve::NOSUB_VERTEX, 0);
                (*begin)->vertex(next_[i])->setTag(subdiv::curve::NEXT_NOSUB_VERTEX, 0);
            }
        }
    }

    return;
}

//------------------------------------------------------------------------------
/// Read tag file
template <typename V, typename E, typename ET>
void Mesh<V, E, ET>::readTags( std::istream & tg )
{
    typedef typename Vertex::VertexTag     VertexTag;

    std::string name;    
    tg >> name;

    FTL_VERIFY( name == "TG" );

    unsigned numVertexTags, numEdgeTags;
    tg >> numVertexTags >> numEdgeTags;

    FTL_VERIFY(numEdgeTags == 0);

    // set vertex tags
    for ( unsigned vt = 0; vt < numVertexTags; vt ++ ) {
        unsigned v;
        int tagInt;
        tg >> v >> tagInt;

        const VertexTag tag = static_cast< VertexTag >( tagInt );
        vertices_[ v ]->setTag( tag );
    }

    return;
}

//! Subdivide mesh
template <typename V, typename E, typename ET>
template< typename SUBDIV >
void Mesh<V, E, ET>::subdivide( SUBDIV * subdiv  ) 
{
    refinementLevel_++;
    
    // subdivide all elements
    typedef std::function<void (ETree*, int, SUBDIV*) > SubdivFn;
    SubdivFn  subdivFn = &ETree::template subdivide<SUBDIV>;
    this -> iterateOverEdgeTrees( std::bind(subdivFn, std::placeholders::_1, refinementLevel_, subdiv));

    return;
}


//! Write mesh to output file
template <typename V, typename E, typename ET>
void Mesh<V, E, ET>::writeMesh(std::ostream & smf,std::ostream * tg,
                               const NamedAddCoord addCoord, std::ostream * vdat )

{    
    VecVPtr vVec;
    VecETreePtr eVec;
    this->collectETreesVertices_(vVec, eVec);
    unsigned numVerticesMesh = vVec.size();

    // write header
    corlib::SmfHead smfHead;
    smfHead.setElementShape( myShape );
    smfHead.setElementNumPoints( numVertices );
    smfHead.write( smf );

    smf << vVec.size() << " " << eVec.size() << std::endl;

    VertexIterator vBegin = vVec.begin(), vEnd = vVec.end();
    unsigned number=0;
    for(; vBegin != vEnd; ++vBegin){
        (*vBegin)->setIndex(number++);
        smf << *vBegin  <<std::endl;
    }

    EdgeTreeVecIterator aBegin = eVec.begin(), aEnd=eVec.end();    
    for(; aBegin!= aEnd; ++aBegin){
        for (unsigned i = 0; i < numVertices; i++)
            smf << (*aBegin)->vertex(i)->index() << " ";
        smf << std::endl;
    }

    // write vertex data on additional co-ordinates
    if ( not addCoord.empty() and vdat ) {
        // write header
        *vdat <<  numVerticesMesh << " " << addCoord.size() << std::endl;
        unsigned numTotalDataPerVertex = 0;
        for ( NamedAddCoord::const_iterator ac = addCoord.begin(); ac != addCoord.end(); ++ac ) {
            *vdat << ac->first << " " << ac->second << std::endl;
            numTotalDataPerVertex += ac->second;
        }
        FTL_VERIFY( numTotalDataPerVertex <= ( Vertex::numVariables - Vertex::dim ) );

        // write coord
        unsigned offset = Vertex::dim;
        for ( NamedAddCoord::const_iterator ac = addCoord.begin(); ac != addCoord.end(); ++ac ) {
            const unsigned numDataPerVertex = ac->second;
            vBegin = vVec.begin();
            for(; vBegin != vEnd; ++vBegin){
                const typename Vertex::VecP & coord = (*vBegin)->point();
                for ( unsigned a = 0; a < numDataPerVertex; ++a ) {
                    *vdat << coord[ offset + a ] << " ";
                }
                *vdat << std::endl;
            }
            // update
            offset += numDataPerVertex;
        }
    }


    return;
}

/// Collect all vertices
template <typename V, typename E, typename ET>
void Mesh<V, E, ET>::collectETreesVertices_(VecVPtr & vVec, VecETreePtr & eVec)
{
    vVec.clear();
    eVec.clear();

    //collect all edgeTrees in a vector
    EdgeIterator  eBegin=edges_.begin(), eEnd=edges_.end();
    for ( ; eBegin != eEnd; ++eBegin)
        edgeTrees_[*eBegin]->addLeafs(std::back_inserter(eVec)); 
    
    //add existing vertices
    vVec.insert(vVec.begin(), this->vBegin(), this->vEnd());

    //add new vertices 
    eBegin = this->eBegin();
    eEnd   = this->eEnd();
    for ( ; eBegin != eEnd; ++eBegin )
        this->edgeTree( *eBegin )->addNewVertices( std::back_inserter( vVec ) );
}


//! Update edge tree and vertex indices
template <typename V, typename E, typename ET>
void Mesh<V, E, ET>::updateIndices_()
{    
    VecVPtr vVec;
    VecETreePtr eVec;
    this->collectETreesVertices_(vVec, eVec);

    //update vertex indices
    VertexIterator vBegin = vVec.begin(), vEnd = vVec.end();
    unsigned number=0;
    for(; vBegin != vEnd; ++vBegin)
        (*vBegin)->setIndex(number++);
    
    // update edge tree indices
    EdgeTreeVecIterator aBegin = eVec.begin(), aEnd = eVec.end();    
    number = 0;
    for(; aBegin!= aEnd; ++aBegin)
        (*aBegin)->setIndex(number++);

    return;
}


