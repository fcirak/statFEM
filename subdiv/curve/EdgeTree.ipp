using namespace subdiv::curve;

/// Local index of next vertex
template < corlib::shape SHAPE, typename V>
const std::array<int, subdiv::curve::ShapeDim<SHAPE>::numVertices>
EdgeTree<SHAPE, V>::next_ = ShapeProp<SHAPE>::next;

/// Local index of previous vertex
template < corlib::shape SHAPE, typename V>
const std::array<int, subdiv::curve::ShapeDim<SHAPE>::numVertices>
EdgeTree<SHAPE, V>::previous_ = ShapeProp<SHAPE>::previous;


/// Constructor
template < corlib::shape SHAPE, typename V>
EdgeTree<SHAPE, V>::EdgeTree(VecVPtrNV vec, int depth, ETree* parent) 
    : vertices_(vec), parent_(parent), depth_(depth)
{
    for (unsigned i = 0; i < numChildren; i++)
        children_[i] = NULL;
    
    for (unsigned i = 0; i < numNeighbors; i++)
        neighbors_[i] = NULL;

    //exSupport_ = NULL;
    index_ = 0;
    newVertices_.clear();
}	


/// Destructor
template < corlib::shape SHAPE, typename V>
EdgeTree<SHAPE, V>::~EdgeTree()
{
    //delete children
    if ( not this->isLeaf( ) ) {        
        std::for_each( children_.begin( ), children_.end( ),
                       corlib::deleteFunctor( ) );
    }

    // delete vertices
    std::for_each( newVertices_.begin( ), newVertices_.end( ),
                   corlib::deleteFunctor( ) );
    newVertices_.clear( );
}

template < corlib::shape SHAPE, typename V>
unsigned EdgeTree<SHAPE, V>::localVertexIndex(const Vertex* v) const
{
    typename VecVPtrNV::const_iterator vIter =
        std::find( vertices_.begin(), vertices_.end(), v );
    FTL_VERIFY_DESCRIPTIVE( vIter != vertices_.end( ),
                            "Vertex is not in facet of facettree\n" );
    return std::distance( vertices_.begin(), vIter );
}

/// Get neighbor opposite to vertex v
template < corlib::shape SHAPE, typename V>
EdgeTree<SHAPE, V>*&
EdgeTree<SHAPE, V>::opposite(const Vertex* v) 
{
    FTL_VERIFY(v != NULL);
    
    unsigned indexOpposite = numNeighbors;
    for (unsigned i = 0; i < numVertices; i++){
        if (v == vertices_[i]){
            indexOpposite = i;
            break;
        }
    }
    FTL_VERIFY(indexOpposite!= numNeighbors);
    return neighbors_[indexOpposite];
}

/// Check if is a leaf (has no children)
template < corlib::shape SHAPE, typename V>
bool EdgeTree<SHAPE, V>::isLeaf() 
{
    bool leaf = true;
    for (unsigned i = 0; i < numChildren; i++)
        if (children_[i]){
            leaf = false;
            break;
        }

    return leaf;
}

/// Subdivide edge tree up to a specific level
template < corlib::shape SHAPE, typename V>
template< typename SUBDIV >
void EdgeTree<SHAPE, V>::subdivide(const int level, SUBDIV * subdiv) 
{
    FTL_VERIFY(level > depth_);
    
    // check if facet tree leaf
    if ( (this->isLeaf()) and 
         (level == depth_+1)) {
        // subdivide the leaf (no children)
        this->subdivideLeaf_( subdiv );
    }

    else if (level > depth_+1){
        //subdivide each children
        for (unsigned i = 0; i < numChildren; i++)
            children_[i]->subdivide(level, subdiv);
    }
    
    return;
}


/// Find extordinary vertex if any
template < corlib::shape SHAPE, typename V>
V* EdgeTree<SHAPE, V>::extVertex()
{
    Vertex* v = NULL;
    for (unsigned i = 0; numNeighbors; i++){
        if (neighbor(i)==NULL)
            v = vertex(i);
    }
    return v;
}

/// Collect active leafs to a vector
template < corlib::shape SHAPE, typename V>
template <typename OutputIterator>
void EdgeTree<SHAPE, V>::addLeafs(OutputIterator iter) 
{
    if ( this->isLeaf() ) {
        *(iter++) = this;
    }
    else {
        for ( unsigned i = 0; i < numChildren; i++ )
            children_[ i ]->addLeafs( iter );
    }

    return;
}

/// Collect all leafs to a vector
template < corlib::shape SHAPE, typename V>
template <typename OutputIterator>
void EdgeTree<SHAPE, V>::addLeafsAll(OutputIterator iter) 
{
    *(iter++) = this;

    if (!isLeaf()){
        for (unsigned i = 0; i < numChildren; i++)
            children_[i]->addLeafsAll(iter);
    }
 
    return;
}

/// Collect all new vertices up to a particular level
template < corlib::shape SHAPE, typename V>
template <typename OutputIterator>
void EdgeTree<SHAPE, V>::addNewVertices(const int level, OutputIterator iter) 
{

    if (level > depth_)
        std::copy( newVertices_.begin( ), newVertices_.end( ), iter );

    if ( level > depth_+1 ) {
        FTL_VERIFY(not this->isLeaf());
        for ( unsigned i = 0; i < numChildren; i++ )
            children_[ i ]->addNewVertices(level,  iter );
    }
    return;
}

/// Collect all new vertices
template < corlib::shape SHAPE, typename V>
template <typename OutputIterator>
void EdgeTree<SHAPE, V>::addNewVertices(OutputIterator iter) 
{
    std::copy( newVertices_.begin( ), newVertices_.end( ), iter );

    if ( not this->isLeaf( ) )
        for ( unsigned i = 0; i < numChildren; i++ )
            children_[ i ]->addNewVertices( iter );

    return;
}


/// Refine the edge vertex
template < corlib::shape SHAPE, typename V>
template< typename SUBDIV >
void EdgeTree<SHAPE, V>::newCoordEdgeVertices_(Vertex*& v,  SUBDIV * subdiv )
{
    v  = new Vertex(depth_+1); 
    typename Vertex::VecP newPos;
    newPos.setZero();

    subdiv->refineEdge(this, depth_, newPos);

    v ->addNewLevel(newPos, depth_+1); 	
    return;
}


/// Refine the existing vertices
template < corlib::shape SHAPE, typename V>
template< typename SUBDIV >
void EdgeTree<SHAPE, V>::newCoordExistingVertices_( SUBDIV * subdiv )
{
    for (unsigned i=0; i<numVertices; ++i) {

	const int maxLev = vertices_[i]->maxActiveLevel();
        if ((depth_+1)==maxLev)  continue;

        typename Vertex::VecP newPos;
        newPos.setZero();
	subdiv->refineVertex(this, i, depth_, newPos);     
	vertices_[i]->addNewLevel(newPos, depth_+1);	    
    }
 
    return;
}

/// Subdivide edge tree leaf
template < corlib::shape SHAPE, typename V>
template< typename SUBDIV >
void EdgeTree<SHAPE, V>::subdivideLeaf_( SUBDIV * subdiv )
{
    // compute new edge midpoint coordinates
    Vertex* v;
    this->newCoordEdgeVertices_(v, subdiv);
    newVertices_.push_back(v);

    // recompute existing coordinates
    this->newCoordExistingVertices_(subdiv);

    //build the children and their relationships
    subdiv::curve::BuildChildren<ETree>()(this, v);    
    
    return;
}

