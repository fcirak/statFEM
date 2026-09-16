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

//! @author Kosala Bandara
//! @date   2011
//------------------------------------------------------------------------------
/// Local index of next vertex
template< typename ET >
const std::array< int, ET::numVertices >
subdiv::curve::BuildChildren< ET >::next_ =
    subdiv::curve::ShapeProp< corlib::LINE >::next;

/// Local index of previous vertex
template< typename ET >
const std::array< int, ET::numVertices >
subdiv::curve::BuildChildren< ET >::previous_ =
    subdiv::curve::ShapeProp< corlib::LINE >::previous;


//------------------------------------------------------------------------------
template< typename ET >
void subdiv::curve::BuildChildren< ET >::operator()(ETree * et, Vertex * ev)
{

    // make children     
    typename ET::VecVPtrNV vtxs;
    
    // build children 
    for (unsigned i = 0; i < ET::numChildren; ++i) {
        vtxs[i] = et->vertex(i);    
        vtxs[next_[i]] = ev;    
        
        et->child(i) = new ET(vtxs, et->depth()+1, et); 
        assert(et->child(i) != NULL);
            
    }
    
    //update vertex tags
    const int level = et->depth();
    for (unsigned i = 0; i < ET::numVertices; ++i) {
        // set the next_nosub tag correctly
        if (et->vertex(i)->tag(level)==subdiv::curve::NOSUB_VERTEX and
            et->vertex(next_[i])->tag(level)==subdiv::curve::NEXT_NOSUB_VERTEX){
            ev->setTag(subdiv::curve::NEXT_NOSUB_VERTEX, level+1);                    
            break;
        } else 
            ev->setTag(subdiv::curve::NOTAG_VERTEX, level+1);                
    }

    for (unsigned i = 0; i < ET::numVertices; ++i) {
        //set other tags
        if (et->vertex(i)->tag(level)==subdiv::curve::NOSUB_VERTEX)
            et->vertex(i)->setTag(subdiv::curve::NOSUB_VERTEX, level+1);
        else
            et->vertex(i)->setTag(subdiv::curve::NOTAG_VERTEX, level+1);
    }
    
    // neighbor relations within children
    for (unsigned i = 0; i < ET::numNeighbors; ++i)
        et->child(i)->neighbor(next_[i]) = et->child(next_[i]);
    
    // amongst other elements    
    ET *en;
    for (unsigned n = 0; n < ET::numNeighbors; n++){
        en = et->neighbor(n);
        if (en and !en->isLeaf()){
            int j = localVtxNumber(en, et->vertex(n));
            en->child(j)->neighbor(j)= et->child(n);
            et->child(n)->neighbor(n) = en->child(j);
        }
    } 
    
    return;
}
