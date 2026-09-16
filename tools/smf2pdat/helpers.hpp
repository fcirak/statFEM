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

//! @author Kosala Bandara
//! @date   12/2012

/** \brief Set of I/O functions common to smf2vtu and psmf2vtu convertors.
 */
//------------------------------------------------------------------------------


#ifndef tools_helpers_kkmb2_smf2tg_helpers_h
#define tools_helpers_kkmb2_smf2tg_helpers_h

#include <set>
#include <vector>

#include <corlib/SmfHead.hpp>
#include <subdiv/surf/ShapeDim.hpp>
#include <boost/numeric/ublas/vector.hpp>

//==============================================================================
// declarations
namespace tools{
    namespace helpers{
        namespace smf2tg{

            namespace ublas = boost::numeric::ublas;

            class CompareNodes;
            class CompareEdges;

            void readSmfFiles( std::istream & full,
                               std::istream & region,
                               std::vector< ublas::bounded_vector<double,3> > & fullCoord,
                               std::vector< ublas::bounded_vector<double,3> > & regionCoord,
                               std::set<std::pair<unsigned, unsigned>, CompareEdges> & fullEdgeSet);

            void writeOutput( const std::string & outFilename, 
                              const unsigned numNodes,
                              std::set<unsigned> & taggedNodes,
                              std::set<std::pair<unsigned, unsigned>, CompareEdges> & taggedEdges,
                              const bool writeTag, 
                              const bool writeETag, 
                              const bool writeFix);
        }
    }
}

//==============================================================================
// definitions

// class to compare nodes
class tools::helpers::smf2tg::CompareNodes
{
public:
    typedef ublas::bounded_vector<double,3> Vec;
    
    CompareNodes(std::vector< ublas::bounded_vector<double,3> > & compVec, const double tol)
        : compVec_(compVec), tol_(tol){}
    
    unsigned findNode(const Vec & x) const 
    { 
        unsigned match = compVec_.size();
#if 1
        for (unsigned i=0; i < compVec_.size(); i++){
            Vec v = compVec_[i];
            if ( (corlib::fuzzyEqual(v[0], x[0], tol_) ) and 
                 (corlib::fuzzyEqual(v[1], x[1], tol_) ) and 
                 (corlib::fuzzyEqual(v[2], x[2], tol_) ) ){
                match = i;
                break;
            }
        }
#else
        // Fix grid points one plane (x-p)n =0
        boost::numeric::ublas::bounded_vector<double,3> p, n;
        p[0] = -.27;     p[1] = -.12;     p[2] = -.02;
        n[0] = .06;      n[1] = -0.54;    n[2] = 0.83;
        
        const double distFromPlane = boost::numeric::ublas::inner_prod((x - p) , n);
        const double tol = 0.05;
        if ( ( corlib::fuzzyEqual(distFromPlane, 0., tol) ) and  
             (x[0] > -.59 ) and
             (x[0] <  .08 ) and
             (x[1] > -.3) and
             (x[1] < .22)){
            match = 0;
        }
            
        
#endif
        return match;
    }
    
private:    
    std::vector< ublas::bounded_vector<double,3> > compVec_;
    double                                         tol_;
};

// class to compare edges
class tools::helpers::smf2tg::CompareEdges
{
public:
    typedef std::pair<unsigned, unsigned> Edge;
    
    bool operator()(const Edge & lhs, const Edge & rhs)const
    { 
        const unsigned m1 = std::min( lhs.first, lhs.second );
        const unsigned m2 = std::min( rhs.first, rhs.second );
        return ( ( m1 < m2 ) or
                 ( ( m1 == m2 ) and
                   ( std::max( lhs.first, lhs.second ) <
                     std::max( rhs.first, rhs.second ) ) ) );        
    }
};

// read input files
void tools::helpers::smf2tg::readSmfFiles( 
    std::istream & full,
    std::istream & region,
    std::vector< ublas::bounded_vector<double,3> > & fullCoord,
    std::vector< ublas::bounded_vector<double,3> > & regionCoord,
    std::set<std::pair<unsigned, unsigned>, CompareEdges> & fullEdgeSet)

{
    // try to read header
    corlib::SmfHead smfHead1, smfHead2;
    smfHead1.read( full );
    smfHead2.read( region );
    // FTL_VERIFY(smfHead1.giveElementNumPoints() == smfHead2.giveElementNumPoints() );
    
    unsigned numNodes, numElements, numNodesPerElem;
    numNodesPerElem = smfHead1.giveElementNumPoints();
    corlib::shape myShape = smfHead1.giveElementShape();

    // read full mesh
    full >> numNodes >> numElements;
    full.ignore( std::numeric_limits<std::streamsize>::max(), '\n' );

    // read coordinates
    for ( unsigned n = 0; n < numNodes; n ++ ) {
        boost::numeric::ublas::bounded_vector<double,3> node;
        for ( unsigned d = 0; d < 3; d ++ ) {
            full >> node[d];
        }
        fullCoord.push_back( node );
        full.ignore( std::numeric_limits<std::streamsize>::max(), '\n' );
    }

    // create edge set
    for ( unsigned e = 0; e < numElements; e ++ ) {
        // read connectivity 
        boost::numeric::ublas::vector<unsigned> connectivity(numNodesPerElem);
        for (unsigned i = 0; i < numNodesPerElem; i++)
            full >> connectivity[i];

        // add edges
        for (unsigned i = 0; i < numNodesPerElem; i++) {
            unsigned v1 = connectivity[i];
            unsigned v2;
            if (myShape == corlib::TRIANGLE) 
                v2=connectivity[ subdiv::surf::ShapeProp< corlib::TRIANGLE >::next[i]];
            else if (myShape ==corlib::QUADRILATERAL)
                v2=connectivity[ subdiv::surf::ShapeProp< corlib::QUADRILATERAL >::next[i]];
            else {
                std::cout << "unknown shape" <<std::endl;
                assert(false);
            }
            fullEdgeSet.insert(std::make_pair(v1, v2));
        }
    }

    // read region mesh
    region >> numNodes >> numElements;
    region.ignore( std::numeric_limits<std::streamsize>::max(), '\n' );

    // read coordinates
    for ( unsigned n = 0; n < numNodes; n ++ ) {
        boost::numeric::ublas::bounded_vector<double,3> node;
        for ( unsigned d = 0; d < 3; d ++ ) {
            region >> node[d];
        }
        regionCoord.push_back( node );
        region.ignore( std::numeric_limits<std::streamsize>::max(), '\n' );
    }
    return;
}

void tools::helpers::smf2tg::writeOutput( const std::string & outFilename, 
                                          const unsigned numNodes,
                                          std::set<unsigned> & taggedNodes, 
                                          std::set<std::pair<unsigned, unsigned>, 
                                                   CompareEdges> & taggedEdges,
                                          const bool writeTag, 
                                          const bool writeETag, 
                                          const bool writeFix)
{
    FTL_VERIFY( taggedNodes.size() );

    // write out pdat file
    std::ofstream pdat(outFilename.c_str());
    pdat << numNodes << " 1"<<std::endl;
    pdat << "tagged 1"<<std::endl;
    
    std::vector<unsigned> fix;

    for(unsigned i =0; i < numNodes; i ++){
        if (taggedNodes.find(i) == taggedNodes.end())
            pdat <<"0" <<std::endl;
        else {
            pdat <<"1" <<std::endl;
            fix.push_back(i);
        }
    }
    pdat.close();

    // write out fix file
    if (writeFix){
        std::string fixFilename = outFilename.substr( 
            0, outFilename.find( ".pdat" ) ) + ".fix";
        std::ofstream fixOut(fixFilename.c_str());
        fixOut << fix.size()<<std::endl;
        
        for(unsigned i = 0; i < fix.size(); i++) fixOut << fix[i] <<" "; 
        fixOut <<std::endl;
        fixOut.close();
    }
    

    // write out tag file
    if (writeTag){
        std::string tagFilename = outFilename.substr( 
            0, outFilename.find( ".pdat" ) ) + ".tg";
        std::ofstream tagOut(tagFilename.c_str());
        tagOut << "TG" <<std::endl;
        tagOut << fix.size()<< " 0" <<std::endl;
        
        for(unsigned i = 0; i < fix.size(); i++){
#if 1
            tagOut << fix[i] <<" 1" <<std::endl; 
#else
            tagOut << fix[i] <<" -0.001 0.001 -0.001 0.001 0.0 12.0" <<std::endl; 
#endif
        }
        tagOut.close();
    }

    // write out edge tag file
    if (writeETag){
        std::string tagFilename = outFilename.substr( 
            0, outFilename.find( ".pdat" ) ) + ".tg";
        std::ofstream tagOut(tagFilename.c_str());
        tagOut << "TG" <<std::endl;
        tagOut << fix.size()<< " " << taggedEdges.size() <<std::endl;
        
        // write vertex tags
        for(unsigned i = 0; i < fix.size(); i++) 
            tagOut << fix[i] <<" 1" <<std::endl; 

        // write edge tags
        typedef std::set<std::pair<unsigned, unsigned>, CompareEdges>::iterator EdgeSetIter;
        EdgeSetIter sBegin =taggedEdges.begin();
        EdgeSetIter sEnd =taggedEdges.end();
        for (; sBegin != sEnd; sBegin++)
            tagOut << (*sBegin).first <<" "<< (*sBegin).second << " 1" <<std::endl; 
        tagOut.close();
    }
    
}

#endif
