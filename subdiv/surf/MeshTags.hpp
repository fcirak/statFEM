//
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

#ifndef subdiv_surf_meshtags_h
#define subdiv_surf_meshtags_h

#include <iostream>
#include <string>
#include <array>

namespace subdiv{
    namespace surf{

        //======================================================================
        /// Enumerator for description of vertex tags
        enum vertexTag {
            VERTEX_NOTAG,               ///< Default 'untagged' vertex
            VERTEX_CORNER,              ///< (Convex) corner, i.e do not sub-divide (fix-vertex)
            VERTEX_CREASE,              ///< Vertex on crease OR boundary OR non-manifold edge
            VERTEX_DART,                ///< End point of a crease with nicely blended
                                        ///< continuity to surrounding surface 
            VERTEX_SMOOTHCORNER,        ///< Vertex in between a pair-wise smooth connection
                                        ///< of two manifolds. (Please excuse oxymoron.)
            VERTEX_SMOOTHCREASE         ///< "Smooth corner" at end of pair-wise smooth connection
                                        ///< of two manifolds. (Please excuse oxymoron.)
        };

        /// Out-stream operator for vertex tags
        ///
        /// \param[in,out]   os      Output stream
        /// \param[in]       eTag    Vertex tag
        /// \return                  Output stream
        std::ostream & operator<<( std::ostream & os, const enum vertexTag vTag )
        {
            std::array< std::string, 6 > vNames = { { "VERTEX_NOTAG",
                                                      "VERTEX_CORNER",
                                                      "VERTEX_CREASE",
                                                      "VERTEX_DART",
                                                      "VERTEX_SMOOTHCORNER",
                                                      "VERTEX_SMOOTHCREASE" } };
            os << vNames[ vTag ];
            return os; 
        }
        
        //======================================================================
        /// Enumerator for description of EdgeTags
        enum edgeTag {
            EDGE_NOTAG,                 ///< Default 'untagged' edge
            EDGE_CREASE,                ///< Crease
            EDGE_SMOOTHCREASE,          ///< Pair-wise smooth connection of
                                        ///< two manifolds. Maybe deactivated
                                        ///< in subdivision. 
                                        ///< (Please excuse oxymoron.) 
            // NOTE: Keep the following internal tags at the back
            EDGE_ASKNEIGHBOR,           ///< Ask neighbour facet at edge
            EDGE_ASKPARENT              ///< Ask parent facet at edge

        };

        /// Out-stream operator for edge tags
        ///
        /// \param[in,out]   os      Output stream
        /// \param[in]       eTag    Edge tag
        /// \return                  Output stream
        std::ostream & operator<<( std::ostream & os, const enum edgeTag eTag )
        {
            std::array< std::string, 5 > eNames = { { "EDGE_NOTAG",
                                                      "EDGE_CREASE",
                                                      "EDGE_SMOOTHCREASE",
                                                      "EDGE_ASKNEIGHBOR",
                                                      "EDGE_ASKPARENT" } };
            os << eNames[ eTag ];
            return os; 
        }

        //======================================================================
        /// Edge category. This is a topological information which can be
        /// unambiguously deduced from the mesh itself.
        enum edgeCategory {
            ECAT_UNKNOWN,             ///< Undetermined
            ECAT_BOUNDARY,            ///< Boundary edge, i.e. is only
                                      ///< attached to single facet (singular)
            ECAT_MANIFOLD,            ///< Manifold edge, i.e. is shared
                                      ///< by exactly two facets (non-singular)
            ECAT_NONMANIFOLD          ///< Non-manifold edge, i.e. is shared
                                      ///< by more than two facets (singular)
        };

        /// Out-stream operator for edge categories
        ///
        /// \param[in,out]   os      Output stream
        /// \param[in]       eCat    Edge category
        /// \return                  Output stream
        std::ostream & operator<<( std::ostream & os, const enum edgeCategory eCat )
        {
            std::array< std::string, 4 > eNames = { { "ECAT_UNKNOWN",
                                                      "ECAT_BOUNDARY",
                                                      "ECAT_MANIFOLD",
                                                      "ECAT_NONMANIFOLD" } };
            os << eNames[ eCat ];
            return os; 
        }

    } // end namespace twoD
} // end namespace geom

#endif
