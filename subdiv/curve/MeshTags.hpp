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

#ifndef subdiv_curve_mesh_tags_h
#define subdiv_curve_mesh_tags_h

#include <iostream>
#include <string>

namespace subdiv{
    namespace curve{
        //--------------------------------------------------------------------------
	//! \brief Enumerator for description of vertex tags
	enum vertexTag {
	    NOTAG_VERTEX, 
	    NOSUB_VERTEX,
            NEXT_NOSUB_VERTEX
        };

        //--------------------------------------------------------------------------
	//! \brief Enumerator for description of refinement status
        enum refinementStatus {
            INACTIVE,           //!< inactive
            ACTIVE,             //!< active
            DEPENDENT           //!< active but linearly dependent
        };
    } // end namespace curve
} // end namespace subdiv
#endif
