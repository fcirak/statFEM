// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//
//                 Computational Structural Mechanics Lab
//                         University of Cambridge
//
// This software is copyrighted and all rights are retained by the University of
// Cambridge. All use, disclosure, and/or reproduction of any part not expressly
// authorised by the University of Cambridge is prohibited. (C) 2010
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

//! @file eigenX.hpp

#ifndef corlib_eigenX_h
#define corlib_eigenX_h

//------------------------------------------------------------------------------

#include <Eigen/Core>

namespace corlib {
	namespace eigenX {

		template<unsigned dim>
		using VectorSd = Eigen::Matrix<double, dim, 1>;

		template<unsigned dim>
		using VectorSi = Eigen::Matrix<int, dim, 1>;

		template<unsigned dim1, unsigned dim2 = dim1>
		using MatrixSd = Eigen::Matrix<double, dim1, dim2>;

	}
}


#endif

