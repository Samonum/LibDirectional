// Copyright (C) 2016 Amir Vaxman <avaxman@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
#ifndef SINGULARITIES_H
#define SINGULARITIES_H
#include <igl/igl_inline.h>
#include <igl/gaussian_curvature.h>
#include <igl/local_basis.h>
#include <igl/triangle_triangle_adjacency.h>

#include <Eigen/Core>
#include <vector>
#include <cmath>


namespace directional
{
	// Takes a vector-set field, and creates a principal matching: that is, an order-preserving matching with effort between
	// [-pi, pi).
	// Inputs:
	//  V: #V X 3 vertex coordinates
	//  F: #F by 3 face vertex indices
	//  EV: #E x 2 edges 2 vertices indices
	//  EF: #E X 2 edges 2 faces indices
	//  vectorSetField: the vector set field, assumed to be ordered CCW, and in xyzxyzxyz...xyz (3*N cols) form. The degree is inferred by the size.
	//  INPORTANT: if the vector set field is not CCW ordered, the result in unpredictable
	// Outputs:
	//  matching: #E matching differences from EF(i,0) to EF(i,1). an integer "k" indicates that vector label j goes to (j+k)%N
	//  effort: #E matching efforts. One can compute singularity indices by a (properly aligned) d0 operator as (d0'*Effort+K)/(2*pi*N)
	// note: in some (probably extreme) cases, matching vector to vector and taking the individual closest angles is not the actual effort, since it can break the order. That is, principal effort != every matching is shortest angle. For now, we disregard this, but this is TODO.
	IGL_INLINE void singularities(const Eigen::SparseMatrix<double, Eigen::RowMajor>& basisCycleMat,
		const Eigen::VectorXi& matching,
		Eigen::VectorXd& singularities)
	{
		singularities = basisCycleMat * matching.cast<double>();
	}
}




#endif


