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
#include <igl/parallel_transport_angles.h>

#include <Eigen/Core>
#include <vector>
#include <cmath>


namespace directional
{
	// Computes a matrix containing the singularity values of all cycles
	// Inputs:
	//   basisCycleMat: #basisCycles by #E the oriented basis cycles around which the singularities are measured
	//   matching: The matching between the representative edge in one face and its closest match in the other.
	// Outputs:
	//   singularities: #basisCycles the index around each cycle.
	IGL_INLINE void singularities(const Eigen::SparseMatrix<double, Eigen::RowMajor>& basisCycleMat,
		const Eigen::VectorXi& matching,
		Eigen::VectorXd& singularities)
	{
		singularities = basisCycleMat * matching.cast<double>();
	}	

	// Computes a matrix containing the singularity values of all cycles
	// Inputs:
	//   basisCycleMat: #basisCycles by #E the oriented basis cycles around which the singularities are measured
	//   adjustAngles: #E angles that encode deviation from parallel transport.
	//   parallelTransportAngles: #E angles used ffor parallel transport between each edge.
	//   N: the degree of the field
	// Outputs:
	//   singularities: #basisCycles the index around each cycle.
	IGL_INLINE void singularities(const Eigen::SparseMatrix<double, Eigen::RowMajor>& basisCycleMat,
		const Eigen::VectorXd& adjustAngles,
		const Eigen::VectorXd& parallelTransportAngles,
		int N,
		Eigen::VectorXd& singularities)
	{
		Eigen::VectorXd cycleHolonomy = basisCycleMat*parallelTransportAngles;
		for (int i = 0; i < cycleHolonomy.size(); i++) {
			while (cycleHolonomy(i) >= M_PI) cycleHolonomy(i) -= 2.0*M_PI;
			while (cycleHolonomy(i) < -M_PI) cycleHolonomy(i) += 2.0*M_PI;
		}

		singularities = ((basisCycleMat * angles + cycleHolonomy).array() / (2.*igl::PI / N));
	}
}




#endif


