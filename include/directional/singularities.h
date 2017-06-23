// Copyright (C) 2016 Amir Vaxman <avaxman@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
#ifndef SINGULARITIES_H
#define SINGULARITIES_H
#include <igl/igl_inline.h>
#include <igl/edge_topology.h>
#include <igl/parallel_transport_angles.h>
#include <igl/per_face_normals.h>
#include <igl/parallel_transport_angles.h>
#include <directional/cycle_holonomy.h>

#include <Eigen/Core>
#include <vector>
#include <cmath>


namespace directional
{
	// Computes a matrix containing the singularity values of all cycles
	// Inputs:
	//   basisCycleMat: #basisCycles by #E the oriented basis cycles around which the singularities are measured
	//   effort: The effort required to match between two faces, equal to the adjustment angles for n-rosies.
	//   cycleHolonomy: #basisCycles angle defect for each basis cycle e.g. from directional::cycle_holonomy.
	//   N: the degree of the field
	// Outputs:
	//   singularities: #basisCycles the index around each cycle.
	IGL_INLINE void singularities(const Eigen::SparseMatrix<double, Eigen::RowMajor>& basisCycleMat,
		const Eigen::VectorXd& adjustAngles,
		const Eigen::VectorXd& cycleHolonomy,
		int N,
		Eigen::VectorXi& singularities)
	{

		Eigen::VectorXd unrounded = ((basisCycleMat * adjustAngles + cycleHolonomy).array() / (2.*igl::PI / N));
		singularities.resize(unrounded.size());
		for (int i = 0; i < unrounded.size(); i++)
			singularities(i) = round(unrounded(i));
	}

	// Computes a matrix containing the singularity values of all cycles
	// Inputs:
	//   V: #V by 3 vertex coordinates
	//   F: #F by 3 face vertex indices
	//   basisCycleMat: #basisCycles by #E the oriented basis cycles around which the singularities are measured
	//   effort: The effort required to match between two faces, equal to the adjustment angles for n-rosies.
	//   N: the degree of the field
	// Outputs:
	//   singularities: #basisCycles the index around each cycle.
	IGL_INLINE void singularities(const Eigen::MatrixXd& V,
		const Eigen::MatrixXi& F,
		const Eigen::SparseMatrix<double, Eigen::RowMajor>& basisCycleMat,
		const Eigen::VectorXd& effort,
		int N,
		Eigen::VectorXi& singularities)
	{
		Eigen::VectorXd cycleHolonomy;
		cycle_holonomy(V, F, basisCycleMat, cycleHolonomy);
		directional::singularities(basisCycleMat, effort, cycleHolonomy, N, singularities);
	}

}




#endif


