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
	//   cycleHolonomy: #basisCycles angle defect for each basis cycle e.g. from directional::cycle_holonomy.
	//   N: the degree of the field
	// Outputs:
	//   singularities: #basisCycles the index around each cycle.
	IGL_INLINE void singularities(const Eigen::SparseMatrix<double, Eigen::RowMajor>& basisCycleMat,
		const Eigen::VectorXd& adjustAngles,
		const Eigen::VectorXd& cycleHolonomy,
		int N,
		Eigen::VectorXd& singularities)
	{

		singularities = ((basisCycleMat * adjustAngles + cycleHolonomy).array() / (2.*igl::PI / N));
	}

	// Computes a matrix containing the singularity values of all cycles
	// Inputs:
	//   basisCycleMat: #basisCycles by #E the oriented basis cycles around which the singularities are measured
	//   adjustAngles: #E angles that encode deviation from parallel transport.
	//   parallelTransportAngles: #E angles used ffor parallel transport between each edge.
	//   N: the degree of the field
	// Outputs:
	//   singularities: #basisCycles the index around each cycle.
	IGL_INLINE void singularities(const Eigen::MatrixXd& V,
		const Eigen::MatrixXi& F,
		const Eigen::SparseMatrix<double, Eigen::RowMajor>& basisCycleMat,
		const Eigen::VectorXd& adjustAngles,
		int N,
		Eigen::VectorXd& singularities)
	{
		Eigen::VectorXd cycleHolonomy;
		cycle_holonomy(V, F, basisCycleMat, cycleHolonomy);
		directional::singularities(basisCycleMat, adjustAngles, cycleHolonomy, N, singularities);
	}

}




#endif


