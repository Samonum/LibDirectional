// Copyright (C) 2016 Amir Vaxman <avaxman@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
#ifndef TRIVIAL_CONNECTIONS_H
#define TRIVIAL_CONNECTIONS_H
#include <igl/igl_inline.h>
#include <igl/gaussian_curvature.h>
#include <igl/local_basis.h>
#include <igl/parallel_transport_angles.h>
#include <igl/edge_topology.h>
#include <igl/gaussian_curvature.h>
#include <igl/boundary_loop.h>
#include <directional/cycle_holonomy.h>
#include <Eigen/Core>
#include <vector>
#include <cmath>


namespace directional
{    
    // Computes the adjustment angles to form a trivial connection according to given cone curvatures (or singularity indices) around basis cycles. 
	// In case the sum of curvature is not consistent with the topology, the system is solved in least squares and unexpected singularities may appear elsewhere. 
	// The output is the modification to the parallel transport.
    // Inputs:
    //  V: #V X 3 vertex coordinates
    //  F: #F by 3 face vertex indices
	//  basisCycles: #basisCycles X #E the oriented (according to EF) basis cycles around which the curvatures are measured
	//               the basis cycles must be arranged so that the first |V| are the vertex cycles (for instance, the result of igl::basis_cycles())
    //  indices: #basisCycles the index around each cycle. They should add up to N*Euler_characteristic of the mesh.
	//  cycleHolonomy: the angle defect for each basis cycle.
	//  solver: The Simplicial LDLT solver used to calculate the trivial connections. If initialized the solve step will be skipped when calculating the field.
	//			The state of  the solver solely depends on the basisCycles, therefore it only needs to be reset if the basisCycles matrix changed.
	//			If the solver is not yet set the solver will be called to prepare the basisCycles.
    //  N: the degree of the field. The curvature of a cycle is measured by (singIndex/N)*(2*pi) (can be negative)
    // Outputs:
    //  adjustAngles: the difference between the parallel transport and the modified one.
	//  error: gives the total error of the field. If this is not approximately 0 your singularities probably don't add up properly.
	IGL_INLINE void trivial_connection(const Eigen::MatrixXd& V,
		const Eigen::MatrixXi& F,
		const Eigen::SparseMatrix<double, Eigen::RowMajor>& basisCycles,
		const Eigen::VectorXd& indices,
		const Eigen::VectorXd& cycleHolonomy,
		Eigen::SimplicialLDLT<Eigen::SparseMatrix<double> >& solver,
		const int N,
		Eigen::VectorXd& adjustAngles,
		double &error)
	{
		using namespace Eigen;
		using namespace std;

		VectorXi rows = ArrayXi::Zero(basisCycles.rows());
		VectorXi columns = ArrayXi::Zero(basisCycles.cols());

		//Filter out empty rows and columns
		for (int k = 0; k < basisCycles.outerSize(); ++k)
			for (SparseMatrix<double, Eigen::RowMajor>::InnerIterator it(basisCycles, k); it; ++it)
			{
				if (!it.value())
					continue;
				rows[it.row()] = 1;
				columns[it.col()] = 1;
			}

		for (int i = 1; i < rows.size(); i++)
			rows[i] += rows[i - 1];
		for (int i = 1; i < columns.size(); i++)
			columns[i] += columns[i - 1];

		SparseMatrix<double, Eigen::RowMajor> reducedCycles(rows[rows.size() - 1], columns[columns.size() - 1]);
		VectorXd reducedIndices(rows[rows.size() - 1]);
		reducedCycles.reserve(VectorXi::Constant(columns[columns.size() - 1], 2));

		for (int k = 0; k < basisCycles.outerSize(); ++k)
			for (SparseMatrix<double, Eigen::RowMajor>::InnerIterator it(basisCycles, k); it; ++it)
				if (it.value())
					reducedCycles.insert(rows(it.row()) - 1, columns(it.col()) - 1) = it.value();

		for (int i = 0; i < indices.size(); i++)
			if((i == 0 && rows(i) == 1) || (i > 0 && rows(i) > rows(i-1)))
				reducedIndices(rows(i) - 1) = indices(i);

		VectorXd  reducedCycleHolonomy(rows(rows.size()-1));

		for (int i = 0; i < cycleHolonomy.size(); i++)
			if ((i == 0 && rows(i) == 1) || (i > 0 && rows(i) > rows(i - 1)))
				reducedCycleHolonomy(rows(i) - 1) = cycleHolonomy(i);


		VectorXd cycleNewCurvature = reducedIndices*(2.0*M_PI / (double)N);

		//Run Solver if uninitialised
		if (!solver.rows())
		{
			SparseMatrix<double> bbt = reducedCycles*reducedCycles.transpose();
			solver.compute(bbt);
		}

		VectorXd reducedAngles = reducedCycles.transpose()*solver.solve((-reducedCycleHolonomy + cycleNewCurvature));

		adjustAngles.resize(columns.size());

		if (columns[0] == 1)
			adjustAngles[0] = reducedAngles[0];

		for (int i = 1; i < columns.size(); i++)
			if (columns[i] != columns[i - 1])
				adjustAngles[i] = reducedAngles[columns[i - 1]];
		error = (reducedCycles*reducedAngles - (-reducedCycleHolonomy + cycleNewCurvature)).lpNorm<Infinity>();
	}

	// Computes the adjustment angles to form a trivial connection according to given cone curvatures (or singularity indices) around basis cycles. 
	// In case the sum of curvature is not consistent with the topology, the system is solved in least squares and unexpected singularities may appear elsewhere. 
	// The output is the modification to the parallel transport.
	// Inputs:
	//  V: #V X 3 vertex coordinates
	//  F: #F by 3 face vertex indices
	//  basisCycles: #basisCycles X #E the oriented (according to EF) basis cycles around which the curvatures are measured
	//               the basis cycles must be arranged so that the first |V| are the vertex cycles (for instance, the result of igl::basis_cycles())
	//  indices: #basisCycles the index around each cycle. They should add up to N*Euler_characteristic of the mesh.
	//  cycleHolonomy: the angle defect for each basis cycle.
	//  N: the degree of the field. The curvature of a cycle is measured by (singIndex/N)*(2*pi) (can be negative)
	// Outputs:
	//  adjustAngles: the difference between the parallel transport and the modified one.
	//  error: gives the total error of the field. If this is not approximately 0 your singularities probably don't add up properly.
	//TODO: work with boundaries
	IGL_INLINE void trivial_connection(const Eigen::MatrixXd& V,
		const Eigen::MatrixXi& F,
		const Eigen::SparseMatrix<double, Eigen::RowMajor>& basisCycles,
		const Eigen::VectorXd& indices,
		const Eigen::VectorXd& cycleHolonomy,
		const int N,
		Eigen::VectorXd& adjustAngles,
		double &error)
	{
		Eigen::SimplicialLDLT<Eigen::SparseMatrix<double> > solver;
		trivial_connection(V, F, basisCycles, indices, cycleHolonomy, solver, N, adjustAngles, error);
	}

	// Computes the adjustment angles to form a trivial connection according to given cone curvatures (or singularity indices) around basis cycles. 
	// In case the sum of curvature is not consistent with the topology, the system is solved in least squares and unexpected singularities may appear elsewhere. 
	// The output is the modification to the parallel transport.
	// Inputs:
	//  V: #V X 3 vertex coordinates
	//  F: #F by 3 face vertex indices
	//  basisCycles: #basisCycles X #E the oriented (according to EF) basis cycles around which the curvatures are measured
	//               the basis cycles must be arranged so that the first |V| are the vertex cycles (for instance, the result of igl::basis_cycles())
	//  indices: #basisCycles the index around each cycle. They should add up to N*Euler_characteristic of the mesh.
	//  solver: The Simplicial LDLT solver used to calculate the trivial connections. If initialized the solve step will be skipped when calculating the field.
	//			The state of  the solver solely depends on the basisCycles, therefore it only needs to be reset if the basisCycles matrix changed.
	//			If the solver is not yet set the solver will be called to prepare the basisCycles.
	//  N: the degree of the field. The curvature of a cycle is measured by (singIndex/N)*(2*pi) (can be negative)
	// Outputs:
	//  adjustAngles: the difference between the parallel transport and the modified one.
	//  error: gives the total error of the field. If this is not approximately 0 your singularities probably don't add up properly.
	IGL_INLINE void trivial_connection(const Eigen::MatrixXd& V,
		const Eigen::MatrixXi& F,
		const Eigen::SparseMatrix<double, Eigen::RowMajor>& basisCycles,
		const Eigen::VectorXd& indices,
		Eigen::SimplicialLDLT<Eigen::SparseMatrix<double> >& solver,
		const int N,
		Eigen::VectorXd& adjustAngles,
		double &error)
	{
		Eigen::MatrixXi EV, x, EF;
		igl::edge_topology(V, F, EV, x, EF);
		Eigen::MatrixXd B1, B2, B3;
		igl::local_basis(V, F, B1, B2, B3);
		Eigen::VectorXd cycleHolonomy;
		cycle_holonomy(V, F, EV, EF, B1, B2, basisCycles, cycleHolonomy);
		trivial_connection(V, F, basisCycles, indices, cycleHolonomy, solver, N, adjustAngles, error);
	}

	// Computes the adjustment angles to form a trivial connection according to given cone curvatures (or singularity indices) around basis cycles. 
	// In case the sum of curvature is not consistent with the topology, the system is solved in least squares and unexpected singularities may appear elsewhere. 
	// The output is the modification to the parallel transport.
	// Inputs:
	//  V: #V X 3 vertex coordinates
	//  F: #F by 3 face vertex indices
	//  basisCycles: #basisCycles X #E the oriented (according to EF) basis cycles around which the curvatures are measured
	//               the basis cycles must be arranged so that the first |V| are the vertex cycles (for instance, the result of igl::basis_cycles())
	//  indices: #basisCycles the index around each cycle. They should add up to N*Euler_characteristic of the mesh.
	//  N: the degree of the field. The curvature of a cycle is measured by (singIndex/N)*(2*pi) (can be negative)
	// Outputs:
	//  adjustAngles: the difference between the parallel transport and the modified one.
	//  error: gives the total error of the field. If this is not approximately 0 your singularities probably don't add up properly.
	IGL_INLINE void trivial_connection(const Eigen::MatrixXd& V,
		const Eigen::MatrixXi& F,
		const Eigen::SparseMatrix<double, Eigen::RowMajor>& basisCycles,
		const Eigen::VectorXd& indices,
		const int N,
		Eigen::VectorXd& adjustAngles,
		double &error)
	{
		Eigen::MatrixXi EV, x, EF;
		igl::edge_topology(V, F, EV, x, EF);
		Eigen::MatrixXd B1, B2, B3;
		igl::local_basis(V, F, B1, B2, B3);
		Eigen::VectorXd cycleHolonomy;
		cycle_holonomy(V, F, EV, EF, B1, B2, basisCycles, cycleHolonomy); 
		trivial_connection(V, F, basisCycles, indices, cycleHolonomy, N, adjustAngles, error);
	}
}

    


#endif


