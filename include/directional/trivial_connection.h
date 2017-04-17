// Copyright (C) 2016 Amir Vaxman <avaxman@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
#ifndef IGL_TRIVIAL_CONNECTIONS_H
#define IGL_TRIVIAL_CONNECTIONS_H
#include <igl/igl_inline.h>
#include <igl/gaussian_curvature.h>
#include <igl/local_basis.h>
#include <igl/parallel_transport_angles.h>
#include <igl/edge_topology.h>
#include <igl/gaussian_curvature.h>
#include <igl/boundary_loop.h>
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
    //  EF: #E X 2 edges 2 faces indices
	//  EV: #E by 2 indices of the vertices surrounding each edge
    //  B1, B2: #F by 3 matrices representing the local base of each face.
	//  basisCycles: #basisCycles X #E the oriented (according to EF) basis cycles around which the curvatures are measured
	//               the basis cycles must be arranged so that the first |V| are the vertex cycles (for instance, the result of igl::basis_cycles())
    //  indices: #basisCycles the index around each cycle. They should add up to N*Euler_characteristic of the mesh.
    //  N: the degree of the field. The curvature of a cycle is measured by (singIndex/N)*(2*pi) (can be negative)
    // Outputs:
    //  adjustAngles: the difference between the parallel transport and the modified one.
	//  error: gives the total error of the field. If this is not approximately 0 your singularities probably don't add up properly.
    //TODO: work with boundaries
	IGL_INLINE void trivial_connection(const Eigen::MatrixXd& V,
		const Eigen::MatrixXi& F,
		const Eigen::MatrixXi& EV,
		const Eigen::MatrixXi& EF,
		const Eigen::MatrixXd& B1,
		const Eigen::MatrixXd& B2,
		const Eigen::SparseMatrix<double, Eigen::RowMajor>& basisCycles,
		const Eigen::VectorXd& indices,
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

		// Start of to be replaced section

		VectorXd edgeParallelAngleChange(columns(columns.size() - 1));  //the difference in the angle representation of edge i from EF(i,0) to EF(i,1)
		//MatrixXd edgeVectors(columns(columns.size() - 1), 3);

		for (int i = 0, j = 0; i < EF.rows(); i++) {
			//skip border edges
			if (EF(i,0) == -1 || EF(i,1) ==  -1)
				continue;

			RowVectorXd edgeVectors = (V.row(EV(i, 1)) - V.row(EV(i, 0))).normalized();
			double x1 = edgeVectors.dot(B1.row(EF(i, 0)));
			double y1 = edgeVectors.dot(B2.row(EF(i, 0)));
			double x2 = edgeVectors.dot(B1.row(EF(i, 1)));
			double y2 = edgeVectors.dot(B2.row(EF(i, 1)));
			edgeParallelAngleChange(j) = atan2(y2, x2) - atan2(y1, x1);
			j++;
		}
		
		
		VectorXd cycleHolonomy = reducedCycles*edgeParallelAngleChange;
		for (int i = 0; i < cycleHolonomy.size(); i++) {
			while (cycleHolonomy(i) >= M_PI) cycleHolonomy(i) -= 2.0*M_PI;
			while (cycleHolonomy(i) < -M_PI) cycleHolonomy(i) += 2.0*M_PI;
		}

		// End of to be replaced section

		//VectorXd cycleHolonomy(reducedCycles.rows());
		VectorXd gc;
		igl::gaussian_curvature(V, F, gc);

		//Copy over for contractible cycles
		for (int i = 0; i < gc.size(); i++)
			if ((i == 0 && rows(i) == 1) || (i > 0 && rows(i) > rows(i - 1)))
				cycleHolonomy(rows(i) - 1) = gc(i);
		Eigen::VectorXi border(V.rows());

		//Found border vertices
		border(0) = 1-rows(0);
		for (int i = 1; i < border.rows(); i++)
			border(i) = rows(i) == rows(i - 1);

		//Sum holonomy for boundaries
		int i = 0, added = 0, total = border.sum();
		while (added < total)
		{
			double temp = cycleHolonomy(V.rows() - total + i);
			cycleHolonomy(V.rows() - total + i) = 0;
			//Loop over boundary cycles
			for (SparseMatrix<double, Eigen::RowMajor>::InnerIterator it(basisCycles, V.rows() + i); it; ++it)
				if (it.value())
					if (border(EV(it.col(), 0)) || border(EV(it.col(), 1)))
					{
						int b = border(EV(it.col(), 1));
						// Substract Pi, as igl::gaussian_curvature calculates gc as 2 * PI - sum(angles)
						cycleHolonomy(V.rows() - total + i) += gc(EV(it.col(), b)) - igl::PI;
						cout << "Curvature at border vertex: " << gc(EV(it.col(), b)) - igl::PI << endl;
						//Ensure vertices won't be counted twice
						border(EV(it.col(), b)) = 0;
						added++;
					}
			cout << "Old holonomy: " << temp << " " << "    New holonomy: " << cycleHolonomy(V.rows() - total + i) << endl;
			i++;
		}



		VectorXd cycleNewCurvature = reducedIndices*(2.0*M_PI / (double)N);


		SimplicialLDLT<SparseMatrix<double> > solver;

		SparseMatrix<double> bbt = reducedCycles*reducedCycles.transpose();
		solver.compute(bbt);
		VectorXd reducedAngles = reducedCycles.transpose()*solver.solve((-cycleHolonomy + cycleNewCurvature));

		adjustAngles.resize(columns.size());

		if (columns[0] == 1)
			adjustAngles[0] = reducedAngles[0];

		for (int i = 1; i < columns.size(); i++)
			if (columns[i] != columns[i - 1])
				adjustAngles[i] = reducedAngles[columns[i - 1]];
		std::cout << "Total holonomy: " << cycleHolonomy.sum() << " Prescribed holonomy:" << cycleNewCurvature.sum();
		error = (reducedCycles*reducedAngles - (-cycleHolonomy + cycleNewCurvature)).lpNorm<Infinity>();
	}

	// Computes the adjustment angles to form a trivial connection according to given cone curvatures (or singularity indices) around basis cycles. 
	// In case the sum of curvature is not consistent with the topology, the system is solved in least squares and unexpected singularities may appear elsewhere. 
	// The output is the modification to the parallel transport.
	// Inputs:
	//  V: #V X 3 vertex coordinates
	//  F: #F by 3 face vertex indices
	//  EF: #E X 2 edges 2 faces indices
	//  EV: #E by 2 indices of the vertices surrounding each edge
	//  basisCycles: #basisCycles X #E the oriented (according to EF) basis cycles around which the curvatures are measured
	//               the basis cycles must be arranged so that the first |V| are the vertex cycles (for instance, the result of igl::basis_cycles())
	//  indices: #basisCycles the index around each cycle. They should add up to N*Euler_characteristic of the mesh.
	//  N: the degree of the field. The curvature of a cycle is measured by (singIndex/N)*(2*pi) (can be negative)
	// Outputs:
	//  adjustAngles: the difference between the parallel transport and the modified one.
	//  error: gives the total error of the field. If this is not approximately 0 your singularities probably don't add up properly.
	//TODO: work with boundaries
	IGL_INLINE void trivial_connection(const Eigen::MatrixXd& V,
		const Eigen::MatrixXi& F,
		const Eigen::MatrixXi& EV,
		const Eigen::MatrixXi& EF,
		const Eigen::SparseMatrix<double, Eigen::RowMajor>& basisCycles,
		const Eigen::VectorXd& indices,
		const int N,
		Eigen::VectorXd& adjustAngles,
		double &error)
	{
		Eigen::MatrixXd B1, B2, B3;
		igl::local_basis(V, F, B1, B2, B3);
		trivial_connection(V, F, EV, EF, B1, B2, basisCycles, indices, N, adjustAngles, error);
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
	//TODO: work with boundaries
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
		trivial_connection(V, F, EV, EF, basisCycles, indices, N, adjustAngles, error);
	}
}

    


#endif


