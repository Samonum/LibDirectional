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
#include <Eigen/Core>
#include <vector>
#include <cmath>


namespace directional
{    
    // Computes the adjustment angles to form a trivial connection according to given cone curvatures (or singularity indices) around basis cycles. In case the sum of curvature is not consistent with the topology, the system is solved in least squares and unexpected singularities may appear elsewhere. The output is the modification to the parallel transport.
    // Inputs:
    //  V: #V X 3 vertex coordinates
    //  F: #F by 3 face vertex indices
    //  EF: #E X 2 edges 2 faces indices
    //  basisCycles: #basisCycles X #E the oriented (according to EF) basis cycles around which the curvatures are measured
    //  the basis cycles must be arranged so that the first |V| are the vertex cycles (for instance, the result of igl::basis_cycles())
    //  indices: #basisCycles the index around each cycle. They should add up to N*Euler_characteristic of the mesh.
    //  N: the degree of the field. The curvature of a cycle is measured by (singIndex/N)*(2*pi) (can be negative)
    // Outputs:
    //  adjustAngles: the difference between the parallel transport and the modified one.
    //TODO: work with boundaries
	IGL_INLINE void trivial_connection(const Eigen::MatrixXd& V,
		const Eigen::MatrixXi& F,
		const Eigen::MatrixXi& EV,
		const Eigen::MatrixXi& EF,
		const Eigen::SparseMatrix<double, Eigen::RowMajor>& basisCycles,
		const Eigen::SparseVector<int>& indices,
		const int N,
		Eigen::VectorXd& adjustAngles)
	{
		using namespace Eigen;
		using namespace std;
		VectorXi rows = ArrayXi::Zero(basisCycles.rows());
		VectorXi columns = ArrayXi::Zero(basisCycles.cols());


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
		SparseVector<double> reducedIndices(rows[rows.size() - 1]);
		reducedCycles.reserve(VectorXi::Constant(columns[columns.size() - 1], 2));

		for (int k = 0; k < basisCycles.outerSize(); ++k)
			for (SparseMatrix<double, Eigen::RowMajor>::InnerIterator it(basisCycles, k); it; ++it)
				if (it.value())
                    reducedCycles.insert(rows(it.row()) - 1, columns(it.col()) - 1) = it.value();

		for (SparseVector<int>::InnerIterator it(indices); it; ++it)
			reducedIndices.insert(rows(it.row()) - 1) = (double)it.value();


		VectorXd VK;  //just for comparison
		igl::gaussian_curvature(V, F, VK);
		MatrixXd B1, B2, B3;
		igl::local_basis(V, F, B1, B2, B3);
		VectorXd edgeParallelAngleChange(columns(columns.size() - 1));  //the difference in the angle representation of edge i from EF(i,0) to EF(i,1)
		MatrixXd edgeVectors(columns(columns.size() - 1), 3);

		// Same as igl::parallel_transport_angles?
		for (int i = 0, j = 0; i < EF.rows(); i++) {
			if (columns(i) == 0 || (i > 0 && columns(i) == columns(i - 1)))
				continue;

			edgeVectors.row(j) = (V.row(EV(i, 1)) - V.row(EV(i, 0))).normalized();
			double x1 = edgeVectors.row(j).dot(B1.row(EF(i, 0)));
			double y1 = edgeVectors.row(j).dot(B2.row(EF(i, 0)));
			double x2 = edgeVectors.row(j).dot(B1.row(EF(i, 1)));
			double y2 = edgeVectors.row(j).dot(B2.row(EF(i, 1)));
			edgeParallelAngleChange(j) = atan2(y2, x2) - atan2(y1, x1);
			j++;
		}

		//TODO: reduce extra cycles to the holonomy
		VectorXd cycleHolonomy = reducedCycles*edgeParallelAngleChange;
		for (int i = 0; i < cycleHolonomy.size(); i++) {
			while (cycleHolonomy(i) >= M_PI) cycleHolonomy(i) -= 2.0*M_PI;
			while (cycleHolonomy(i) < -M_PI) cycleHolonomy(i) += 2.0*M_PI;
		}

		VectorXd cycleNewCurvature = reducedIndices*(2.0*M_PI / (double)N);
		//MatrixXd Test(cycleHolonomy.rows(),3); Test<<cycleHolonomy, VK, cycleNewCurvature;
		//cout<<"basis cycle differences"<<(cycleHolonomy-VK).lpNorm()<<endl;
		//SparseQR<SparseMatrix<double>, COLAMDOrdering<int> > solver;
		//solver.compute(basisCycles.transpose());
		//adjustAngles=solver.solve(-cycleHolonomy+cycleNewCurvature);
		//SparseMatrix<double> Q, R, P;
		//Q = SparseQR<SparseMatrix<double>, COLAMDOrdering<int> >(basisCycles).matrixQ();
		//R=solver.matrixR().transpose();

		//[Q, R, E] = qr(A', 0);
		//x = Q * (R' \ (E' \ b));

	   // cout<<"solver.colsPermutation().rows(), solver.colsPermutation().cols(): "<<solver.colsPermutation().rows()<<", "<<solver.colsPermutation().cols()<<endl;

		//VectorXd ETrhs=solver.colsPermutation().transpose()*(-cycleHolonomy+cycleNewCurvature);
		//VectorXd y =R.triangularView<Lower>().solve(ETrhs);
		//adjustAngles=Q*y;

		//x = A' * ((A*A')\ b);
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

		std::cout << "Error of adjustment angles computation: " << (reducedCycles*reducedAngles - (-cycleHolonomy + cycleNewCurvature)).lpNorm<Infinity>() << std::endl;
	}
}

    


#endif


