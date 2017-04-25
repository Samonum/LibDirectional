#ifndef CYCLE_HOLONOMY_H
#define CYCLE_HOLONOMY_H
#include <igl/igl_inline.h>
#include <igl/gaussian_curvature.h>
#include <igl/local_basis.h>
#include <igl/parallel_transport_angles.h>
#include <igl/edge_topology.h>
#include <igl/is_border_vertex.h>
#include <igl/gaussian_curvature.h>
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
	IGL_INLINE void cycle_holonomy(const Eigen::MatrixXd& V,
		const Eigen::MatrixXi& F,
		const Eigen::MatrixXi& EV,
		const Eigen::MatrixXi& EF,
		const Eigen::MatrixXd& B1,
		const Eigen::MatrixXd& B2,
		const Eigen::SparseMatrix<double, Eigen::RowMajor>& basisCycles,
		Eigen::VectorXd& cycleHolonomy)
	{
		using namespace Eigen;
		using namespace std;

		vector<bool> border = igl::is_border_vertex(V, F);
		// Start of to be replaced section
		VectorXd edgeParallelAngleChange(basisCycles.cols());  //the difference in the angle representation of edge i from EF(i,0) to EF(i,1)
																		//MatrixXd edgeVectors(columns(columns.size() - 1), 3);

		for (int i = 0, j = 0; i < EF.rows(); i++) {
			//skip border edges
			if (EF(i, 0) == -1 || EF(i, 1) == -1)
				continue;

			RowVectorXd edgeVectors = (V.row(EV(i, 1)) - V.row(EV(i, 0))).normalized();
			double x1 = edgeVectors.dot(B1.row(EF(i, 0)));
			double y1 = edgeVectors.dot(B2.row(EF(i, 0)));
			double x2 = edgeVectors.dot(B1.row(EF(i, 1)));
			double y2 = edgeVectors.dot(B2.row(EF(i, 1)));
			edgeParallelAngleChange(i) = atan2(y2, x2) - atan2(y1, x1);
		}


		cycleHolonomy = basisCycles*edgeParallelAngleChange;
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
			if (!border[i])
				cycleHolonomy(i) = gc(i);

		//Sum holonomy for boundaries
		int i = 0, added = 0, total = count(border.begin(), border.end(), true);
		while (added < total)
		{
			double temp = cycleHolonomy(V.rows() - total + i);
			cycleHolonomy(V.rows() - total + i) = 0;
			//Loop over boundary cycles
			for (SparseMatrix<double, Eigen::RowMajor>::InnerIterator it(basisCycles, V.rows() + i); it; ++it)
				if (it.value())
				{
					int b = (it.value() + 1) / 2;
					if (border[EV(it.col(), b)])
					{
						// Substract Pi, as igl::gaussian_curvature calculates gc as 2 * PI - sum(angles)
						cycleHolonomy(V.rows() - total + i) += gc(EV(it.col(), b)) - igl::PI;
						cout << "Curvature at border vertex: " << gc(EV(it.col(), b)) - igl::PI << endl;
						//Ensure vertices won't be counted twice
						border[EV(it.col(), b)] = false;
						added++;
					}
				}
			cout << "Old holonomy: " << temp << " " << "    New holonomy: " << cycleHolonomy(V.rows() - total + i) << endl;
			i++;
		}
	}
}




#endif