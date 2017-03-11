#ifndef REPRESENTATIVE_TO_ADJUSTMENT_H
#define REPRESENTATIVE_TO_ADJUSTMENT_H
#include <directional/adjustment_to_representative.h>
#include <directional/representative_to_raw.h>
#include <igl/local_basis.h>
#include <igl/parallel_transport_angles.h>
#include <Eigen/Core>

namespace directional
{
	// Computes the raw vector field given the adjustment angles.
	// Inputs::
	//  V: #V X 3 vertex coordinates.
	//  F: #F by 3 face vertex indices.
	//  EV: #E x 2 edges 2 vertices indices.
	//  EF: #E X 2 edges 2 faces indices.
	//  norm: #F normals for each face.
	//  adjustAngles: #E angles that encode deviation from parallel transport.
	//  N: the degree of the field.
	//  globalRotation: The angle between the vector on the first face and its basis in radians.
	// Outputs:
	//  raw: #F by 3*N matrix with all N explicit vectors of each directional.
	IGL_INLINE void representative_to_adjustment(const Eigen::MatrixXd& V,
		const Eigen::MatrixXi& F,
		const Eigen::MatrixXi& EV,
		const Eigen::MatrixXi& EF,
		const Eigen::MatrixXi& FE,
		const Eigen::MatrixXd &representative,
		int N,
		Eigen::VectorXd& adjustAngles,
		Eigen::VectorXi& matching,
		double &globalRotation)
	{
		Eigen::MatrixXd B1, B2, B3;
		Eigen::VectorXd local(F.rows()), parallel;
		igl::local_basis(V, F, B1, B2, B3);
		adjustAngles.resize(EV.rows(), 1);
		matching.resize(EV.rows(), 1);

		//Obtain angles for parallel transport
		igl::parallel_transport_angles(V, F, B3, EF, FE, parallel);

		//Calculate translations in local base
		for (size_t i = 0; i < F.rows(); i++)
		{
			double dot = representative.row(i).dot(B1.row(i));
			Eigen::Matrix3d mat;
			mat << B1.row(i),
				representative.row(i),
				B3.row(i);
			double det = mat.determinant();
			local(i) = atan2(det, dot);
		}

		globalRotation = local(0);

		for (size_t i = 0; i < EF.rows(); i++)
		{
			adjustAngles(i) = local(EF(i,1)) - local(EF(i,0)) -parallel(i);
			matching(i) = adjustAngles(i) / (igl::PI / N);
			adjustAngles(i) -= (igl::PI / N) * matching(i);
		}



		/*
		Eigen::VectorXd local(F.rows());
		int end = 1;
		//Calculate local rotations
		for (size_t i = 0; i < F.rows(); i++)
		{
			int cur = queue[i];
			for (int j = 0; j < 3; j++)
			{
				if (flags[TT(cur, j)])
					continue;
				Eigen::MatrixXd edge = V.row(F(cur, (j + 1) % 3)) - V.row(F(cur, j));
				edge.normalize();
				double T1 = acos(edge.dot(B1.row(cur)));
				double T2 = acos(edge.dot(B1.row(TT(cur, j))));
				local[TT(cur, j)] = local[cur] - T1 + T2 + (E(TTi(cur, j), 0) == cur ? 1 : -1) * trivial[TTi(cur, j)];
				for (int k = 0; k < 3; k++)
				{
					int x = TT(TT(cur, j), k);
					if (!flags[x])
					{
						flags[x] = 1;
						queue[cur++] = x;
					}
				}
			}
		}

		//Calculate representative vectors
		representative.resize(F.rows(), 3);
		for (size_t i = 0; i < F.rows(); i++)
			representative.row(i) << (Eigen::AngleAxisd(local[i], B3.row(i)).toRotationMatrix() * B1.row(i).transpose()).transpose();
	*/
	}
}

#endif
