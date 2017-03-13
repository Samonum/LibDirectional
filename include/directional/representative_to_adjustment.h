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
		Eigen::VectorXd local, parallel;
		igl::local_basis(V, F, B1, B2, B3);
		adjustAngles.resize(EV.rows(), 1);
		matching.resize(EV.rows(), 1);

		local.resize(F.rows());

		//Obtain angles for parallel transport
		igl::parallel_transport_angles(V, F, B3, EF, FE, parallel);

		//Calculate translations in local base
		for (int i = 0; i < F.rows(); i++)
		{
			double dot = representative.row(i).dot(B1.row(i));
			Eigen::Matrix3d mat;
			mat << representative.row(i),
				B1.row(i),
				B3.row(i);
			double det = mat.determinant();
			local[i] = atan2(det, dot);
		}

		globalRotation = local(0);

		for (int i = 0; i < EF.rows(); i++)
		{
			adjustAngles(i) = local(EF(i,0)) - local(EF(i,1)) -parallel(i);
			matching(i) =(int)(adjustAngles(i) / (2*igl::PI / N));
			adjustAngles(i) -= (2 * igl::PI / N) * matching(i);
			if (adjustAngles(i) > igl::PI/N)
			{
				adjustAngles(i) -= 2 * igl::PI/N;
				matching(i) -= N;
			}
		}
	}
}

#endif
