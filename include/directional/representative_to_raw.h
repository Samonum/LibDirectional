#ifndef REPRESENTATIVE_TO_RAW_H
#define REPRESENTATIVE_TO_RAW_H
#include <igl/per_face_normals.h>
#include <igl/igl_inline.h>
#include <igl/PI.h>
#include <Eigen/Core>
#include <Eigen/Geometry>
#include <iostream>


namespace directional
{
	// Computes the raw vector field given a set of representative vectors.
	// Inputs:
	//  norm: #F by 3 coordinates of the normals of each face.
	//  representative: #F by 3 coordinates of representative vectors.
	//  N: the degree of the field.
	// Outputs:
	//  raw: #F by 3*N matrix with all N explicit vectors of each directional.
	IGL_INLINE void representative_to_raw(const Eigen::MatrixXd& norm,
		const Eigen::MatrixXd& representative,
		const int N,
		Eigen::MatrixXd& raw)
	{
		raw.resize(representative.rows(), 3 * N);

		for (int i = 0; i < representative.rows(); i++)
		{
			raw.block<1, 3>(i, 0) << representative.row(i);

			Eigen::MatrixXd rot = Eigen::AngleAxisd(2. / N*igl::PI,
				norm.row(i)).toRotationMatrix();
			for (int j = 1; j < N; j++)
				raw.block<1, 3>(i, j * 3) << (rot*raw.block<1, 3>(i, (j - 1) * 3).transpose()).transpose();
		}
	}

	// Computes the raw vector field given a set of representative vectors.
	// This version recalculates the face normals every time it's called.
	// Inputs:
	//  V: #V X 3 vertex coordinates.
	//  F: #F by 3 face vertex indices.
	//  representative: #F by 3 coordinates of representative vectors.
	//  N: the degree of the field.
	// Outputs:
	//  raw: #F by 3*N matrix with all N explicit vectors of each directional.
	IGL_INLINE void representative_to_raw(const Eigen::MatrixXd& V,
		const Eigen::MatrixXi& F,
		const Eigen::MatrixXd& representative,
		const int N,
		Eigen::MatrixXd& raw)
	{
		Eigen::MatrixXd norm;
		igl::per_face_normals(V, F, norm);
		representative_to_raw(norm, representative, N, raw);
	}
}

#endif
