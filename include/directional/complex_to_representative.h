// Copyright (C) 2016 Daniele Panozzo <daniele.panozzo@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
#ifndef COMPLEX_TO_REPRESENTATIVE_H
#define COMPLEX_TO_REPRESENTATIVE_H

#include <Eigen/Geometry>
#include <Eigen/Sparse>
#include <Eigen/SparseCholesky>
#include <Eigen/Eigenvalues>
#include <igl/local_basis.h>
#include <iostream>

namespace directional
{
	// Takes a complex field and extracts a representative vector.
	// Inputs:
	//  B1, B2: #F by 3 matrices representing the local base of each face.
	//  complex: Representation of the field as complex double.
	//  N: The degree of the field.
	// Outputs:
	//  representative: #F x 3 representative vectors within the tangent space (supporting plane) of the respetive face.
	IGL_INLINE void complex_to_representative
	(
		const Eigen::MatrixXd& B1,
		const Eigen::MatrixXd& B2,
		const Eigen::MatrixXcd& complex,
		const int N,
	    Eigen::MatrixXd& representative
	)
	{
		// Convert the interpolated polyvector into Euclidean vectors
		representative.resize(B1.rows(), 3);
		for (int f = 0; f < B1.rows(); ++f)
		{
			// Find the roots of p(t) = (t - c0)^n using
			// https://en.wikipedia.org/wiki/Companion_matrix
			Eigen::MatrixXcd M = Eigen::MatrixXcd::Zero(N, N);
			for (int i = 1; i < N; ++i)
				M(i, i - 1) = std::complex<double>(1, 0);
			M(0, N - 1) = complex(f, 0);
			std::complex<double> root = M.eigenvalues()(0);
			representative.row(f) = B1.row(f) * root.real() + B2.row(f) * root.imag();
		}
	}


	// Takes a complex field and extracts a representative vector.
	// on a mesh.
	// Inputs:
	//  V: #V X 3 vertex coordinates.
	//  F: #F by 3 face vertex indices.
	//  complex: Representation of the field as complex double.
	//  N: The degree of the field.
	// Outputs:
	//  representative: #F x 3 representative vectors within the tangent space (supporting plane) of the respetive face.
	IGL_INLINE void complex_to_representative
	(
		const Eigen::MatrixXd& V,
		const Eigen::MatrixXi& F,
		const Eigen::MatrixXcd& complex,
		const int N,
		Eigen::MatrixXd& representative
	)
	{
		Eigen::MatrixXd B1, B2, x;
		igl::local_basis(V, F, B1, B2, x);
		complex_to_representative(B1, B2, complex, N, representative);
	}
}
#endif