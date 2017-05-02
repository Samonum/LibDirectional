// Copyright (C) 2016 Daniele Panozzo <daniele.panozzo@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
#ifndef POLY_TO_RAW_H
#define POLY_TO_RAW_H

#include <Eigen/Geometry>
#include <Eigen/Sparse>
#include <Eigen/SparseCholesky>
#include <Eigen/Eigenvalues>
#include <unsupported/Eigen/Polynomials>
#include <igl/triangle_triangle_adjacency.h>
#include <igl/local_basis.h>
#include <iostream>

namespace directional
{
	// Returns a list of faces, vertices and colour values that can be used to draw a vector field 
	// on a mesh.
	// Inputs:
	//  V: #V X 3 vertex coordinates.
	//  F: #F by 3 face vertex indices.
	//  TT: #F by 3 Triangle-triangle adjecencies.
	//  B1, B2: #F by 3 matrices representing the local base of each face.
	//  N: The degree of the field..
	// Outputs:
	//  complex: Representation of the field as complex double
	IGL_INLINE void poly_to_raw
	(
		const Eigen::MatrixXd& B1,
		const Eigen::MatrixXd& B2,
		const Eigen::MatrixXcd& poly,
		const int N,
		Eigen::MatrixXd& raw)
	{
		raw.resize(B1.rows(), 3 * N);
		for (int f = 0; f < B1.rows(); f++)
		{
			// Find the roots of p(t) = (t - c0)^n using
			// https://en.wikipedia.org/wiki/Companion_matrix
			Eigen::MatrixXcd M = Eigen::MatrixXcd::Zero(N, N);
			for (int i = 1; i < N; ++i)
				M(i, i - 1) = std::complex<double>(1, 0);
			M.col(N - 1) = poly.row(f).transpose();
			Eigen::MatrixXcd roots = M.eigenvalues();
			for (int i = 0; i < N; i++)
			{
				std::complex<double> root = roots(i);
				raw.block<1, 3>(f, 3 * i) = B1.row(f) * root.real() + B2.row(f) * root.imag();
			}
		}
	}

	// Takes a complex field and extracts a representative vector.
	// on a mesh.
	// Inputs:
	//  V: #V X 3 vertex coordinates.
	//  F: #F by 3 face vertex indices.
	//  N: The degree of the field..
	// Outputs:
	//  complex: Representation of the field as complex double
	IGL_INLINE void poly_to_raw
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
