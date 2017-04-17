// Copyright (C) 2016 Amir Vaxman <avaxman@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
#ifndef ADJUSTMENT_TO_REPRESENTATIVE_H
#define ADJUSTMENT_TO_REPRESENTATIVE_H
#include <igl/igl_inline.h>
#include <igl/gaussian_curvature.h>
#include <igl/local_basis.h>
#include <igl/edge_topology.h>
#include <Eigen/Core>
#include <Eigen/Sparse>
#include <vector>
#include <cmath>

namespace directional
{
	// Takes a (trivial) connection, in the form of adjustment angles that encode the change from the parallel transport, 
	// and computes a unit-norm N-rosy that conforms to this connection. The adjustment angles must
	// respect the triviality constraints within the dgree N or the results would be unpredictable
	// Inputs:
	//  V: #V X 3 vertex coordinates
	//  F: #F by 3 face vertex indices
	//  EV: #E x 2 edges 2 vertices indices
	//  EF: #E X 2 edges 2 faces indices
    //  B1, B2: #F by 3 matrices representing the local base of each face.
	//  adjustAngles: #E angles that encode deviation from parallel transport
	//  N: the degree of the field
	//  globalRotation: The angle between the vector on the first face and its basis in radians
	// Outputs:
	//  representative: #F x 3 representative vectors within the tangent space (supporting plane) of the respetive face.

	IGL_INLINE void adjustment_to_representative(const Eigen::MatrixXd& V,
		const Eigen::MatrixXi& F,
		const Eigen::MatrixXi& EV,
		const Eigen::MatrixXi& EF,
		const Eigen::MatrixXd& B1,
		const Eigen::MatrixXd& B2,
		const Eigen::VectorXd& adjustAngles,
		const int N,
		double globalRotation,
		Eigen::MatrixXd& representative)
	{
		typedef std::complex<double> Complex;
		using namespace Eigen;
		using namespace std;

		Complex globalRot = exp(Complex(0, globalRotation));
		MatrixXcd edgeRep(EF.rows(), 2);
		for (int i = 0; i<EF.rows(); i++) {
			for (int j = 0; j<2; j++) {
				if (EF(i, j) == -1)
					continue;
				VectorXd edgeVector = (V.row(EV(i, 1)) - V.row(EV(i, 0))).normalized();
				edgeRep(i, j) = pow(Complex(edgeVector.dot(B1.row(EF(i, j))), edgeVector.dot(B2.row(EF(i, j)))), (double)N);
			}
		}

		SparseMatrix<Complex> aP1Full(EF.rows(), F.rows());
		SparseMatrix<Complex> aP1(EF.rows(), F.rows() - 1);
		vector<Triplet<Complex> > aP1Triplets, aP1FullTriplets;
		for (int i = 0; i<EF.rows(); i++) {
			if (EF(i, 0) == -1 || EF(i, 1) == -1)
				continue;

			aP1FullTriplets.push_back(Triplet<Complex>(i, EF(i, 0), conj(edgeRep(i, 0))*exp(Complex(0, (double)N*adjustAngles(i)))));
			aP1FullTriplets.push_back(Triplet<Complex>(i, EF(i, 1), -conj(edgeRep(i, 1))));
			if (EF(i, 0) != 0)
				aP1Triplets.push_back(Triplet<Complex>(i, EF(i, 0) - 1, conj(edgeRep(i, 0))*exp(Complex(0, (double)N*adjustAngles(i)))));
			if (EF(i, 1) != 0) 
				aP1Triplets.push_back(Triplet<Complex>(i, EF(i, 1) - 1, -conj(edgeRep(i, 1))));
		}
		aP1Full.setFromTriplets(aP1FullTriplets.begin(), aP1FullTriplets.end());
		aP1.setFromTriplets(aP1Triplets.begin(), aP1Triplets.end());
		VectorXcd torhs = VectorXcd::Zero(F.rows()); torhs(0) = globalRot;  //global rotation
		VectorXcd rhs = -aP1Full*torhs;

		Eigen::SimplicialLDLT<Eigen::SparseMatrix<Complex> > solver;
		solver.compute(aP1.adjoint()*aP1);
		VectorXcd complexPowerField(F.rows());
		complexPowerField(0) = globalRot;
		complexPowerField.tail(F.rows() - 1) = solver.solve(aP1.adjoint()*rhs);
		VectorXcd complexField = pow(complexPowerField.array(), 1.0 / (double)N);

		//std::cout << "Error of field integration: " << (aP1*complexPowerField.tail(F.rows() - 1) - rhs).lpNorm<Infinity>() << std::endl;
		//std::cout << "Full Error of field integration: " << (aP1Full*complexPowerField).lpNorm<Infinity>() << std::endl;

		//constructing representative
		representative.conservativeResize(F.rows(), 3);
		for (int i = 0; i < F.rows(); i++) {
			RowVector3d currVector = B1.row(i)*complexField(i).real() + B2.row(i)*complexField(i).imag();
			//std::cout<<"currVector: "<<currVector<<std::endl;
			representative.row(i) = currVector;
		}
	}


	// Takes a (trivial) connection, in the form of adjustment angles that encode the change from the parallel transport, 
	// and computes a unit-norm N-rosy that conforms to this connection. The adjustment angles must
	// respect the triviality constraints within the dgree N or the results would be unpredictable
	// Inputs:
	//  V: #V X 3 vertex coordinates
	//  F: #F by 3 face vertex indices
	//  EV: #E x 2 edges 2 vertices indices
	//  EF: #E X 2 edges 2 faces indices
	//  adjustAngles: #E angles that encode deviation from parallel transport
	//  N: the degree of the field
	//  globalRotation: The angle between the vector on the first face and its basis in radians
	// Outputs:
	//  representative: #F x 3 representative vectors within the tangent space (supporting plane) of the respetive face.
	IGL_INLINE void adjustment_to_representative(const Eigen::MatrixXd& V,
		const Eigen::MatrixXi& F,
		const Eigen::MatrixXi& EV,
		const Eigen::MatrixXi& EF,
		const Eigen::VectorXd& adjustAngles,
		const int N,
		double globalRotation,
		Eigen::MatrixXd& representative)
	{
		Eigen::MatrixXd B1, B2, x;
		igl::local_basis(V, F, B1, B2, x);
		directional::adjustment_to_representative(V, F, EV, EF, B1, B2, adjustAngles, N, globalRotation, representative);
	}

	// Takes a (trivial) connection, in the form of adjustment angles that encode the change from the parallel transport, 
	// and computes a unit-norm N-rosy that conforms to this connection. The adjustment angles must
	// respect the triviality constraints within the dgree N or the results would be unpredictable
	// Inputs:
	//  V: #V X 3 vertex coordinates
	//  F: #F by 3 face vertex indices
	//  adjustAngles: #E angles that encode deviation from parallel transport
	//  N: the degree of the field
	//  globalRotation: The angle between the vector on the first face and its basis in radians
	// Outputs:
	//  representative: #F x 3 representative vectors within the tangent space (supporting plane) of the respetive face.
	IGL_INLINE void adjustment_to_representative(const Eigen::MatrixXd& V,
		const Eigen::MatrixXi& F,
		const Eigen::VectorXd& adjustAngles,
		const int N,
		double globalRotation,
		Eigen::MatrixXd& representative)
	{
		Eigen::MatrixXi EV, x, EF;
		igl::edge_topology(V, F, EV, x, EF);
		directional::adjustment_to_representative(V, F, EV, EF, adjustAngles, N, globalRotation, representative);
	}
}

#endif
