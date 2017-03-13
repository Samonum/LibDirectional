#ifndef ADJUSTMENT_TO_RAW_H
#define ADJUSTMENT_TO_RAW_H
#include <directional/adjustment_to_representative.h>
#include <directional/representative_to_raw.h>


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
	IGL_INLINE void adjustment_to_raw(const Eigen::MatrixXd& V,
		const Eigen::MatrixXi& F,
		const Eigen::MatrixXi& EV,
		const Eigen::MatrixXi& EF,
		const Eigen::MatrixXd& norm,
		const Eigen::MatrixXd& adjustAngles,
		int N,
		double globalRotation,
		Eigen::MatrixXd& raw)
	{
		Eigen::MatrixXd representative;
		adjustment_to_representative(V, F, EV, EF, adjustAngles, N, globalRotation, representative);
		representative_to_raw(norm, representative, N, raw);
	}
}

#endif
