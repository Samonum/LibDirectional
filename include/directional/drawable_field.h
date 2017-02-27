#ifndef DRAWABLE_FIELD
#define DRAWABLE_FIELD
#include <igl/igl_inline.h>
#include <igl/barycenter.h>
#include <directional/representative_to_raw.h>
#include <directional/point_spheres.h>
#include <directional/line_cylinders.h>
#include <Eigen/Core>

namespace directional
{
	void IGL_INLINE drawable_field(const Eigen::MatrixXd &V, 
		const Eigen::MatrixXi &F, 
		const Eigen::MatrixXd &field, 
		int N, 
		const Eigen::MatrixXd &color, 
		Eigen::MatrixXd &fieldV, 
		Eigen::MatrixXi &fieldF, 
		Eigen::MatrixXd &fieldC)
	{
		
		Eigen::MatrixXd rawField, barycenters, vectorColors, P1, P2;
		igl::barycenter(V, F, barycenters);

		P1.resize(F.rows() * N, 3);
		P2.resize(F.rows() * N, 3);
		vectorColors.resize(F.rows() * N, 3);
		P1 = barycenters.replicate(N, 1);
		
		if (field.cols() == 3)
			representative_to_raw(V, F, field, rawField, N);
		else
			rawField = field;

		for (int i = 0; i < N; i++)
			P2.middleRows(F.rows()*i, F.rows()) = rawField.middleCols(3*i, 3);

		P2.array() *= .01;
		P2 += P1;
		if (color.rows() == 1)
			vectorColors = color.replicate(P1.rows(), 1);
		else if (color.rows() == F.rows())
			vectorColors = color.replicate(N, 1);
		else
			vectorColors = color;

		line_cylinders(P1, P2, .001, vectorColors, 4, false, fieldV, fieldF, fieldC);
	}

}

#endif