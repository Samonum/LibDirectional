#ifndef DRAWABLE_FIELD
#define DRAWABLE_FIELD
#include <igl/igl_inline.h>
#include <igl/barycenter.h>
#include <igl/per_face_normals.h>
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
		Eigen::MatrixXd normals;
		double width = 0.001;
		double length = .01;
		igl::per_face_normals(V, F, normals);

		Eigen::MatrixXd rawField, barycenters, vectorColors, P1, P2;
		igl::barycenter(V, F, barycenters);

		P1.resize(F.rows() * N, 3);
		P2.resize(F.rows() * N, 3);
		vectorColors.resize(F.rows() * N, 3);
		
		if (field.cols() == 3)
			representative_to_raw(normals, field, rawField, N);
		else
			rawField = field;

		normals.array() *= width;
		barycenters += normals;
		P1 = barycenters.replicate(N, 1);

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

		Eigen::MatrixXd Vc, Cc, Vs, Cs;
		Eigen::MatrixXi Fc, Fs;
		line_cylinders(P1, P2, width, vectorColors, 6, false, Vc, Fc, Cc);
		
		point_spheres(barycenters, width*2, vectorColors.topRows(barycenters.rows()), 10, false, Vs, Fs, Cs);

		Fs.array() += Fc.rows();

		fieldV.resize(Vc.rows() + Vs.rows(), Vc.cols());
		fieldC.resize(Cc.rows() + Cs.rows(), Cc.cols());
		fieldF.resize(Fc.rows() + Fs.rows(), Fc.cols());
		fieldV << Vc, Vs;
		fieldC << Cc, Cs;
		fieldF << Fc, Fs;
	}

}

#endif