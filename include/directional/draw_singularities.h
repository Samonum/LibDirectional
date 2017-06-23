#ifndef DRAW_SINGULARITIES
#define DRAW_SINGULARITIES
#include <cmath>
#include <igl/igl_inline.h>
#include <igl/barycenter.h>
#include <igl/per_face_normals.h>
#include <directional/representative_to_raw.h>
#include <directional/point_spheres.h>
#include <Eigen/Core>
#include <igl/avg_edge_length.h>


namespace directional
{
	// Returns a list of faces, vertices and colour values that can be used to draw a vector field on a mesh.
	// Inputs:
	//  V: #V X 3 vertex coordinates.
	//  F: #F by 3 face vertex indices.
	//  field: Either a representative or a raw vector field.
	//  color: An array of either 1 by 3 color values for each vector, #F by 3 colors for each individual directional or #F*N by 3 colours for each individual vector, ordered by #F times vector 1, followed by #F times vector 2 etc.
	//  N: The degree of the field.
	//  width: The width of each vector in the vector field.
	//  length: The width of each vector in the vector field, respective to each vector's length.
	//  colorPerVertex: The coloring mode used to draw the meash.
	// Outputs:
	//  fieldV: The vertices of the field.
	//  fieldF: The faces of the field.
	//  fieldC: The colors of  the field.
	void IGL_INLINE draw_singularities(const Eigen::MatrixXd& V,
		const Eigen::VectorXi& singularities,
		const Eigen::MatrixXd& positiveColors,
		const Eigen::MatrixXd& negativeColors,
		const double r,
		Eigen::MatrixXd &singV, 
		Eigen::MatrixXi &singF,
		Eigen::MatrixXd &singC)
	{
		std::vector<int> indices;
		for (int i = 0; i < V.rows(); i++)
		{
			if (singularities(i))
				indices.push_back(i);
		}

		Eigen::MatrixXd points(indices.size(), 3);
		Eigen::MatrixXd colors(indices.size(), 3);
		for (int i = 0; i < indices.size(); i++)
		{
			points.row(i) = V.row(indices[i]);
			if (singularities(indices[i]) > 0)
				colors.row(i) = positiveColors.row(std::min(std::round(singularities(indices[i]))-1, (double)positiveColors.rows()-1));
			else
				colors.row(i) = negativeColors.row(std::min(std::abs(std::round(singularities(indices[i])))-1, (double)negativeColors.rows()-1));
		}
		directional::point_spheres(points, r, colors, 16, false, singV, singF, singC);
	}

}

#endif