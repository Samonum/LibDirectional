#ifndef DRAW_CYCLES
#define DRAW_CYCLES
#include <igl/igl_inline.h>
#include <igl/barycenter.h>
#include <igl/per_face_normals.h>
#include <directional/representative_to_raw.h>
#include <directional/point_spheres.h>
#include <directional/line_cylinders.h>
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
	void IGL_INLINE draw_cycles(const Eigen::MatrixXi &EF, 
		const Eigen::SparseMatrix<double, Eigen::RowMajor>& basisCycleMat,
		const Eigen::VectorXd &cycleColor,
		int index,
		Eigen::MatrixXd &color)
	{
		for (Eigen::SparseMatrix<double, Eigen::RowMajor>::InnerIterator it(basisCycleMat, index); it; ++it)
		{
			color.row(EF(it.col(), 0)) = cycleColor;
			color.row(EF(it.col(), 1)) = cycleColor;
		}
	}

}

#endif