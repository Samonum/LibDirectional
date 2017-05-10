#include <iostream>
#include <directional/drawable_field.h>
#include <directional/complex_field.h>
#include <directional/complex_to_representative.h>
#include <directional/complex_to_raw.h>
#include <directional/poly_to_raw.h>
#include <directional/poly_vector.h>
#include <Eigen/Core>
#include <igl/viewer/Viewer.h>
#include <igl/read_triangle_mesh.h>
#include <igl/per_face_normals.h>
#include <igl/unproject_onto_mesh.h>
#include <igl/boundary_loop.h>
#include "Main.h"


Eigen::VectorXi cIDs;
Eigen::MatrixXi F, fieldF, meshF;
Eigen::MatrixXd V, C, meshV, meshC, fieldV, fieldC, raw, cValues;
Eigen::MatrixXcd complex;
igl::viewer::Viewer viewer;

//Degree of the field
int N = 2;


//User input variables
int cur = 0;
bool drag = false;
bool normalized = false;

void ConcatMeshes(const Eigen::MatrixXd &VA, const Eigen::MatrixXi &FA, const Eigen::MatrixXd &VB, const Eigen::MatrixXi &FB, Eigen::MatrixXd &V, Eigen::MatrixXi &F)
{
	V.resize(VA.rows() + VB.rows(), VA.cols());
	V << VA, VB;
	F.resize(FA.rows() + FB.rows(), FA.cols());
	F << FA, (FB.array() + VA.rows());
}

void draw_field()
{
	directional::poly_vector(meshV, meshF, cIDs, cValues, N, complex);
	directional::poly_to_raw(meshV, meshF, complex, N, raw);
	if (normalized)
		for(int n = 0; n < N; n++)
			raw.middleCols(n*3, 3).rowwise().normalize();
	directional::drawable_field(meshV, meshF, raw, Eigen::RowVector3d(0, 0, 1), N, false, fieldV, fieldF, fieldC);
	meshC = Eigen::RowVector3d(1, 1, 1).replicate(meshF.rows(), 1);

	for (int i = 0; i < cIDs.rows(); i++)
		meshC.row(cIDs(i)) = Eigen::RowVector3d(1, 0, 0);
	
	C.resize(meshC.rows() + fieldC.rows(), 3);
	C << meshC, fieldC;

	ConcatMeshes(meshV, meshF, fieldV, fieldF, V, F);
	viewer.data.clear();
	viewer.data.set_face_based(true);
	viewer.data.set_mesh(V, F);
	viewer.data.set_colors(C);
}


bool key_down(igl::viewer::Viewer& viewer, int key, int modifiers)
{
	int borders;
	switch (key)
	{
	case '1':
		cur = 0;
		break;
	case '2':
		cur = std::max(1, N - 1);
		break;
	case '3':
		cur = std::max(2, N - 1);
		break;
	case '4':
		cur = std::max(3, N - 1);
		break;
	case '5':
		cur = std::max(4, N - 1);
		break;
	case '6':
		cur = std::max(5, N - 1);
		break;
	case 'C':
		draw_field();
		break;
	case 'D':
		drag = !drag;
		break;
	case 'R':
		cIDs.resize(0);
		cValues.resize(0, 6);
		draw_field();
		break;
	case 'N':
		normalized = !normalized;
		draw_field();
		break;
	/*case 'S':
		if (directional::write_trivial_field("test", meshV, meshF, indices, N, 0))
			std::cout << "Saved mesh" << std::endl;
		else
			std::cout << "Unable to save mesh. Error: " << errno << std::endl;
		break;
	case 'L':
		double x;
		directional::read_trivial_field("test", meshV, meshF, indices, N, x);
		update_mesh();
		calculate_field();
		draw_field();
		break;*/
	}
	return true;
}

//Select vertices using the mouse
bool mouse_down(igl::viewer::Viewer& viewer, int key, int modifiers)
{
	if (drag || key != 0)
		return false;
	int fid;
	Eigen::Vector3d bc;

	// Cast a ray in the view direction starting from the mouse position
	double x = viewer.current_mouse_x;
	double y = viewer.core.viewport(3) - viewer.current_mouse_y;
	if (igl::unproject_onto_mesh(Eigen::Vector2f(x, y), viewer.core.view * viewer.core.model,
		viewer.core.proj, viewer.core.viewport, meshV, meshF, fid, bc))
	{
		int i;
		for (i = 0; i < cIDs.rows(); i++)
			if (cIDs(i) == fid)
				break;
		if (i == cIDs.rows())
		{
			cIDs.conservativeResize(cIDs.rows() + 1);
			cIDs(i) = fid;
			cValues.conservativeResize(cValues.rows() + 1, 3*N);
			cValues.row(i).fill(0);
		}
		cValues.block<1, 3>(i, cur*3) =
			 (meshV.row(meshF(fid, 0)) * bc(0) + 
			 meshV.row(meshF(fid, 1)) * bc(1) + 
			 meshV.row(meshF(fid, 2)) * bc(2) - 
			(meshV.row(meshF(fid, 0)) + 
				meshV.row(meshF(fid, 1)) + 
				meshV.row(meshF(fid, 2))) / 3).normalized();
		draw_field();
		return true;
	}
	return false;
};

int main()
{
	viewer.callback_key_down = &key_down;
	viewer.callback_mouse_down = &mouse_down;
	igl::readOBJ("../../data/torus.obj", meshV, meshF);

	cIDs.resize(0);
	cValues.resize(0, 3*N);


	draw_field();
	viewer.launch();
}
