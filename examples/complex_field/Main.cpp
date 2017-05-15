#include <iostream>
#include <directional/drawable_field.h>
#include <directional/complex_field.h>
#include <directional/complex_to_representative.h>
#include <directional/complex_to_raw.h>
#include <Eigen/Core>
#include <igl/viewer/Viewer.h>
#include <igl/read_triangle_mesh.h>
#include <igl/per_face_normals.h>
#include <igl/triangle_triangle_adjacency.h>
#include <igl/unproject_onto_mesh.h>
#include <igl/boundary_loop.h>
#include "Main.h"


Eigen::VectorXi cIDs;
Eigen::MatrixXi F, fieldF, meshF;
Eigen::MatrixXd V, C, meshV, meshC, fieldV, fieldC, representative, cValues;
Eigen::MatrixXcd complex;
igl::viewer::Viewer viewer;


int N = 4;
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
	// Calculate the field
	directional::complex_field(meshV, meshF, cIDs, cValues, N, complex);

	// Convert it so it can be drawn
	directional::complex_to_representative(meshV, meshF, complex, N, representative);

	// Normalize if wanted
	if (normalized)
		representative.rowwise().normalize();

	// Calculate the vectors, faces and colors of the field representation
	directional::drawable_field(meshV, meshF, representative, Eigen::RowVector3d(0, 0, 1), N, 0, fieldV, fieldF, fieldC);
	meshC = Eigen::RowVector3d(1, 1, 1).replicate(meshF.rows(), 1);

	// Merge meshes
	for (int i = 0; i < cIDs.rows(); i++)
		meshC.row(cIDs(i)) = Eigen::RowVector3d(1, 0, 0);
	
	C.resize(meshC.rows() + fieldC.rows(), 3);
	C << meshC, fieldC;

	ConcatMeshes(meshV, meshF, fieldV, fieldF, V, F);

	// Update the viewer
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
	// Toggle field drawing for easier rotation
	case 'D':
		drag = !drag;
		break;

	// Reset the constraints
	case 'R':
		cIDs.resize(0);
		cValues.resize(0, 6);
		draw_field();
		break;

	// Toggle normalization
	case 'N':
		normalized = !normalized;
		draw_field();
		break;
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
			cValues.conservativeResize(cValues.rows() + 1, 3);
		}

		// Calculate direction from the center of the face to the mouse
		cValues.row(i) = 
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
	std::cout <<
		"  R       Reset the constraints" << std::endl <<
		"  N       Toggle field normalization" << std::endl <<
		"  L-bttn  place constraint" << std::endl <<
		"  D       Toggle constraint placement" << std::endl;

	// Load mesh
	igl::readOBJ("../../data/torus.obj", meshV, meshF);

	cIDs.resize(0);
	cValues.resize(0, 3);


	draw_field();
	viewer.launch();
}
