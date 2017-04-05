#include <iostream>
#include <directional/drawable_field.h>
#include <directional/draw_singularities.h>
#include <directional/dual_cycles.h>
#include <directional/trivial_connection.h>
#include <directional/adjustment_to_representative.h>
#include <directional/draw_cycles.h>
#include <Eigen/Core>
#include <igl/viewer/Viewer.h>
#include <igl/read_triangle_mesh.h>
#include <igl/per_face_normals.h>
#include <igl/edge_topology.h>
#include <igl/unproject_onto_mesh.h>
#include "Main.h"


Eigen::MatrixXi F, fieldF, meshF, singF, EV, FE, EF;
Eigen::MatrixXd V, C, meshV, meshC, fieldV, fieldC, singV, singC, norm, representative;
Eigen::VectorXd adjustmentField, other;
Eigen::VectorXi match; 
Eigen::VectorXd indices;
Eigen::SparseMatrix<double, Eigen::RowMajor> cycles;
igl::viewer::Viewer viewer;

Eigen::MatrixXd positiveIndices(4, 3),
				negativeIndices(4,3);

int N = 4;
int ring1 = 265, ring2 = 10;

int genusLoop;

bool drag = false;

void ConcatMeshes(const Eigen::MatrixXd &VA, const Eigen::MatrixXi &FA, const Eigen::MatrixXd &VB, const Eigen::MatrixXi &FB, Eigen::MatrixXd &V, Eigen::MatrixXi &F)
{
	V.resize(VA.rows() + VB.rows(), VA.cols());
	V << VA, VB;
	F.resize(FA.rows() + FB.rows(), FA.cols());
	F << FA, (FB.array() + VA.rows());
}

void drawField()
{
	double e;
	directional::trivial_connection(meshV, meshF, EV, EF, cycles, indices, N, adjustmentField, e);
	std::cout << "error: " << e << std::endl;
	directional::adjustment_to_representative(meshV, meshF, EV, EF, adjustmentField, N, 0, representative);

	Eigen::MatrixXd color;
	color.resize(representative.rows()*N, 3);
	color << Eigen::RowVector3d(0, 1, 0).replicate(representative.rows(), 1),
		Eigen::RowVector3d(0, 0, 1).replicate((N - 1) *representative.rows(), 1);
	directional::drawable_field(meshV, meshF, representative, color, N, false, fieldV, fieldF, fieldC);

	meshC = Eigen::RowVector3d(1, 1, 1).replicate(meshF.rows(), 1);
	directional::draw_cycles(EF, cycles, Eigen::Vector3d(1, 0, 0), ring1, meshC);
	directional::draw_cycles(EF, cycles, Eigen::Vector3d(0, 1, 0), ring2, meshC);

	Eigen::MatrixXd spheres;
	directional::draw_singularities(meshV, indices, positiveIndices, negativeIndices, .015, singV, singF, singC);

	Eigen::MatrixXd a;
	Eigen::MatrixXi b;
	ConcatMeshes(meshV, meshF, fieldV, fieldF, a, b);
	ConcatMeshes(a, b, singV, singF, V, F);

	C.resize(F.rows(), 3);
	C << meshC, fieldC, singC;
	viewer.data.clear();
	viewer.data.set_face_based(true);
	viewer.data.set_mesh(V, F);
	viewer.data.set_colors(C);
}

void draw_rings()
{
	//Clear out the mesh
	meshC = Eigen::RowVector3d(1, 1, 1).replicate(meshF.rows(), 1);
	//Draw the new cycles
	directional::draw_cycles(EF, cycles, Eigen::Vector3d(1, 0, 0), ring1, meshC);
	directional::draw_cycles(EF, cycles, Eigen::Vector3d(0, 1, 0), ring2, meshC);
	C << meshC, fieldC, singC;
	viewer.data.set_colors(C);
}

bool key_down(igl::viewer::Viewer& viewer, int key, int modifiers)
{
	switch (key)
	{
	case '-':
	case '_':
		indices(ring1)--;
		indices(ring2)++;
		drawField();
		break;
	case '+':
	case '=':
		indices(ring1)++;
		indices(ring2)--;
		drawField();
		break;
	case 'D':
		drag = !drag;
		break;
	case '<':
	case ',':
		indices(genusLoop)--;
		drawField();
		break;
	case '>':
	case '.':
		indices(genusLoop)++;
		drawField();
		break;
	case '1':
		ring1 = cycles.rows() - 3;
		draw_rings();
		break;
	case '9':
		genusLoop = cycles.rows() - 2;
		//Clear out the mesh
		meshC = Eigen::RowVector3d(1, 1, 1).replicate(meshF.rows(), 1);
		//Draw the generator cycle
		directional::draw_cycles(EF, cycles, Eigen::Vector3d(0, 0, 1), genusLoop, meshC);
		C << meshC, fieldC, singC;
		viewer.data.set_colors(C);
		break;
	case '0':
		genusLoop = cycles.rows() - 1;
		//Clear out the mesh
		meshC = Eigen::RowVector3d(1, 1, 1).replicate(meshF.rows(), 1);
		//Draw the generator cycle
		directional::draw_cycles(EF, cycles, Eigen::Vector3d(0, 0, 1), genusLoop, meshC);
		C << meshC, fieldC, singC;
		viewer.data.set_colors(C);
		break;
	}
	return true;
}

//Select vertices using the mouse
bool mouse_down(igl::viewer::Viewer& viewer, int key, int modifiers)
{
	if (drag || (key != 0 && key != 2))
		return false;
	int fid;
	Eigen::Vector3d bc;

	// Cast a ray in the view direction starting from the mouse position
	double x = viewer.current_mouse_x;
	double y = viewer.core.viewport(3) - viewer.current_mouse_y;
	if (igl::unproject_onto_mesh(Eigen::Vector2f(x, y), viewer.core.view * viewer.core.model,
		viewer.core.proj, viewer.core.viewport, meshV, meshF, fid, bc))
	{
		viewer.data.set_colors(C);
		double d = 0;
		for (int i = 0; i < 3; i++)
		{
			//Skip border vertices
			if (cycles.row(meshF(fid, i)).squaredNorm() == 0)
				continue;
			double cur = bc(i);
			//Save closest vertex
			if (cur > d)
			{
				d = cur;
				//Left button for the red cycle (ring1)
				if (key == 0)
					ring1 = meshF(fid, i);
				//Right button for the green cycle (ring2)
				else
					ring2 = meshF(fid, i);
			}
		}
		draw_rings();
		return true;
	}
	return false;
};

int main()
{
	viewer.callback_key_down = &key_down;
	viewer.callback_mouse_down = &mouse_down;
	igl::readOBJ("../../data/chipped-torus.obj", meshV, meshF);

	igl::edge_topology(meshV, meshF, EV, FE, EF);
	igl::per_face_normals(meshV, meshF, norm);

	positiveIndices << .25, 0, 0,
					   .5,  0, 0,
					   .75, 0, 0,
					   1,   0, 0;

	negativeIndices << 0, .25, 0,
					   0, .5,  0,
					   0, .75, 0,
					   0, 1,   0;



	directional::dual_cycles(meshV, meshF, EV, EF, cycles);
	indices = Eigen::VectorXd::Zero(cycles.rows());
	genusLoop = cycles.rows() - 1;
	drawField();
	viewer.launch();
}
