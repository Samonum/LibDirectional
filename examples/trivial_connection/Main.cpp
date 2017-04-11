#include <iostream>
#include <directional/drawable_field.h>
#include <directional/draw_singularities.h>
#include <directional/dual_cycles.h>
#include <directional/trivial_connection.h>
#include <directional/adjustment_to_representative.h>
#include <directional/draw_cycles.h>
#include <directional/read_trivial_field.h>
#include <Eigen/Core>
#include <igl/viewer/Viewer.h>
#include <igl/read_triangle_mesh.h>
#include <igl/per_face_normals.h>
#include <igl/edge_topology.h>
#include <igl/unproject_onto_mesh.h>
#include <igl/boundary_loop.h>
#include "Main.h"


Eigen::MatrixXi F, fieldF, meshF, singF, EV, FE, EF;
Eigen::MatrixXd V, C, meshV, meshC, fieldV, fieldC, singV, singC, representative;
Eigen::VectorXd adjustmentField, other;
Eigen::VectorXi match; 
Eigen::VectorXd indices;
Eigen::SparseMatrix<double, Eigen::RowMajor> cycles;
igl::viewer::Viewer viewer;

Eigen::MatrixXd positiveIndices(4, 3),
				negativeIndices(4,3);

int euler;
int generators;

int N = 4;
int ring1 = 200;

bool drag = false;

void ConcatMeshes(const Eigen::MatrixXd &VA, const Eigen::MatrixXi &FA, const Eigen::MatrixXd &VB, const Eigen::MatrixXi &FB, Eigen::MatrixXd &V, Eigen::MatrixXi &F)
{
	V.resize(VA.rows() + VB.rows(), VA.cols());
	V << VA, VB;
	F.resize(FA.rows() + FB.rows(), FA.cols());
	F << FA, (FB.array() + VA.rows());
}

void calculate_field()
{
	double e;
	directional::trivial_connection(meshV, meshF, EV, EF, cycles, indices, N, adjustmentField, e);
	directional::adjustment_to_representative(meshV, meshF, EV, EF, adjustmentField, N, 0, representative);

	std::cout << e << std::endl;

	int sum = round(indices.head(indices.size() - generators).sum());
	if (euler*N != sum)
	{
		std::cout << "Warning: All non-generator singularities should add up to N * the Euler characteristic."<<std::endl;
		std::cout << "Total indices: " << sum << std::endl;
		std::cout << "Expected: " << euler *N << std::endl;
	}
	Eigen::MatrixXd color;
	color.resize(representative.rows()*N, 3);
	color << Eigen::RowVector3d(0, 1, 0).replicate(representative.rows(), 1),
		Eigen::RowVector3d(0, 0, 1).replicate((N - 1) *representative.rows(), 1);
	directional::drawable_field(meshV, meshF, representative, color, N, false, fieldV, fieldF, fieldC);
}

void draw_field()
{
	meshC = Eigen::RowVector3d(1, 1, 1).replicate(meshF.rows(), 1);
	directional::draw_cycles(EF, cycles, Eigen::Vector3d(1, 0, 0), ring1, meshC);

	Eigen::MatrixXd spheres;
	directional::draw_singularities(meshV, indices, positiveIndices, negativeIndices, .015, singV, singF, singC);

	Eigen::MatrixXd a;
	Eigen::MatrixXi b;
	ConcatMeshes(meshV, meshF, fieldV, fieldF, a, b);
	if (singV.rows())
	{
		ConcatMeshes(a, b, singV, singF, V, F);
		C.resize(F.rows(), 3);
		C << meshC, fieldC, singC;
	}
	else
	{
		V = a;
		F = b;
		C.resize(F.rows(), 3);
		C << meshC, fieldC;
	}

	viewer.data.clear();
	viewer.data.set_face_based(true);
	viewer.data.set_mesh(V, F);
	viewer.data.set_colors(C);
}

void update_mesh()
{

	igl::edge_topology(meshV, meshF, EV, FE, EF);

	directional::dual_cycles(meshF, EV, EF, cycles);
	std::vector<std::vector<int>> boundaryLoops;
	igl::boundary_loop(meshF, boundaryLoops);
	euler = meshV.rows() - EV.rows() + meshF.rows();
	generators = cycles.rows() - meshV.rows() - boundaryLoops.size();

}

void paint_ring()
{
	//Clear out the mesh
	meshC = Eigen::RowVector3d(1, 1, 1).replicate(meshF.rows(), 1);
	//Draw the new cycles
	directional::draw_cycles(EF, cycles, Eigen::Vector3d(1, 0, 0), ring1, meshC);
	if (singC.rows())
		C << meshC, fieldC, singC;
	else
		C << meshC, fieldC;
	viewer.data.set_colors(C);
}

bool key_down(igl::viewer::Viewer& viewer, int key, int modifiers)
{
	int borders;
	switch (key)
	{
	case '-':
	case '_':
		indices(ring1)--;
		draw_field();
		break;
	case '+':
	case '=':
		indices(ring1)++;
		draw_field();
		break;
	case 'C':
		calculate_field();
		draw_field();
		break;
	case 'D':
		drag = !drag;
		break;
	case 'B':
		borders = cycles.rows()- (generators + meshV.rows());
		if (borders)
		{
			//Loop through the border cycles.
			if (ring1 >= meshV.rows() && ring1 < meshV.rows() + borders - 1)
				ring1++;
			else
				ring1 = meshV.rows();
			paint_ring();
		}
		break;
	case 'G':
		if (generators)
		{
			//Loop through the generators cycles.
			if (ring1 >= cycles.rows() - generators && ring1 < cycles.rows() - 1)
				ring1++;
			else
				ring1 = cycles.rows() - generators;
			paint_ring();
		}
		break;
	case 'S':
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
				if (key == 0)
					ring1 = meshF(fid, i);
			}
		}
		paint_ring();
		return true;
	}
	return false;
};

int main()
{
	viewer.callback_key_down = &key_down;
	viewer.callback_mouse_down = &mouse_down;
	igl::readOBJ("../../data/half-torus.obj", meshV, meshF);

	positiveIndices << .25, 0, 0,
					   .5,  0, 0,
					   .75, 0, 0,
					   1,   0, 0;

	negativeIndices << 0, .25, 0,
					   0, .5,  0,
					   0, .75, 0,
					   0, 1,   0;

	update_mesh();

	indices = Eigen::VectorXd::Zero(cycles.rows());

	calculate_field();
	draw_field();
	viewer.launch();
}
