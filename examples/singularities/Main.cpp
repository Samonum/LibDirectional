#include <iostream>
#include <directional/drawable_field.h>
#include <directional/draw_singularities.h>
#include <directional/dual_cycles.h>
#include <directional/trivial_connection.h>
#include <directional/adjustment_to_representative.h>
#include <directional/representative_to_adjustment.h>
#include <directional/singularities.h>
#include <directional/principal_matching.h>
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
Eigen::VectorXd indices, calculatedIndices;
Eigen::SparseMatrix<double, Eigen::RowMajor> cycles;
igl::viewer::Viewer viewer;

Eigen::MatrixXd positiveIndices(4, 3),
				negativeIndices(4,3);

int sing_mode = 0;

int N = 4;

void ConcatMeshes(const Eigen::MatrixXd &VA, const Eigen::MatrixXi &FA, const Eigen::MatrixXd &VB, const Eigen::MatrixXi &FB, Eigen::MatrixXd &V, Eigen::MatrixXi &F)
{
	V.resize(VA.rows() + VB.rows(), VA.cols());
	V << VA, VB;
	F.resize(FA.rows() + FB.rows(), FA.cols());
	F << FA, (FB.array() + VA.rows());
}

void draw_singularities()
{
	Eigen::MatrixXd spheres;
	if(sing_mode)
		directional::draw_singularities(meshV, calculatedIndices, positiveIndices, negativeIndices, .015, singV, singF, singC);
	else
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

void drawField()
{
	double e;
	directional::trivial_connection(meshV, meshF, EV, EF, cycles, indices, N, adjustmentField, e);
	std::cout << "error: " << e << std::endl;
	directional::adjustment_to_representative(meshV, meshF, EV, EF, adjustmentField, N, 0, representative);
	double r;
	directional::representative_to_adjustment(meshV, meshF, EV, FE, EF, representative, N, other, r);

	Eigen::VectorXi matching;
	directional::principal_matching(meshV, meshF, EV, EF, FE, representative, N, matching);
	directional::singularities(cycles, matching, calculatedIndices);

	Eigen::MatrixXd color;
	directional::drawable_field(meshV, meshF, representative, Eigen::RowVector3d(0, 0, 1), N, false, fieldV, fieldF, fieldC);

	meshC = Eigen::RowVector3d(1, 1, 1).replicate(meshF.rows(), 1);
	draw_singularities();
}

bool key_down(igl::viewer::Viewer& viewer, int key, int modifiers)
{
	switch (key)
	{
	case '1':
		sing_mode = 0;
		draw_singularities();
		break;
	case '2':
		sing_mode = 1;
		draw_singularities();
		break;
	}
	return true;
}


int main()
{
	viewer.callback_key_down = &key_down;
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
	indices[20] = N;
	indices[120] = -2*N;
	indices[220] = N;
	drawField();
	viewer.launch();
}
