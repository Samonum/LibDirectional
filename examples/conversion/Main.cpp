#include <iostream>
#include <directional/representative_to_raw.h>
#include <directional/adjustment_to_representative.h>
#include <directional/point_spheres.h>
#include <directional/drawable_field.h>
#include <directional/dual_cycles.h>
#include <directional/trivial_connection.h>
#include <directional/adjustment_to_raw.h>
#include <directional/representative_to_adjustment.h>
#include <Eigen/Core>
#include <igl/viewer/Viewer.h>
#include <igl/read_triangle_mesh.h>
#include <igl/per_face_normals.h>
#include <igl/edge_topology.h>
#include "Main.h"


Eigen::MatrixXi F, fieldF, meshF, EV, FE, EF;
Eigen::MatrixXd V, C, meshV, meshC, fieldV, fieldC, norm, representative, raw;
Eigen::VectorXd adjustmentField, other;
Eigen::VectorXi match;
Eigen::SparseMatrix<double, Eigen::RowMajor> cycles;
bool mode = true;

int N = 1;


void ConcatMeshes(const Eigen::MatrixXd &VA, const Eigen::MatrixXi &FA, const Eigen::MatrixXd &VB, const Eigen::MatrixXi &FB, Eigen::MatrixXd &V, Eigen::MatrixXi &F)
{
		V.resize(VA.rows() + VB.rows(), VA.cols());
		V << VA, VB;
		F.resize(FA.rows() + FB.rows(), FA.cols());
		F << FA, (FB.array() + VA.rows());
}

void UpdateViewer(igl::viewer::Viewer &viewer, Eigen::VectorXd &field)
{
	directional::adjustment_to_raw(meshV, meshF, EV, EF, norm, field, N, 0, raw);
	directional::drawable_field(meshV, meshF, raw, Eigen::RowVector3d(0, 0, 1), N, false, fieldV, fieldF, fieldC);

	ConcatMeshes(meshV, meshF, fieldV, fieldF, V, F);
	C.resize(F.rows(), 3);
	C << meshC, fieldC;

	viewer.data.set_mesh(V, F);
	viewer.data.set_colors(C);
	viewer.data.set_face_based(true);

}

bool key_down(igl::viewer::Viewer& viewer, unsigned char key, int modifiers)
{
	switch (key)
	{
	case 'Q':
		if (mode)
			UpdateViewer(viewer, other);
		else
			UpdateViewer(viewer, adjustmentField);
		break;
	}
	return true;
}

int main()
{
	igl::viewer::Viewer viewer;
	viewer.callback_key_down = &key_down;
	igl::readOBJ("../../data/spherers.obj", meshV, meshF);
	meshC = Eigen::RowVector3d(.8, .8, .8).replicate(meshF.rows(),1);

	igl::edge_topology(meshV, meshF, EV, FE, EF);
	igl::per_face_normals(meshV, meshF, norm);

	directional::dual_cycles(meshF, EV, EF, cycles);
	Eigen::VectorXi indices = Eigen::VectorXi::Zero(cycles.rows());

	indices(0) = N;
	indices(40) = N;
	
	directional::trivial_connection(meshV, meshF, EV, EF, cycles, indices, N, adjustmentField);
	
	UpdateViewer(viewer, adjustmentField);

	directional::adjustment_to_representative(meshV, meshF, EV, EF, adjustmentField, N, 0, representative);
	double angle;
	directional::representative_to_adjustment(meshV, meshF, EV, EF, FE, representative, N, other, match, angle);


	
	viewer.launch();
}