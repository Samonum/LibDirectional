#include <iostream>
#include <directional/representative_to_raw.h>
#include <directional/adjustment_to_representative.h>
#include <directional/point_spheres.h>
#include <directional/drawable_field.h>
#include <directional/dual_cycles.h>
#include <directional/trivial_connection.h>
#include <Eigen/Core>
#include <igl/viewer/Viewer.h>
#include <igl/read_triangle_mesh.h>
#include <igl/local_basis.h>
#include <igl/edge_topology.h>


Eigen::MatrixXi F, fieldF, meshF, EV, FE, EF;
Eigen::MatrixXd V, C, meshV, meshC, fieldV, fieldC, norm(5,3), representative(5,3), raw, B1, B2, B3;
Eigen::VectorXd adjustmentField;
Eigen::SparseMatrix<double, Eigen::RowMajor> cycles;

int N = 1;


void UpdateCurrentView(igl::viewer::Viewer& viewer)
{

}

void ConcatMeshes(const Eigen::MatrixXd &VA, const Eigen::MatrixXi &FA, const Eigen::MatrixXd &VB, const Eigen::MatrixXi &FB, Eigen::MatrixXd &V, Eigen::MatrixXi &F)
{
		V.resize(VA.rows() + VB.rows(), VA.cols());
		V << VA, VB;
		F.resize(FA.rows() + FB.rows(), FA.cols());
		F << FA, (FB.array() + VA.rows());
}

int main()
{
	for(int i = 0; i < norm.rows(); i++)
		norm.row(i) << ((i % 3) == 0), ((i % 3) == 1), ((i % 3) == 2);
	representative.row(0) << 0, 1, 0;
	representative.row(1) << 1, 0, 0;
	representative.row(2) << 0, -1, 0;
	representative.row(3) << 0, 0, 1;
	representative.row(4) << 0, 0, -1;
	directional::representative_to_raw(norm, representative, raw, N);

	std::cout << "normals: " << std::endl;
	std::cout << norm << std::endl;
	std::cout << "representative: " << std::endl;
	std::cout << representative << std::endl;
	std::cout << "raw: " << std::endl;
	std::cout << raw << std::endl;

	std::cin.get();
	igl::viewer::Viewer viewer;
	//viewer.callback_key_down = &key_down;
	UpdateCurrentView(viewer);
	igl::readOBJ("../../data/spherers.obj", meshV, meshF);
	meshC = Eigen::RowVector3d(.8, .8, .8).replicate(meshF.rows(),1);

	//igl::local_basis(meshV, meshF, B1, B2, B3);
	igl::edge_topology(meshV, meshF, EV, FE, EF);


	directional::dual_cycles(meshF, EV, EF, cycles);
	Eigen::VectorXi indices = Eigen::VectorXi::Zero(cycles.rows());

	indices(0) = N;
	indices(20) = N;

	directional::trivial_connection(meshV, meshF, EV, EF, cycles, indices, N, adjustmentField);

	directional::adjustment_to_representative(meshV, meshF, EV, EF, adjustmentField, N, 0, B1);

	directional::drawable_field(meshV, meshF, B1, N, Eigen::RowVector3d(0, 0, 1), fieldV, fieldF, fieldC);
	ConcatMeshes(meshV, meshF, fieldV, fieldF, V, F);
	C.resize(F.rows(), 3);
	C << meshC, fieldC;

	

	viewer.data.set_mesh(V, F);
	viewer.data.set_colors(C);
	viewer.data.set_face_based(true);
	viewer.launch();
}