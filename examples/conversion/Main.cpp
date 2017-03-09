#include <iostream>
#include <directional/representative_to_raw.h>
#include <directional/point_spheres.h>
#include <directional/drawable_field.h>
#include <Eigen/Core>
#include <igl/viewer/Viewer.h>
#include <igl/read_triangle_mesh.h>
#include <igl/local_basis.h>

Eigen::MatrixXi F, fieldF, meshF;
Eigen::MatrixXd V, C, meshV, meshC, fieldV, fieldC, norm(5,3), representative(5,3), raw, B1, B2, B3;

int N = 4;


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

	igl::local_basis(meshV, meshF, B1, B2, B3);
	directional::drawable_field(meshV, meshF, B1, Eigen::RowVector3d(0, 0, 1), N, false,  fieldV, fieldF, fieldC);
	ConcatMeshes(meshV, meshF, fieldV, fieldF, V, F);
	C.resize(F.rows(), 3);
	C << meshC, fieldC;

	viewer.data.set_mesh(V, F);
	viewer.data.set_colors(C);
	viewer.data.set_face_based(true);
	viewer.launch();
}