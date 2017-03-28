#include <iostream>
#include <directional/representative_to_raw.h>
#include <directional/adjustment_to_representative.h>
#include <directional/point_spheres.h>
#include <directional/drawable_field.h>
#include <directional/dual_cycles.h>
#include <directional/trivial_connection.h>
#include <directional/adjustment_to_raw.h>
#include <directional/representative_to_adjustment.h>
#include <directional/point_spheres.h>
#include <directional\principal_matching.h>
#include <directional/singularities.h>
#include <directional/draw_cycles.h>
#include <Eigen/Core>
#include <igl/viewer/Viewer.h>
#include <igl/read_triangle_mesh.h>
#include <igl/per_face_normals.h>
#include <igl/edge_topology.h>
#include "Main.h"


Eigen::MatrixXi F, fieldF, meshF, singF, EV, FE, EF;
Eigen::MatrixXd V, C, meshV, meshC, fieldV, fieldC, singV, singC, norm, representative, raw, B1, B2, B3;
Eigen::VectorXd adjustmentField, other;
Eigen::VectorXi match; 
Eigen::SparseVector<int> indices;
Eigen::SparseMatrix<double, Eigen::RowMajor> cycles;
bool mode = true;

int N = 3;
int ring;



// define the format you want, you only need one instance of this...
const static Eigen::IOFormat CSVFormat(Eigen::StreamPrecision, Eigen::DontAlignCols, ", ", "\n");
template<typename _Scalar, int _Rows, int _Cols, int _Options, int _MaxRows, int _MaxCols>
void writeToCSVfile(std::string name, Eigen::Matrix<_Scalar, _Rows, _Cols, _Options, _MaxRows, _MaxCols> matrix)
{
	std::ofstream file(name.c_str());
	file << matrix.format(CSVFormat);
}

void ConcatMeshes(const Eigen::MatrixXd &VA, const Eigen::MatrixXi &FA, const Eigen::MatrixXd &VB, const Eigen::MatrixXi &FB, Eigen::MatrixXd &V, Eigen::MatrixXi &F)
{
		V.resize(VA.rows() + VB.rows(), VA.cols());
		V << VA, VB;
		F.resize(FA.rows() + FB.rows(), FA.cols());
		F << FA, (FB.array() + VA.rows());
}

void UpdateViewer(igl::viewer::Viewer &viewer, Eigen::MatrixXd &raw)
{
	Eigen::MatrixXd color;
	color.resize(raw.rows()*N, 3);
	color << Eigen::RowVector3d(0, 1, 0).replicate(raw.rows(), 1),
		Eigen::RowVector3d(0, 0, 1).replicate((N - 1) *raw.rows(), 1);
	raw.rowwise().normalize();
	directional::drawable_field(meshV, meshF, raw, color, N, false, fieldV, fieldF, fieldC);
	
	meshC = Eigen::RowVector3d(.8, .8, .8).replicate(meshF.rows(),1);
	directional::draw_cycles(EF, cycles, Eigen::Vector3d(1, 0, 0), ring, meshC);

	Eigen::MatrixXd spheres;
	//spheres.resize(2, 3);
	//spheres << meshV.row(30), meshV.row(200);
	//directional::point_spheres(spheres, .008, Eigen::RowVector3d(1, 0, 0).replicate(2, 1), 10, false, singV, singF, singC);
	Eigen::MatrixXd a;
	Eigen::MatrixXi b;
	ConcatMeshes(meshV, meshF, fieldV, fieldF, a, b);
	ConcatMeshes(a, b, singV, singF, V, F);

	C.resize(F.rows(), 3);
	C << meshC, fieldC, singC;
	viewer.data.clear();
	viewer.data.set_mesh(V, F);
	viewer.data.set_colors(C);
	viewer.data.set_face_based(true);
}

void UpdateViewer(igl::viewer::Viewer &viewer, Eigen::VectorXd &field)
{
	directional::adjustment_to_representative(meshV, meshF, EV, EF, field, N, 0, raw);
	UpdateViewer(viewer, raw);
}

bool key_down(igl::viewer::Viewer& viewer, unsigned char key, int modifiers)
{
	switch (key)
	{
	case 'N':
		if (++ring >= cycles.rows())
			ring -= cycles.rows();
		UpdateViewer(viewer, adjustmentField);
		break;
	case 'P':
		if (--ring < 0)
			ring += cycles.rows();
		UpdateViewer(viewer, adjustmentField);
		break;
	default:;
	}
	return true;
}

int main()
{
	igl::viewer::Viewer viewer;
	viewer.callback_key_down = &key_down;
	igl::readOBJ("../../data/loop.obj", meshV, meshF);

	igl::edge_topology(meshV, meshF, EV, FE, EF);
	igl::per_face_normals(meshV, meshF, norm);

	directional::dual_cycles(meshV, meshF, EV, EF, cycles);
	//std::cout << cycles.row(meshV.rows()) << std::endl;
	//std::cout << cycles.row(meshV.rows() + 1) << std::endl;
	indices.resize(cycles.rows());

	ring = cycles.rows() - 1;

	writeToCSVfile("cycles.txt", Eigen::MatrixXd(cycles));
	indices.insert(cycles.rows() - 1) = N;
	indices.insert(cycles.rows() - 2) = -N;
	//indices.insert(30) = N;
	//indices.insert(200) = N;
	//indices.insert(300) = -2*N;
	
	directional::trivial_connection(meshV, meshF, EV, EF, cycles, indices, N, adjustmentField);
	
	directional::adjustment_to_representative(meshV, meshF, EV, EF, adjustmentField, N, 0, representative);
	Eigen::VectorXi matching;
	Eigen::VectorXd sing;
	directional::principal_matching(meshV, meshF, EV, EF, FE, representative, N, matching);
	directional::singularities(cycles, matching, sing);
	writeToCSVfile("sing.txt", sing);
	

	//double angle;
	//directional::representative_to_adjustment(meshv, meshf, ev, ef, fe, representative, n, other, angle);

	UpdateViewer(viewer, adjustmentField);
	viewer.launch();
}