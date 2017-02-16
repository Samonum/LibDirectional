#include <iostream>
#include <directional/representative_to_raw.h>
#include <Eigen/Core>
#include <Eigen/Geometry> 

Eigen::MatrixXi F;
Eigen::MatrixXd V, norm(5,3), representative(5,3), raw;

int N = 4;
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
}