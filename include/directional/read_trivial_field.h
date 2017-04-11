#ifndef READ_TRIVIAL_FIELD
#define READ_TRIVIAL_FIELD
#include <cmath>
#include <igl/read_triangle_mesh.h>
#include <igl/writeOFF.h>
#include <igl/boundary_loop.h>
#include <directional/dual_cycles.h>
#include <Eigen/Core>
#include <string>
#include <sys/types.h>
#include <sys/stat.h>
#include <fstream>

#include <filesystem>
#include <direct.h>  

namespace directional
{
	// Writes a list of singularities to a file.
	// Inputs:
	//   fileName: The file name to which the singularities should be saved.
	//   singularities: The vector containing the singularities
	//   N: The degree of the field
	//   globalRotation: The angle of rotation between the vector on the first face and its basis in radians
	//   round: Whether or not singularities should be rounded before saving. 
	//          In most cases singularities should add up to whole integers.
	// Return:
	//   Whether or not the file was written successfully
	bool IGL_INLINE write_singularities(const std::string &fileName, const Eigen::VectorXd &singularities, const int N, const double globalRotation, const bool round = true)
	{
		try
		{
			std::ofstream f(fileName, std::ios::trunc);
			f << N << " " << singularities.size() << " " << globalRotation << std::endl;
			for (int i = 0; i < singularities.rows(); i++)
			{
				double s = round ? std::round(singularities(i)) : singularities(i);
				if (s)
					f << i << " " << s << std::endl;
			}
			f.close();
			return !f.fail();
		}
		catch(std::exception e)
		{
			return false;
		}
	}

	// Reads a list of singularities from a file
	// Inputs:
	//   fileName: The to be loaded file.
	// Outputs:
	//   singularities: The vector containing the singularities
	//   N: The degree of the field
	//   globalRotation: The angle of rotation between the vector on the first face and its basis in radians
	// Return:
	//   Whether or not the file was written successfully
	bool IGL_INLINE read_singularities(const std::string &fileName, Eigen::VectorXd &singularities, int& N, double& globalRotation)
	{
		try
		{
			std::ifstream f(fileName);
			int s;
			f >> N;
			f >> s;
			f >> globalRotation;

			singularities = Eigen::VectorXd::Zero(s);

			while(f)
			{
				int i;
				f >> i;
				f >> singularities.coeffRef(i);
			}
			f.close();
			return f.fail();
		}
		catch (std::exception e)
		{
			return false;
		}
	}

	// Writes a list of singularities and the mesh to files. The method will create a folder containing a "mesh.off" file with the mesh and "singularities.sing" file with the singularities.
	// Inputs:
	//   folder: The folder in which the singularities and mesh should be saved.
	//   V: List of Vertices
	//   F: List of faces
	//   singularities: The vector containing the singularities
	//   N: The degree of the field
	//   globalRotation: The angle of rotation between the vector on the first face and its basis in radians
	//   round: Whether or not singularities should be rounded before saving. 
	//          In most cases singularities should add up to whole integers.
	// Return:
	//   Whether or not the file was written successfully
	bool IGL_INLINE write_trivial_field(const std::string &folder, const Eigen::MatrixXd& V, const Eigen::MatrixXi& F, const Eigen::VectorXd &singularities, const int N, const double globalRotation, const bool round = true)
	{
		std::tr2::sys::path p(folder);

#if defined(_WIN32)
			_mkdir(p.generic_string().c_str()); // can be used on Windows
#else 
			mkdir(p.generic_string().c_str(), 0733); // can be used on non-Windows
#endif
		return igl::writeOFF((p / "mesh.off").generic_string(), V, F) &&
			directional::write_singularities((p / "singularities.sing").generic_string(), singularities, N, globalRotation, round);
	}

	// Will search for a mesh file "mesh.off" and a singularity file "singularities.sing" in the given folder and loads their data.
	// Inputs:
	//   folder: The folder which the singularities and mesh can be found.
	// Return:
	//   V: List of Vertices
	//   F: List of faces
	//   singularities: The vector containing the singularities
	//   N: The degree of the field
	//   globalRotation: The angle of rotation between the vector on the first face and its basis in radians
	//   Whether or not the file was written successfully
	bool IGL_INLINE read_trivial_field(const std::string &folder, Eigen::MatrixXd& V, Eigen::MatrixXi& F, Eigen::VectorXd &singularities, int N, double globalRotation)
	{
		std::tr2::sys::path p(folder);
		return igl::readOFF((p / "mesh.off").generic_string(), V, F) &&
			directional::read_singularities((p / "singularities.sing").generic_string(), singularities, N, globalRotation);
	}


}

#endif