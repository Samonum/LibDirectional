# libdirectional

### representative_to_raw

#### IGL_INLINE void representative_to_raw(const Eigen::MatrixXd& norm, const Eigen::MatrixXd& representative, Eigen::MatrixXd& raw, const int N)
Computes the raw vector field given a set of representative vectors.
Inputs:
norm: #F by 3 coordinates of the normals of each face.
representative: #F by 3 coordinates of representative vectors.
N: the degree of the field.
Outputs:
raw: #F by 3*N matrix with all N explicit vectors of each directional.
	
#### IGL_INLINE void representative_to_raw(const Eigen::MatrixXd& V, const Eigen::MatrixXi& F, const Eigen::MatrixXd& representative, Eigen::MatrixXd& raw, const int N)
Computes the raw vector field given a set of representative vectors.
This version recalculates the face normals every time it's called.
Inputs:
V: #V X 3 vertex coordinates.
F: #F by 3 face vertex indices.
representative: #F by 3 coordinates of representative vectors.
N: the degree of the field.
Outputs:
raw: #F by 3*N matrix with all N explicit vectors of each directional.