# libdirectional



## Representations

Libdirectional uses several different representations to describe vector fields.

| Method            | Representation                                                                                                            |
|-------------------|---------------------------------------------------------------------------------------------------------------------------|
| **Raw**               | \|F\| by 3*\|N\| double matrix, representing each of the N vectors representing the full directional in the form X<sub>1</sub>, Y<sub>1</sub>, Z<sub>1</sub>, X<sub>2</sub>, Y<sub>2</sub>, Z<sub>2</sub> ... X<sub>N</sub>, Y<sub>N</sub>, Z<sub>N</sub>,  |
| **Representative**    | \|F\| by 3 double matrix representing the first vector in a directional, only available for N-rosies                                     |
| **Adjustment** Angles | \|E\| by 1 double matrix representing the rotation between vectors on two neighbouring triangles, used in combination with a global rotation to uniquely define the field, only available for N-rosies               |

Libdirectional provides a number of conversion functions to switch between different forms of representation. Each of the functions is of the form \<method 1>\_to\_\<method 2>, where \<method 1> and \<method 2> are the bold parts of the method name in the above table in lower cacse. e.g. adjustment_to_raw()

For N-rosies you will most likely work primarily with the representative and adjustment angle representation, using the raw representation mostly for drawing.




### IGL_INLINE void representative_to_raw(const Eigen::MatrixXd& norm, const Eigen::MatrixXd& representative, Eigen::MatrixXd& raw, const int N)
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