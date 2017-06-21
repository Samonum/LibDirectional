# libdirectional
Libdirectional is a C++ library for creating, manipulating and drawing vectorfields on 3D surfaces build on [libigl](www.github.com/libigl/libigl)<sup>[1](#fn1)</sup>. Libdirectional allows for n-rosy fields up to an arbitrary degree. Currently three types of vectorfields are supported, the trivial field, which attempts to create an as smooth as possible field given a set of singularities, the complex or globally optimal field, which tries to create a field that is as parallel as possible to a set of given example directionals and the polyvector field which generalizes the complex field to work for indepemdent 1<sup>n</sup> vector fields.

Parts of the code are based on the 2016 SIGGRAPH course on [directional field design](https://github.com/avaxman/DirectionalFieldSynthesis)<sup>[2](#fn2)</sup>. Furthermore some code is borrowed from the [libhedra](https://github.com/avaxman/libhedra) library<sup>[3](#fn3)</sup>.


## Installation
Libdirectional is a header only library where each file generally includes one function. To use the library simply add the _include_ directory to your include path and make sure libigl and its prerequisites are set up properly. After that you can include any files you need normally, using for example `#include <directional/trivial_field.h>`.

To use the examples simply clone the repository using:
```git
git clone --recursive https://github.com/Samonum/libdirectional.git
```

Then open a shell in the folder containing the example you wish to run and call:
```shell
mkdir build
cd build
cmake ..
```

This should properly set up the example including all dependencies.


## Representations

Libdirectional uses several different representations to describe vector fields.

| Method            | Representation                                                                                                            |
|-------------------|---------------------------------------------------------------------------------------------------------------------------|
| **Raw**               | \|F\| by 3\*N double matrix, representing each of the N vectors representing the full directional in the form X<sub>1</sub>, Y<sub>1</sub>, Z<sub>1</sub>, X<sub>2</sub>, Y<sub>2</sub>, Z<sub>2</sub> ... X<sub>N</sub>, Y<sub>N</sub>, Z<sub>N</sub>. Vectors are ordered in counter clockwise order.|
| **Representative**    | \|F\| by 3 double matrix representing the first vector in a directional. Only available for N-rosies.|
| **Adjustment** Angles | \|E\| by 1 double matrix representing the rotation between vectors on two neighbouring triangles, used in combination with a global rotation to uniquely define the field. Only available for N-rosies and does not encode vector length.|
| **complex**           | \|F\| by 1 complex double matrix, used by the complex field method. Allows for different length N-rosies, as long as each vector within a rosy is of the same length.              |
| **poly**vector        | \|F\| by N complex double matrix. Polyvectors are represented as a complex polynomial which uniquely devines each vector individually. Because each vector in a polyvector is defined individually they can not be converted to representative or Adjustment angle form.|

Libdirectional provides a number of conversion functions to switch between different forms of representation. Each of the functions is of the form \<method 1>\_to\_\<method 2>, where \<method 1> and \<method 2> are the bold parts of the method name in the above table in lower cacse. e.g. `adjustment_to_raw()`

For N-rosies you will most likely work primarily with the representative and adjustment angle representation, using the more verbose raw representation mostly for drawing. As polyvectors can not be described using a single representative vector it is required more often to use their raw representation.


## Trivial Field

The Trivial Field method attempts to create an as smooth as possible field covering the mesh. One of the main advantages of the trivial field is that generally all vectors in the field are of equal length.

In order to create the field it is generally needed to assign certain points around which the field will circle or from where the field will diverge. Even when this is not necessary it might still be desired to do so in order to control the flow of the field. These singularities form the primary input of the algorithm needed to calculate the field along with the basis cycle matrix describing the possible locations for singularities.<sup>[4](#fn4)</sup>

### Basis Cycles

The basis cycles form the cycles around which singularities are described on the mesh. The basis cycles are described using a \|cycles\| by \|E\| sparse matrix. Each row in the matrix describes one cycle and contains a 1 or -1 depending on the orrientation of the edge for each edge contained in the current cycle and a 0 for all other values. Cycles come in 3 types which should be described in order within the basis cycle matrix. The first \|V\| cycles are the contractible cycles around each vertex in ther mesh. Following these cycles are the border cycles surrounding the edges of the mesh. Finally are the generator cycles which describe cycles around handles.

The LibDirectional `dual_cycles()` method can calculate the proper basis cycles for a given mesh.

### Singularities

The singularities are described as an `Eigen::VectorXd` containing the singularity index of each basis cycle represents a 1/N rotation when following a vector around a particular cycle. In order to create a smooth field it is required that the indices of all singularities add up to N * the euler characteristic of the mesh and the indices are integers.

The *Singularities* example shows how one can calculate singularities given a vector field using the matching between neighbouring faces. Due to sampling issues the found singularities differ from the given singularities in most cases.

### Pre-calculations: Cycle Holonomy and Solver

To speed up field calculations it is possible to pre-calculate the cycle holonomy using the `cycle_holonomy()` function. 

To calculate the field LibDirectional uses the `Eigen::SimplicialLDLT` solver. This solver calculates the field in 2 steps, the first of which being solely dependant on the structure of the mesh, therefore it is possible to reuse this step as long as the mesh stays the same. To do so simply create a solver object and pass it into the `trivial_connection()` function.

### Examples
The *trivial_connection* example contains a small example program which allows picking basis cycles and altering their singularity index to see how this affects the field. It is possible to save fields generated with the *trivial_connection* example and load them in the *singularities* example.

Do all pre-calculations and generate a 4-rosy field given the vertices V and faces F of a mesh. Calculation of cycle holonomy and creation of the solver are optional, but will speed up subsequential field calculations for the same mesh.

Without precalculations:
```cpp
// Calculate basis cycles
Eigen::SparseMatrix<double, Eigen::RowMajor> cycles;
directional::dual_cycles(F, cycles);

// Set indices, should add up to N * the euler characteristic
Eigen::VectorXd indices = Eigen::VectorXd::Zero(cycles.rows());
indices[10] = N;
indices[20] = N;

// Calculate the field
Eigen::MatrixXd adjustAngles; // Field in the form of Adjustment Angles
double error;
directional::trivial_connection(V, F, cycles, indices, N, adjustAngles, error);

```

With precalculations:
```cpp
//Degree of the field (number of vectors within each directional)
int N = 4;

Eigen::MatrixXi EV, FE, EF;
Eigen::MatrixXd B1, B2, B3;

igl::edge_topology(V, F, EV, FE, EF);
igl::local_basis(V, F, B1, B2, B3);

//Calculate the basis cycles
Eigen::SparseMatrix<double, Eigen::RowMajor> cycles;
directional::dual_cycles(F, EV, EF, cycles);

//Calculate cycle holonomy
Eigen::MatrixXd cycleHolonomy;
cycle_holonomy(V, F, EV, EF, B1, B2, cycles, cycleHolonomy);

//Initialise solver
Eigen::SimplicialLDLT<Eigen::SparseMatrix<double> > solver;

//Set indices, should add up to N * the euler characteristic
Eigen::VectorXd indices = Eigen::VectorXd::Zero(cycles.rows());
indices[10] = N;
indices[20] = N;

//Calculate field
Eigen::MatrixXd adjustAngles; // Field in the form of Adjustment Angles
double error;
trivial_connection(V, F, basisCycles, indices, cycleHolonomy, solver, N, adjustAngles, error);
```

## Obtaining Singularities
Libdirectional is able to calculate singularities for a given field using the `singularities()` method. Singularities can be calculated using either the principal matching (obtained by passing a representative or raw field into the `principle_matching()` function) or the adjustment angles. Singularity calculation suffers from sampling issues, so unless calculated using the original adjustment angles you will most likely not obtain the same singularities as used to create the field. 

An illustration of these issues can be found in the *Singularities* example, which allows you to toggle between the original singularities and the calculated singularities. It is possible to save fields generated with the *trivial_connection* example and load them in the *singularities* example.

Singularity calculations work the same way for fields generated from singularities using the Trivial Field method and for other fields that have not been created using singularities.

## Complex Fields
Also known as Globally Optimal or As Parallel As Possible tries to generate a field that is as smooth as possible given a set of user defined soft constraints in the form of example vectors. Unlike the Trivial Field vector sizes for the Complex Field do not need to be equal.

### Defining Constraints
Constraints are defined as a pair of matrices of equal hight, refered to as `soft_ids` and `soft_values`. `Soft_ids` is a 1 wide integer matrix containing the ids of the faces (index in the F matrix) on which the constraints are placed. Meanwhile `soft_values` contains the x, y and z values of the matching vector, representing it in the same way as the representative vectors. Constraints do not need to be of unit length, but their size does matter for how they affect the field.

### Precomputing Solver
It is possible to speed up computations by precomputing the solver used to compute the complex field. This solver can then be reused to compute changes in the field as long as the mesh and **the ids of the constrained faces** remain the same.

### Examples
The *complex_field* example contains a small program which allows setting the soft constraints dynamically and see how it affects the field.

The below code creates a field of degree 3 and sets the first face so that its first vector aligns with the first edge. The other vectors are equally spaced to create a 3-rosy. V and F are the Vertices and Faces of the mesh. To see an example that alligns all vectors on the first face with an edge see the polyvector field example.

Without precalculations:
```cpp
// Set constraints
Eigen::VectorXi soft_ids(1);
Eigen::MatrixXd soft_values(1, N*3);
// Set to all faces that should be constrained
soft_ids(0) = 0;
// Set each matching row to the N vectors on that face
soft_values << V.row(F(0,0)) - V.row(F(0,1));

// Matrix containing the field
Eigen::MatrixXcd complex;

//Calculate the field
directional::complex_field(V, F, soft_ids, soft_values, N, complex);
```

With precalculations:
```cpp
//Degree of the field (number of vectors within each directional)
int N = 3;

Eigen::MatrixXi TT;
igl::triangle_triangle_adjacency(F, TT);
Eigen::MatrixXd B1, B2, x;
igl::local_basis(V, F, B1, B2, x);
    
// Set constraints
Eigen::VectorXi soft_ids(1);
Eigen::MatrixXd soft_values(1, N*3);
// Set to all faces that should be constrained
soft_ids(0) = 0;
// Set each matching row to the N vectors on that face
soft_values << V.row(F(0,0)) - V.row(F(0,1));
    
// Prepare the solver, must be recalculated whenever soft_ids changes (optional)
Eigen::SimplicialLDLT<Eigen::SparseMatrix<std::complex<double>>> solver;
Eigen::SparseMatrix<std::complex<double>> energy;
Complex_field_prepare_solver(V, F, TT, B1, B2, soft_id, N, solver, energy);

// Calculate the field
complex_field(B1, B2, soft_id, soft_value, solver, energy, N, complex);
```


## Polyvector Field
The Polyvector field is a generalisation of the standard complex field method, which allows defining each vector in a directional individually in both direction and length. Besides that it works largely the same way as the complex field.<sup>[5](#fn5)</sup> 

### Defining Constraints
Just like the Complex Field, the Polyvector Field takes a `soft_ids` matrix defining the face indices and a matching `soft_values` matrix defining the directionals on the faces, however the polyvector `soft_values` matrix is 3\*N wide, containing the X, Y, and Z values for each individual vector.

### Precomputing Solvers
It is possible to precompute the solvers for the Polyvector Field. To precompute the solvers one should pass an empty vector of SimplicialLDLT pointers (`std::vector<Eigen::SimplicialLDLT<Eigen::SparseMatrix<std::complex<double>>>*>`) into the `poly_vector_prepare_solvers()` function before passing them along with the energy matrix to the `poly_vector()` method.

The Solvers can be reused as long as the `soft_ids` remain the same, and must be properly `deleted` afterwards.

### Examples
The *poly_vectors* example shows the polyvector field in action, allowing the user to set constraints for each vector on each face individually. 

The below code creates a field of degree 3 and sets the first face so that its vectors each align with one edge of the triangle.  V and F are the Vertices and Faces of the mesh.

Without precalculations:
```cpp
// Set constraints
Eigen::VectorXi soft_ids(1);
Eigen::MatrixXd soft_values(1, N*3);
// Set to all faces that should be constrained
soft_ids(0) = 0;
// Set each matching row to the N vectors on that face
soft_values << V.row(F(0,0)) - V.row(F(0,1)), V.row(F(0,1)) - V.row(F(0,2)), V.row(F(0,2)) - V.row(F(0,0));

// Matrix containing the field
Eigen::MatrixXcd complex;

//Calculate the field
directional::poly_vector(V, F, soft_ids, soft_values, N, complex);
```

With precalculations:
```cpp
//Degree of the field (number of vectors within each directional)
int N = 3;

Eigen::MatrixXi TT;
igl::triangle_triangle_adjacency(F, TT);
Eigen::MatrixXd B1, B2, x;
igl::local_basis(V, F, B1, B2, x);
    
// Set constraints
Eigen::VectorXi soft_ids(1);
Eigen::MatrixXd soft_values(1, N*3);
// Set to all faces that should be constrained
soft_ids(0) = 0;
// Set each matching row to the N vectors on that face
soft_values << V.row(F(0,0)) - V.row(F(0,1)), V.row(F(0,1)) - V.row(F(0,2)), V.row(F(0,2)) - V.row(F(0,0));
    
// Prepare the solvers, must be recalculated whenever soft_ids changes (optional)
std::vector<Eigen::SimplicialLDLT<Eigen::SparseMatrix<std::complex<double>>>*> solvers;
std::vector<Eigen::SparseMatrix<std::complex<double>>> energy;
poly_vector_prepare_solvers(V, F, TT, B1, B2, soft_ids, N, solvers, energy);

// Matrix containing the field
Eigen::MatrixXcd poly;

// Calculate the field
poly_vector(B1, B2, soft_ids, soft_values, solvers, energy, N, poly);

...

// Make sure to properly dispose of all solvers
for (std::vector< Eigen::SimplicialLDLT<Eigen::SparseMatrix<std::complex<double>>>* >::iterator it = solvers.begin(); it != solvers.end(); ++it)
{
    delete (*it);
}
```

## Drawing fields
To draw a vector field first it must be converted into either its representative or raw representation. Then it can be passed to the `drawable_field()` method which will create a list of vertices, colors and faces to represent the field, finally those can be merged with the mesh data and passed to the viewer. If the viewer is set to take colors per vertex the `COLOR_PER_VERTEX` flag should be set.

There are several ways you can set the field color:<br>
The easiest way is to simply pass one color as a single row matrix. This will create one uniform ly colored field.<br>
The second option is to pass a matrix with \|F\| rows, each ro representing the color for the matching face.<br>
It is also possible to pass N\*\|F\| colors to give the color for every single vector on every single directional. In this case colors should be ordered to first give the color for every first vector on each face and only after that the color for every second vertex. e.g. F<sub>1</sub>1, F<sub>2</sub>1 ... F<sub>1</sub>2, F<sub>2</sub>2 ... F<sub>1</sub>N, F<sub>2</sub>N.<br>
Finally it is possible to pass one color per vertexafter setting the `PER_VECTOR_COLOR` flag.

By default the size of each vector is set to be related to the average edge length, as well as the length of the actual vector. The base length and with can be manually set if needed. If you want all vectors to be equal in size you can scale them by normalizing each vector in the field matrix. 

## References
<a name="fn1">1</a>: A. Jacobson and D. Panozzo and others, [libigl: A simple C++ geometry processing library](http://libigl.github.io/libigl/), 2016<br>
<a name="fn2">2</a>: A. Vaxman et al., [Directional Field Synthesis, Design, and Processing](https://github.com/avaxman/DirectionalFieldSynthesis), 2016<br>
<a name="fn3">3</a>: A. Vaxman et al., [libhedra](https://github.com/avaxman/libhedra), 2016<br>
<a name="fn4">4</a>: K. Crane and M. Desbrun and P. Schr&ouml;der, [Trivial Connections on Discrete Surfaces](https://www.cs.cmu.edu/~kmcrane/Projects/TrivialConnections/), 2010<br>
<a name="fn5">5</a>: O. Diamanti et al., [Designing N-PolyVector Fields with Complex Polynomials](http://igl.ethz.ch/projects/complex-roots/n-polyvector-fields.pdf), 2014
