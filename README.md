# libdirectional



## Representations

Libdirectional uses several different representations to describe vector fields.

| Method            | Representation                                                                                                            |
|-------------------|---------------------------------------------------------------------------------------------------------------------------|
| **Raw**               | \|F\| by 3\*\|N\| double matrix, representing each of the N vectors representing the full directional in the form X<sub>1</sub>, Y<sub>1</sub>, Z<sub>1</sub>, X<sub>2</sub>, Y<sub>2</sub>, Z<sub>2</sub> ... X<sub>N</sub>, Y<sub>N</sub>, Z<sub>N</sub>,  |
| **Representative**    | \|F\| by 3 double matrix representing the first vector in a directional, only available for N-rosies                                     |
| **Adjustment** Angles | \|E\| by 1 double matrix representing the rotation between vectors on two neighbouring triangles, used in combination with a global rotation to uniquely define the field, only available for N-rosies               |

Libdirectional provides a number of conversion functions to switch between different forms of representation. Each of the functions is of the form \<method 1>\_to\_\<method 2>, where \<method 1> and \<method 2> are the bold parts of the method name in the above table in lower cacse. e.g. `adjustment\_to\_raw()`

For N-rosies you will most likely work primarily with the representative and adjustment angle representation, using the raw representation mostly for drawing.


## Trivial Field

The Trivial Field method attempts to create an as smooth as possible field covering the mesh. One of the main advantages of the trivial field is that generally all vectors in the field are of equal length.

In order to create the field it is generally needed to assign certain points around which the field will circle or from where the field will diverge. Even when this is not necessary it might still be desired to do so in order to control the flow of the field. These singularities form the primary input of the algorithm needed to calculate the field along with the basis cycle matrix describing the possible locations for singularities.

### Basis Cycles

The basis cycles form the cycles around which singularities are described on the mesh. The basis cycles are described using a \|cycles\| by \|E\| sparse matrix. Each row in the matrix describes one cycle and contains a 1 or -1 depending on the orrientation of the edge for each edge contained in the current cycle and a 0 for all other values. Cycles come in 3 types which should be described in order within the basis cycle matrix. The first \|V\| cycles are the contractible cycles around each vertex in ther mesh. Following these cycles are the border cycles surrounding the edges of the mesh. Finally are the generator cycles which describe cycles around handles.

The LibDirectional `dual_cycles()` method can calculate the proper basis cycles for a given mesh.

### Singularities

The singularities are described as a Eigen VectorXi containing the singularity index of each basis cycle represents a 1/N rotation when following a vector around a particular cycle. In order to create a smooth field it is required that the indices of all singularities add up to N * the euler characteristic of the mesh. 

### Pre-calculations: Cycle Holonomy and Solver

To speed up field calculations it is possible to pre-calculate the cycle holonomy using the `cycle_holonomy()` function. 

To calculate the field LibDirectional uses the Eigen::SimplicialLDLT solver. This solver calculates the field in 2 steps, the first of which being solely dependant on the structure of the mesh, therefore it is possible to reuse this step as long as the mesh stays the same. To do so simply create a solver object and pass it into the `trivial_connection()` function.


### Example
Do all pre-calculations and generate a 4-rosy field given the vertices V and faces F of a mesh. Calculation of cycle holonomy and creation of the solver are optional, but will speed up subsequential field calculations for the same mesh.

```cpp
int N = 4;

Eigen::MatrixXi EV, FE, EF;
Eigen::MatrixXd B1, B2, B3;

igl::edge_topology(V, F, EV, FE, EF);
igl::local_basis(V, F, B1, B2, B3);

//Calculate the basis cycles
directional::dual_cycles(F, EV, EF, cycles);

//Set indices, should add up to N * the euler characteristic
Eigen::VectorXi indices = Eigen::VectorXd::Zero(cycles.rows());
indices[10] = N;
indices[20] = N;

//Calculate cycle holonomy
Eigen::MatrixXd cycleHolonomy;
cycle_holonomy(V, F, EV, EF, B1, B2, cycles, cycleHolonomy);

//Initialise solver
Eigen::SimplicialLDLT<Eigen::SparseMatrix<double> > solver;

//Calculate field
Eigen::MatrixXd adjustAngles; // Field in the form of Adjustment Angles
double error;
trivial_connection(V, F, basisCycles, indices, cycleHolonomy, solver, N, adjustAngles, error);
```

### IGL\_INLINE void representative\_to\_raw(const Eigen::MatrixXd& norm, const Eigen::MatrixXd& representative, Eigen::MatrixXd& raw, const int N)
Computes the raw vector field given a set of representative vectors.
Inputs:
norm: #F by 3 coordinates of the normals of each face.
representative: #F by 3 coordinates of representative vectors.
N: the degree of the field.
Outputs:
raw: #F by 3*N matrix with all N explicit vectors of each directional.
	
#### IGL\_INLINE void representative\_to\_raw(const Eigen::MatrixXd& V, const Eigen::MatrixXi& F, const Eigen::MatrixXd& representative, Eigen::MatrixXd& raw, const int N)
Computes the raw vector field given a set of representative vectors.
This version recalculates the face normals every time it's called.
Inputs:
V: #V X 3 vertex coordinates.
F: #F by 3 face vertex indices.
representative: #F by 3 coordinates of representative vectors.
N: the degree of the field.
Outputs:
raw: #F by 3*N matrix with all N explicit vectors of each directional.
