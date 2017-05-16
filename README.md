# libdirectional



## Representations

Libdirectional uses several different representations to describe vector fields.

| Method            | Representation                                                                                                            |
|-------------------|---------------------------------------------------------------------------------------------------------------------------|
| **Raw**               | \|F\| by 3\*\|N\| double matrix, representing each of the N vectors representing the full directional in the form X<sub>1</sub>, Y<sub>1</sub>, Z<sub>1</sub>, X<sub>2</sub>, Y<sub>2</sub>, Z<sub>2</sub> ... X<sub>N</sub>, Y<sub>N</sub>, Z<sub>N</sub>. Vectors are ordered in counter clockwise order.|
| **Representative**    | \|F\| by 3 double matrix representing the first vector in a directional. Only available for N-rosies.|
| **Adjustment** Angles | \|E\| by 1 double matrix representing the rotation between vectors on two neighbouring triangles, used in combination with a global rotation to uniquely define the field. Only available for N-rosies and does not encode vector length.|
| **complex**           | \|F\| by 1 complex double matrix, used by the complex field method. Allows for different length N-rosies, as long as each vector within a rosy is of the same length.              |
| **poly**vector        | \|F\| by N complex double matrix. Polyvectors are represented as a complex polynomial which uniquely devines each vector individually. Because each vector in a polyvector is defined individually they can not be converted to representative or Adjustment angle form.|

Libdirectional provides a number of conversion functions to switch between different forms of representation. Each of the functions is of the form \<method 1>\_to\_\<method 2>, where \<method 1> and \<method 2> are the bold parts of the method name in the above table in lower cacse. e.g. `adjustment_to_raw()`

For N-rosies you will most likely work primarily with the representative and adjustment angle representation, using the more verbose raw representation mostly for drawing. As polyvectors can not be described using a single representative vector it is required more often to use their raw representation.


## Trivial Field

The Trivial Field method attempts to create an as smooth as possible field covering the mesh. One of the main advantages of the trivial field is that generally all vectors in the field are of equal length.

In order to create the field it is generally needed to assign certain points around which the field will circle or from where the field will diverge. Even when this is not necessary it might still be desired to do so in order to control the flow of the field. These singularities form the primary input of the algorithm needed to calculate the field along with the basis cycle matrix describing the possible locations for singularities.

### Basis Cycles

The basis cycles form the cycles around which singularities are described on the mesh. The basis cycles are described using a \|cycles\| by \|E\| sparse matrix. Each row in the matrix describes one cycle and contains a 1 or -1 depending on the orrientation of the edge for each edge contained in the current cycle and a 0 for all other values. Cycles come in 3 types which should be described in order within the basis cycle matrix. The first \|V\| cycles are the contractible cycles around each vertex in ther mesh. Following these cycles are the border cycles surrounding the edges of the mesh. Finally are the generator cycles which describe cycles around handles.

The LibDirectional `dual_cycles()` method can calculate the proper basis cycles for a given mesh.

### Singularities

The singularities are described as an `Eigen::VectorXd` containing the singularity index of each basis cycle represents a 1/N rotation when following a vector around a particular cycle. In order to create a smooth field it is required that the indices of all singularities add up to N * the euler characteristic of the mesh and the indices are integers.

The *Singularities* example shows how one can calculate singularities given a vector field using the matching between neighbouring faces. Due to sampling issues the found singularities differ from the given singularities in most cases.

### Pre-calculations: Cycle Holonomy and Solver

To speed up field calculations it is possible to pre-calculate the cycle holonomy using the `cycle_holonomy()` function. 

To calculate the field LibDirectional uses the Eigen::SimplicialLDLT solver. This solver calculates the field in 2 steps, the first of which being solely dependant on the structure of the mesh, therefore it is possible to reuse this step as long as the mesh stays the same. To do so simply create a solver object and pass it into the `trivial_connection()` function.

### Examples
The *trivial_connection* example contains a small example program which allows picking basis cycles and altering their singularity index to see how this affects the field. 

Do all pre-calculations and generate a 4-rosy field given the vertices V and faces F of a mesh. Calculation of cycle holonomy and creation of the solver are optional, but will speed up subsequential field calculations for the same mesh.

```cpp
//Degree of the field (number of vectors within each directional)
int N = 4;

Eigen::MatrixXi EV, FE, EF;
Eigen::MatrixXd B1, B2, B3;

igl::edge_topology(V, F, EV, FE, EF);
igl::local_basis(V, F, B1, B2, B3);

//Calculate the basis cycles
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
Libdirectional is able to calculate singularities for a given field using the `singularities()` method. Singularities can be calculated using either the principal matching (obtained by passing a representative or raw field into the `principle_matching()` function) or the adjustment angles. Singularity calculation suffers from sampling issues, so unless calculated using the original adjustment angles you will most likely not obtain the same singularities as used to create the field. An illustration of these issues can be found in the *Singularities* example, which allows you to toggle between the original singularities and the calculated singularities.

Singularity calculations work the same way for fields generated from singularities using the Trivial Field method and for other fields that have not been created using singularities.


## Complex Fields
Also known as Globally Optimal or As Parallel As Possible tries to generate a field that is as smooth as possible given a set of user defined soft constraints in the form of example vectors. Unlike the Trivial Field vector sizes for the Complex Field do not need to be of the ame size.

### Defining constraints
Constraints are defined as a pair of matrices of equal hight, refered to as `soft_ids` and `soft_values`. `Soft_ids` is a 1 wide integer matrix containing the ids of the faces (index in the F matrix) on which the constraints are placed. Meanwhile `soft_values` contains the x, y and z values of the matching vector, representing it in the same way as the representative vectors. Constraints do not need to be of unit length, but their size does matterfor how they affect the field.


