# README for Computational Geophysics
## Introduction
This project is a collection of scripts that regroup the deliverable for the course Computational Geophysics at Utrecht University GEO4-1427, directed by [Asst. Prof. Cedric Thieulot](https://cedricthieulot.net/). 

Normally, in any practical situation, we would not rewrite software from scratch, but I believe there is no better way to learn how something works than to completely rebuild it,. Therefore, this project will aim to build all the basic elements of a geomechanical simulation software as well as compare different schemes, solvers and methods to achieve this goal. In the same vein, the focus while coding these scripts will be readability, so the code will not be optimized for speed as it is often the case for python scripts.

## Structure
- `bin/`: contains the executable files
  - `main.py`: answer of the exercises from [Fieldstone](https://cedricthieulot.net/manual.pdf). Parameters can be set in the beginning to fit different way to solve either diffusion or advection equations.
  - `numerical-integration-rules-ex#`: comparison of two methods of numerical approximation of a solution of integral equations: the midpoint rule or the trapezoidal rule. 
  - `numerical-integration-quadrature-ex#`: experimentation of the quadrature rule which is another way to solve integral equations in an approximate way. It consists in replacing the equation by a sum of integrands evaluated at several quadrature points multiplied by quadrature weights.
- `modules/`: contains the functions and classes
  - `methods.py`: highest level of the tree, contain scheme and methodology.
  - `solvers.py`: mid level of the tree, contain the solvers.
  - `operations.py`: lowest level of the tree, contain simple calculations or function used to present results.
- `output/`: contains the output files
- `test/`: contains the test files
- `environment.yml`: file to set environment with conda

All exercises for FDM and FEM are executable from the `bin/main.py` script, with instructions for each exercise present in [Fieldstone](https://cedricthieulot.net/manual.pdf). In addition, three scripts contain the comparison of different numerical integration methods, located in the `bin/` directory. The `modules/` directory contains all the functions that are not called directly and are organized into three different scripts: `methods.py`, `solvers.py` and `operations.py`. `methods. py` is the highest level of the tree, so it is able to call the other scripts in the same folder, and contains the complete schema and methodology. `solvers.py` is one level lower than `methods.py` and contains the different solvers, this script is able to call the other scripts `operations.py`. This last script is at the lowest level of the tree and will not call any other function in this project, it contains simple calculations or formatting functions such as result plotting.


## Step-by-step evolution
This section will explain step by step each feature of this project. Each step will contain 1D and 2D implementation as well as solving the diffusion and advection equation. 

[//]: # (Talk about all the elements &#40;solvers, scheme, etc.&#41;)

[//]: # (Starting with FDM &#40;Finit Difference Method&#41;)

### Finite Element Method
First, let's consider a 1D problem where we have a rod of 1000m long with constant temperature at the extremity
(Dirichlet Boundary Conditions). We want to calculate the temperature profile at stead-state achieved after
multiples iterations, and we will use the finite difference method to achieve that goal.
TODO: explain the finite difference method.
We will use two approach, the explicit and the implicit method.
Ex2: In the explicit method, we use the node at the current time step to calculate the node at the next time step,
meaning the time derivative is taken forward in time.

Ex3+: In the implicit method, we will have a derivation backward in time by using only one node at the current time
step. In this method, the solution is unconditionally stable but not necessarly accurate, thus there is no
restriction on the time step to obtain a solution. However, it will not necessarily converge, so it's important to
choose the time step wisely and verify whether the convergence is achieved.
Initially, the solver using the implicit did not converge (with the given matrix). This is because the solver is not
able to handle when matrix A is not diagonally dominant. Therefor, we have a function that will verify this and test
all rows permutations to try to find at least one that is diagonally dominant. If it is not possible, then an error
will be raised.

We have now a working script to solve the diffusion heat equation using the implicit method with three different
solvers: Jacobi, Gauss-Seidel and, for verification, SciPy's solver. As expected, SciPy's solver is the fastest as
it is the most optimised, followed by Jacobi and in the last position Gauss-Seidel. On paper, the Gauss-Seidel
should converge faster because there is no need for an initial guess but this is not the case in this implementation.
Additionally, the Crank-Nicolson scheme was tested, which is an hybrid between an implicit and explicit approach, thus
limiting the drawbacks of each methods and beating all other method in term of time while being second order
accurate in time and space.

Ex 6: The next step consist of solving the advection equation using two different methods: FTCS and Lax-Friderichs.
The FTCS method is a first order scheme that is based on the forward time derivative and the backward space with the
main problem being that this method is unconditionally unstable and will explode for any given time step, as we can
see with the results. To palliate this problem, the Lax-Friedrichs method is using the mean of the two adjacent T
nodes to calculate the next iteration instead of the current and, this time, output the correct results for any
dt <= h/u.

Ex 7+: We now will do the same in 2D. The principle remain the same but we will use slightly different formulas and
the matrix will be fill differently. We did this for diffusion explicit, diffusion implicit and advection.

Ex FEM-1 and FEM-2: Implementation of the finite element method for, respectively, the diffusion and the advection equations in 1D.


### Numerical Integration Methods


### Finite Element Method

