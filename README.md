# README for Computational Geophysics

This is a project containing all the deliverable for the course Computational Geophysics at Utrecht University.

Intially, each exercise were contain in a separated python file in `bin/`. But for simplicity, all the answers and the thought process will be in the `bin/main.py` file.
When there is a function or classes that is used elsewhere, it will be in `modules/`.

## Structure
- `bin/`: contains the executable files
- `data/`: contains the data files
- `modules/`: contains the functions and classes
- `output/`: contains the output files
- `test/`: contains the test files
- `environment.yml`: file to set environment with conda
- `requirements.txt`: file to install dependencies with pip

## Explanation
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


