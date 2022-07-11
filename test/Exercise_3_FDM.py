#! /usr/bin/env python3
"""
Exercise 3: Finite Difference Method - Creation of the Jacobi Solver
====================================================================
solvers are located in modules/solvers.py

ANSWERS TO FDM-3 QUESTIONS:
--------------------------
With the given matrix, the solver does not converge (= explode).
Multiplying the matrix diagonal by 10 will make the solver diverge slower.
This is because the solver is not able to handle the case where the matrix A is not diagonally dominant.
The gauss-seidel solver is between 2 to 10 times slower than the jacobi solver.
"""

# import
import sys
import timeit

import numpy as np
from scipy.linalg import solve

from modules import solvers


def main_function():
    """
    Main function.

    to test the jacobi solver
    """
    tol = 1e-4  # tolerance number
    random_number = False  # set if you want to use random numbers or preset
    n = 10  # number of rows and columns in case of random number

    # generate
    if random_number:
        A = np.random.rand(n, n)
        b = np.random.rand(n)
    else:
        A = np.array([[50, 3, -2], [3, 50, 6], [2, 4, 30]])
        b = np.array([5, 7, 8])

    # solve
    for i in range(3):
        start = timeit.default_timer()
        if i==0:
            method = 'scipy'
            x = solve(A, b)
        elif i==1:
            method = 'jacobi'
            x = solvers.jacobi(A, b, tol) # solve with jacobi solver
        elif i==2:
            method = 'gauss_seidel'
            x = solvers.gauss_seidel(A, b, tol)
        stop = timeit.default_timer()
        print("Time for "+ method + " solver: ", stop - start)
        print("Solution: ", x, "\n")

if __name__ == "__main__":
    print('Starting Exercise 3...')
    print('------------------------\n')
    main_function()
    print('\n------------------------')
    print('Done!')
    sys.exit(0)
