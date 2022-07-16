#! /usr/bin/env python3
"""
Module that contain the operations for the projects, meaning function that create array, control proprieties or
print out information.
This script call for no other module in this project and is the lowest level.

The function of this script are used in main.py and solvers.py
"""

# import
from itertools import permutations

import matplotlib.pyplot as plt
import numpy as np


def get_k(i, j, nnx):
    """
    Get the index k for a given i and j (flatten index)
    @param i: int, row index
    @param j: int, column index
    @param nnx: int, size of the domain
    @return: int, k index
    """
    return j * nnx + i


def diag_dominant(A, b, perm=True):
    """
    Check if matrix A is diagonally dominant (and permute b accordingly)
    If permutations is True, then check all permutations of the matrix.
    :param A: np.array 2D to be checked for diagonally dominance
    :param b: np.array 1D
    :param perm: bool with default = True, if True, check all permutations of the matrix
    :return A: np.array 2D, permuted if necessary
    :return b: np.array 1D, permuted the same way as A
    """
    test_diag = 0  # initialize
    test_zero = 0  # initialize
    if not perm:
        # check if matrix is diagonally dominant (without permutations)
        if np.linalg.det(A) >= 0:  # check if matrix is diagonally dominant
            test_diag = 1
        if np.all(np.diag(A) != 0):  # check if matrix has a zero on the diagonal
            test_zero = 1
        if test_zero == 1 and test_diag == 1:
            return A, b
    elif perm and np.linalg.det(A) < 0:
        # check if any permutation of the matrix is diagonally dominant
        # condition: permutation test is ask (perm=True) and matrix is not diagonally dominant already
        Ab_array = np.append(A, b.reshape(-1, 1), axis=1)  # append b to A
        Ab_tuple = tuple(map(tuple, Ab_array))  # convert array to tuple
        p = set(permutations(Ab_tuple))  # create permutations of A
        for item in list(p):
            item_array = np.asarray(item)  # convert tuple to array
            A_perm = item_array[:, :-1]
            b_perm = item_array[:, -1]
            if np.linalg.det(A_perm) >= 0:  # check if matrix is diagonally dominant
                test_diag = 1
            if np.all(np.diag(A_perm) != 0):  # check if matrix has a zero on the diagonal
                test_zero = 1
            if test_zero == 1 and test_diag == 1:
                return A_perm, b_perm
    else:
        return A, b
    if test_diag == 0:  # if all permutation are not diagonally dominant
        raise ValueError("Matrix permutations are not diagonally dominant!")
    elif test_zero == 0:  # if all permutation have
        raise ValueError("Matrix permutations have a zero on the diagonal!")


def matrix_filling_implicit(nnx, s):
    """
    Fill a matrix A from top left to bottom right (line by line) to use in implicit solvers.
    @param nnx: integer, size of the domain used in the solver (also the size of the A matrix)
    @param s: float, parameter for implicit method with the form s = kappa * dt / h²
        kappa is the thermal conductivity [W/mK]
        dt is the time step [s]
        h is the space step [m]
    @return: np.array 2D, A matrix
    """
    A = np.zeros((nnx, nnx))  # initialize matrix
    for i in range(1, nnx - 1):  # fill rows
        for j in range(0, nnx):  # fill columns
            if i == j:
                A[i, j] = 1 + 2 * s
            elif i == j - 1:
                A[i, j] = -s
            elif i == j + 1:
                A[i, j] = -s
    A[0, 0] = 1.0  # left boundary
    A[-1, -1] = 1.0  # right boundary
    return A


def matrix_filling_implicit_2d(nnx, nny, sx, sy):
    """
    Fill a matrix A from top left to bottom right (line by line) to use in implicit solvers.
    @param nnx: integer, size of the domain used in the solver on x-axis
    @param nny: integer, size of the domain used in the solver on y-axis
    @param sx: float, parameter for implicit method with the form sx = kappa * dt / hx²
    @param sy: float, parameter for implicit method with the form sy = kappa * dt / hy²
            kappa is the thermal conductivity [W/mK]
            dt is the time step [s]
            h is the space step [m]
    @return: np.array 2D, A matrix
    """
    nnp = nnx * nny
    A = np.eye(nnp)  # initialize matrix
    for j in range(0, nny):  # fill rows
        for i in range(0, nnx):  # fill columns
            if i == 0 or j == 0 or i == nnx - 1 or j == nny - 1:
                A[get_k(i, j, nnx), get_k(i, j, nnx)] = 1
            else:
                A[get_k(i, j, nnx), get_k(i, j, nnx)] = 1 + 2 * sx + 2 * sy  # diagonal
                A[get_k(i - 1, j, nnx) + 1, get_k(i - 1, j, nnx)] = -sx  # right
                A[get_k(i + 1, j, nnx) - 1, get_k(i + 1, j, nnx)] = -sx  # left
                A[get_k(i, j - 1, nnx) + nnx, get_k(i, j - 1, nnx)] = -sy  # top
                A[get_k(i, j + 1, nnx) - nnx, get_k(i, j + 1, nnx)] = -sy  # bottom
    return A


def matrix_filling_explicit(nnx, s):
    """
    Fill a matrix A from top left to bottom right (line by line) to use in explicit solvers (e.g. RHS of crank-nicolson).
    @param nnx: integer, size of the domain used in the solver (also the size of the A matrix)
    @param s: float, parameter for implicit method with the form s = kappa * dt / h²
        kappa is the thermal conductivity [W/mK]
        dt is the time step [s]
        h is the space step [m]
    @return: np.array 2D, A matrix
    """
    A = np.zeros((nnx, nnx))  # initialize matrix
    for i in range(1, nnx - 1):  # fill rows
        for j in range(0, nnx):  # fill columns
            if i == j:
                A[i, j] = 1 - 2 * s
            elif i == j - 1:
                A[i, j] = s
            elif i == j + 1:
                A[i, j] = s
    A[0, 0] = 1.0  # left boundary
    A[-1, -1] = 1.0  # right boundary
    return A


def results_print(method, err_iter, iter_count, max_iter, time):
    """
    Print out the results of the solver.
    @param method: string, name of the method
    @param err: float, error calculate max(|T_n - T_n+1|)
    @param iter_count: int, number of iterations done
    @param max_iter: int, maximum number of iterations
    """
    print('Method: ' + method + '\n')
    if iter_count == max_iter:
        print('Maximum number of iterations reached! \n')
    else:
        print('Number of iterations: ' + str(iter_count) + '\n')
    print('Final error: ' + str(round(err_iter, 4)) + '\n')
    print('Time: ' + str(round(time, 4)) + '\n')
    print('Time to compute solution: ' + str(round(time, 2)) + 's \n')
    print('- -\n')


def results_plot_1d(path_output, x, T_init, T_final, label_1, T_final_2=None, label_2=None, T_final_3=None, label_3=None):
    """
    Plot the initial and final temperature profiles.
    @param x: np.array 1D, x-coordinates of the domain
    @param T_init: np.array 1D, initial temperature profile
    @param T_final: np.array 1D, final temperature profile
    @param T_final_2: optional np.array 1D, final temperature profile for a 2nd method
    @param T_final_2: optional np.array 1D, final temperature profile for a 3rd method
    """
    plt.figure(figsize=(6, 5))
    plt.title('Tempearture profile of the domain')
    plt.plot(x, T_init, color='black', linewidth=2, label='Initial')
    plt.plot(x, T_final, color='blue', linewidth=2, label=label_1)
    plt.plot(x, T_final_2, color='red', linewidth=2, label=label_2) if T_final_2 is not None else None
    plt.plot(x, T_final_3, color='orange', linewidth=2, label=label_3) if T_final_3 is not None else None
    plt.xlabel('x [m]')
    plt.ylabel('T [C]')
    plt.legend()
    plt.savefig(path_output)
    plt.show()


def results_plot_2d(path_output, T, x, y, nnx, nny, title):
    """
    Plot the temperature profile T in 3D.
    @param path_output:
    @param T: np.array, temperature profile to plot (z-axis)
    @param x: np.array, x-axis
    @param y: np.array, y-axis
    @param nnx: int, number of points in x-direction
    @param nny: int, number of points in y-direction
    """
    plt.figure(figsize=(6, 5))
    ax = plt.axes(projection='3d')
    ax.set_title(title)
    ax.plot_surface(x.reshape((nny, nnx)), y.reshape((nny, nnx)), T.reshape((nny, nnx)), color='blue')
    ax.set_xlabel('x [m]')
    ax.set_ylabel('y [m]')
    ax.set_zlabel('T [K]')
    plt.savefig(path_output)
    plt.show()
