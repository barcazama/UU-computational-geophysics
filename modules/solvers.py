#! /usr/bin/env python3
"""
Module that contain the solver for the projects.
This module call for function in operations.py.
This module is call by methods.py
"""

# import
import numpy as np

from modules.operations import diag_dominant


def jacobi(A, b, tol=0.1, max_iter=1000, perm=True):
    """
    Jacobi solver.
    Equation: x_new = diag_inv * (b - (L+U) * x_old)
    :param A: matrix of size n x n
    :param b: vector of size n (make sure it's not a matrix)
    :param tol: tolerance number
    :param max_iter: maximum number of iterations
    :return: solution vector x
    """

    A, b = diag_dominant(A, b, perm)  # condition for convergence

    # initialize
    D = np.diag(np.diag(A))  # create diagonal matrix
    D_inv = np.linalg.inv(D)  # create inverse of diagonal matrix
    L = np.tril(A, -1)  # create lower triangular matrix
    U = np.triu(A, 1)  # create upper triangular matrix
    err = 10  # initialize error
    iter_count = 0  # initialize iteration count
    x = np.zeros(len(b))  # initialize x

    # iterate
    while err > tol or iter_count < max_iter:
        x_old = x.copy()  # save old x
        x_new = D_inv @ (b - (L + U) @ x_old)  # calculate new x
        x = x_new.copy()  # update x
        err = np.linalg.norm(x - x_old)  # calculate error
        iter_count += 1  # update iteration count
    return x


def gauss_seidel(A, b, tol=0.1, max_iter=1000, perm=True):
    """
    gauss seidel solver.
    Equation: x_new = (L+diag)_inv * (b - U * x_old)
    :param A: matrix of size n x n
    :param b: vector of size n (make sure it's not a matrix)
    :param tol: tolerance number
    :param max_iter: maximum number of iterations
    :param perm: if true, try all matrix permutation
    :return: solution vector x
    """

    A, b = diag_dominant(A, b, perm)  # condition for convergence

    # initialize
    D = np.diag(np.diag(A))  # create diagonal matrix
    L = np.tril(A, -1)  # create lower triangular matrix
    U = np.triu(A, 1)  # create upper triangular matrix
    err = 1  # initialize error
    iter_count = 0  # initialize iteration count
    x = np.zeros(len(b))  # initialize x

    # iterate
    while err > tol or iter_count < max_iter:
        x_old = x.copy()  # save old x
        x_new = np.linalg.inv(L + D) @ (b - U @ x_old)  # calculate new x
        x = x_new.copy()  # update x
        err = np.linalg.norm(x - x_old)  # calculate error
        iter_count += 1  # update iteration count
    return x


def ftcs(b_old, u, dt, h, nnx):
    """
    Solving advection using FTCS solver (Forward time centered space)
    @param b_old: np.array of size nnx containing the initial or previous iteration temperature values
    @param u: float, velocity [m/s]
    @param dt: float, time step [s]
    @param h: foat, space step [m]
    @param nnx: int, number of nodes in space
    @return b_new: np.array of size nnx containing the new iteration temperature values
    """
    b_new = np.zeros(nnx)  # initialise b_new array
    for i in range(1, nnx - 1):
        b_new[i] = b_old[i] - u * dt / 2 / h * (b_old[i + 1] - b_old[i - 1])
    return b_new


def lax(b_old, u, dt, h, nnx):
    """
    Solving advection using Lax-Friedrichs solver
    @param b_old: np.array of size nnx containing the initial or previous iteration temperature values
    @param u: float, velocity [m/s]
    @param dt: float, time step [s]
    @param h: foat, space step [m]
    @param nnx: int, number of nodes in space
    @return b_new: np.array of size nnx containing the new iteration temperature values
    """
    b_new = np.zeros(nnx)  # initialise b_new array
    for i in range(1, nnx - 1):
        b_new[i] = 0.5 * ((b_old[i + 1] + b_old[i - 1]) - (u * dt / h) * (b_old[i + 1] - b_old[i - 1]))
    return b_new


def midpoint_rule(f, a, b, n):
    """
    Midpoint rule for numerical integration
    f: function to be integrated
    a: lower limit
    b: upper limit
    n: number of subintervals
    """
    h = (b - a) / (n - 1)  # compute subinterval length
    x = np.linspace(a, b, n)  # compute subinterval points
    I = 0
    for i in range(n):
        try:
            I += f((x[i] + x[i + 1]) / 2)
        except:
            I += (x[i] + x[i + 1]) / 2
    I *= h
    return I


def trapeze_rule(f, a, b, n):
    """ TODO: unfinished
    Trapeze rule for numerical integration
    f: function to be integrated
    a: lower limit
    b: upper limit
    n: number of subintervals
    """
    h = (b - a) / (n - 1)  # compute subinterval length
    x = np.linspace(a, b, n)  # compute subinterval points
    I = f(x[0]) + f(x[-1])
    for i in range(1, n - 1):
        I += f(x[i]) + f(x[i + 1])
    I *= h / 2
    return I
