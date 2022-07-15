#! /usr/bin/env python3
"""
Module that contain the methods for the projects, meaning the complete methods.
Situated at the highest level of the tree, it will only be called by main.py
"""

# import
import numpy as np
from scipy.linalg import solve
from scipy.sparse.linalg import spsolve
from scipy.sparse import csc_matrix

import operations
import solvers


def fde_explicit(T_init, kappa, nnx, h, dt, tol, max_iter):
    """
    Solving 1D diffusion with the explicit method
    @param T_init: np.array, initial temperature
    @param kappa: float, thermal conductivity [W/m/K]
    @param nnx: int, number of grid point in x-axis
    @param h: float, grid spacing [m]
    @param dt: float, time interval [s]
    @param tol: float, stop condition (error tolerance)
    @param max_iter: int, maximum number of iterations
    @return method_name: str, method used for the solution
    @return T: np.array, temperature at each grid point at the end of the simulation
    @return iter_count: int, number of iterations used for the solution
    @return err: float, error of the solution
    """
    err = 10  # initialize error
    iter_count = 0  # initialize iteration counter
    Told = T_init.copy()  # initial temperature array
    Tnew = T_init.copy()  # initialize temperature array
    while err > tol and iter_count < max_iter:
        for i in range(1, nnx - 1):
            Tnew[i] = Told[i] + dt * kappa * (Told[i - 1] - 2 * Told[i] + Told[i + 1]) / (h ** 2)
        err = np.max(np.abs(Told - Tnew))  # compute error
        Told[:] = Tnew[:]  # update Told
        iter_count += 1  # increment iteration counter
    method_name = 'explicit'
    T = Tnew.copy()
    return method_name, T, iter_count, err


def fde_implicit_scipy_1d(T_init, kappa, dt, h, tol, max_iter):
    """
    Solve the heat equation using the scipy python library.
    @param s: float, parameter for implicit method with the form s = kappa * dt / h²
        kappa is the thermal conductivity [W/mK]
        dt is the time step [s]
        h is the space step [m]
    @param T_init: np.array, initial temperature
    @param tol: float, stop condition (error tolerance)
    @param max_iter: int, maximum number of iterations
    @return method_name: str, method_name used for the solution
    @return T: np.array, temperature at each grid point at the end of the simulation
    @return iter_count: int, number of iterations used for the solution
    @return err: float, error of the solution
    """
    nnx = len(T_init)  # number of grid points in x direction
    s = kappa * dt / (h ** 2)  # parameter for implicit method
    A = operations.matrix_filling_implicit(nnx, s)  # filling matrix A for implicit method
    err = 10  # initialise error to start loop
    iter_count = 0  # initialise iteration counter
    b = T_init.copy()  # initial RHS with initial temperature profile
    while err > tol and iter_count < max_iter:
        b_old = b.copy()
        b_new = solve(A, b)  # solve system of linear equations using jacobi method
        b = b_new.copy()  # update Temperature profile forward in time
        err = np.max(np.abs(b_old - b_new))  # calculate error between b_old and b_new
        iter_count += 1  # iteration counter update
    method_name = 'scipy'  # storing method used for the solution
    T = b  # storing solution
    return method_name, T, iter_count, err


def fde_implicit_jacobi_1d(T_init, kappa, dt, h, tol, max_iter):
    """
    Solve the heat equation using the Jacobi method.
    @param s: float, parameter for implicit method with the form s = kappa * dt / h²
        kappa is the thermal conductivity [W/mK]
        dt is the time step [s]
        h is the space step [m]
    @param T_init: np.array, initial temperature
    @param tol: float, stop condition (error tolerance)
    @param max_iter: int, maximum number of iterations
    @return method_name: str, method_name used for the solution
    @return T: np.array, temperature at each grid point at the end of the simulation
    @return iter_count: int, number of iterations used for the solution
    @return err: float, error of the solution
    """
    nnx = len(T_init)  # number of grid points in x direction
    s = kappa * dt / (h ** 2)  # parameter for implicit method
    A = operations.matrix_filling_implicit(nnx, s)  # filling matrix A for implicit method
    err = 10  # initialise error to start loop
    iter_count = 0  # initialise iteration counter
    b = T_init.copy()  # initial RHS with initial temperature profile
    while err > tol and iter_count < max_iter:
        b_old = b.copy()
        b_new = solvers.jacobi(A, b, tol, max_iter)  # solve system of linear equations using jacobi method
        b = b_new.copy()  # update Temperature profile forward in time
        err = np.max(np.abs(b_old - b_new))  # calculate error between b_old and b_new
        iter_count += 1  # iteration counter update
    method_name = 'jacobi'  # storing method used for the solution
    T = b  # storing solution
    return method_name, T, iter_count, err


def fde_implicit_gauss_1d(T_init, kappa, dt, h, tol, max_iter):
    """
    Solve the heat equation using the gauss-seidel method.
    @param s: float, parameter for implicit method with the form s = kappa * dt / h²
        kappa is the thermal conductivity [W/mK]
        dt is the time step [s]
        h is the space step [m]
    @param T_init: np.array, initial temperature
    @param tol: float, stop condition (error tolerance)
    @param max_iter: int, maximum number of iterations
    @return method_name: str, method_name used for the solution
    @return T: np.array, temperature at each grid point at the end of the simulation
    @return iter_count: int, number of iterations used for the solution
    @return err: float, error of the solution
    """
    nnx = len(T_init)  # number of grid points in x direction
    s = kappa * dt / (h ** 2)  # parameter for implicit method
    A = operations.matrix_filling_implicit(nnx, s)  # filling matrix A for implicit method
    err = 10  # initialise error to start loop
    iter_count = 0  # initialise iteration counter
    b = T_init.copy()  # initial RHS with initial temperature profile
    while err > tol and iter_count < max_iter:
        b_old = b.copy()
        b_new = solvers.gauss_seidel(A, b, tol, max_iter)  # solve system of linear equations using jacobi method
        b = b_new.copy()  # update Temperature profile forward in time
        err = np.max(np.abs(b_old - b_new))  # calculate error between b_old and b_new
        iter_count += 1  # iteration counter update
    method_name = 'gauss-seidel'  # storing method used for the solution
    T = b  # storing solution
    return method_name, T, iter_count, err


def fde_implicit_crank_1d(T_init, kappa, dt, h, tol, max_iter):
    """
    Solve the heat equation using the Crank-Nicolson method.
    @param s: float, parameter for implicit method with the form s = kappa * dt / h²
        kappa is the thermal conductivity [W/mK]
        dt is the time step [s]
        h is the space step [m]
    @param T_init: np.array, initial temperature
    @param tol: float, stop condition (error tolerance)
    @param max_iter: int, maximum number of iterations
    @return method_name: str, method_name used for the solution
    @return T: np.array, temperature at each grid point at the end of the simulation
    @return iter_count: int, number of iterations used for the solution
    @return err: float, error of the solution
    """
    nnx = len(T_init)  # number of grid points in x direction
    s = kappa * dt / (2 * h ** 2)  # parameter for implicit method
    A_implicit = operations.matrix_filling_implicit(nnx, s)  # filling matrix A (implicit)
    A_explicit = operations.matrix_filling_explicit(nnx, s)  # filling matrix A (explicit)
    A_implicit_inv = np.linalg.inv(A_implicit)  # inverse of A_implicit
    err = 10  # initialise error to start loop
    iter_count = 0  # initialise iteration counter
    b = T_init.copy()  # initial RHS with initial temperature profile
    while err > tol and iter_count < max_iter:
        b_old = b.copy()
        b_new = A_implicit_inv @ (A_explicit @ b_old)
        b = b_new.copy()  # update Temperature profile forward in time
        err = np.max(np.abs(b_old - b_new))  # calculate error between b_old and b_new
        iter_count += 1  # iteration counter update
    method_name = 'crank-nicolson'  # storing method used for the solution
    T = b  # storing solution
    return method_name, T, iter_count, err


def fde_advection_ftcs_1d(T_init, u, dt, h, tol, max_iter):
    """
    Solve the advection equation using the FTCS method.
    @param T_init: np.array, initial temperature
    @param u: float, velocity [m/s]
    @param dt: flaot, time step [s]
    @param h: float, space step [m]
    @param tol: float, stop condition (error tolerance)
    @param max_iter: int, maximum number of iterations
    @return method_name: str, method_name used for the solution
    @return T: np.array, temperature at each grid point at the end of the simulation
    @return iter_count: int, number of iterations used for the solution
    @return err: float, error of the solution
    """
    nnx = len(T_init)  # number of grid points in x direction
    err = 10  # initialise error to start loop
    iter_count = 0  # initialise iteration counter
    b = T_init.copy()  # initial RHS with initial temperature profile
    while err > tol and iter_count < max_iter:
        b_old = b.copy()
        b_new = solvers.ftcs(b_old, u, dt, h, nnx)
        b = b_new.copy()
        err = np.max(np.abs(b_old - b_new))  # calculate error between b_old and b_new
        iter_count += 1  # iteration counter update
    method_name = 'FTCS'
    T = b  # storing solution
    return method_name, T, iter_count, err


def fde_advection_lax_1d(T_init, u, dt, h, tol, max_iter):
    """
    Solve the advection equation using the Lax-Friedrichs method.
    @param T_init: np.array, initial temperature
    @param u: float, velocity [m/s]
    @param dt: flaot, time step [s]
    @param h: float, space step [m]
    @param tol: float, stop condition (error tolerance)
    @param max_iter: int, maximum number of iterations
    @return method_name: str, method_name used for the solution
    @return T: np.array, temperature at each grid point at the end of the simulation
    @return iter_count: int, number of iterations used for the solution
    @return err: float, error of the solution
    """
    nnx = len(T_init)  # number of grid points in x direction
    err = 10  # initialise error to start loop
    iter_count = 0  # initialise iteration counter
    b = T_init.copy()  # initial RHS with initial temperature profile
    while err > tol and iter_count < max_iter:
        b_old = b.copy()
        b_new = solvers.lax(b_old, u, dt, h, nnx)
        b = b_new.copy()
        err = np.max(np.abs(b_old - b_new))  # calculate error between b_old and b_new
        iter_count += 1  # iteration counter update
    method_name = 'Lax-Friedrichs'
    T = b  # storing solution
    return method_name, T, iter_count, err


def fde_diffusion_ana_init_2d(T_max, T0, kappa, sigma, xcoords, ycoords, nnx, nny, max_iter):
    """
    Solve the diffusion equation using the analytical for the initial conditions (t=0).
    @param T_max: float, maximum temperature
    @param T0: float, initial temperature
    @param kappa: float, thermal conductivity [W/mK]
    @param sigma: float, surface heat flux [W/m²]
    @param xcoords: np.array, x coordinates of the grid points
    @param ycoords: np.array, y coordinates of the grid points
    @param nnx: int, number of grid points in x direction
    @param nny: int, number of grid points in y direction
    @param max_iter: int, maximum number of iterations
    @return method_name: str, method_name used for the solution
    @return T: np.array, temperature at each grid point at the end of the simulation
    @return iter_count: int, number of iterations used for the solution
    """
    iter_count = 0  # initialise iteration counter
    b_new = np.zeros(nnx * nny)  # initialise b_new
    while iter_count < max_iter:
        t = 0  # time = 0 because we just want initial
        for j in range(0, nny):
            for i in range(0, nnx):
                k = j * nnx + i
                b_new[k] = T0 + T_max / (1 + 4 * t * kappa / (sigma ** 2)) * np.exp(
                    -(xcoords[k] ** 2 + ycoords[k] ** 2) / (sigma ** 2 + 4 * t * kappa))
        b = b_new.copy()
        iter_count += 1  # iteration counter update
    method_name = 'Analytical'
    T = b  # storing solution
    return method_name, T, iter_count


def fde_diffusion_explicit_2d(T_init, hx, hy, kappa, Q, nnx, nny, dt, tol, max_iter):
    """
    Solve the diffusion equation using the explicit method.
    @param T_init: np.array, initial temperature
    @param hx: float, space step in x direction [m]
    @param hy: float, space step in y direction [m]
    @param kappa: float, thermal conductivity [W/(m*K)]
    @param Q: float, heat source [W/m^2]
    @param nnx: int, number of grid points in x direction
    @param nny: int, number of grid points in y direction
    @param dt: float, time step [s]
    @param tol: float, stop condition (error tolerance)
    @param max_iter: int, maximum number of iterations
    @return method_name: str, method_name used for the solution
    @return T: np.array, temperature at each grid point at the end of the simulation
    @return iter_count: int, number of iterations used for the solution
    @return err: float, error of the solution
    """
    sx = kappa * dt / (hx ** 2)  # parameter for explicit method
    sy = kappa * dt / (hy ** 2)  # parameter for explicit method
    err = 10  # initialise error to start loop
    iter_count = 0  # initialise iteration counter
    b = T_init.copy()  # initialise b_old
    b_new = np.zeros(nnx * nny)  # initialise b_new
    # compute solution
    while err > tol or iter_count < max_iter:
        b_old = b.copy()
        for j in range(0, nny):
            for i in range(0, nnx):
                k = j * nnx + i  # index of current grid point for when array is flattened (i and j together)
                if i == 0 or j == 0 or i == nnx - 1 or j == nny - 1:
                    b_new[k] = T_init[k]
                else:
                    b_new[k] = b_old[k] + sx * (
                            b_old[j * nnx + i - 1] - 2 * b_old[j * nnx + i] + b_old[j * nnx + i + 1]) + sy * (
                                       b_old[(j - 1) * nnx + i] - 2 * b_old[j * nnx + i] + b_old[
                                   (j + 1) * nnx + i]) + Q * dt
        b = b_new.copy()
        err = np.max(np.abs(b_old - b_new))  # calculate error between b_old and b_new
        iter_count += 1  # iteration counter update
    method_name = 'Explicit'
    T = b  # storing solution
    return method_name, T, iter_count, err


def fde_diffusion_implicit_2d(T_init, kappa, dt, hx, hy, nnx, nny, tol, max_iter):
    """
    Solve the diffusion equation in 2D using the implicit method.
    @param T_init: np.array, initial temperature
    @param kappa: float, thermal conductivity [W/(m*K)]
    @param dt: float, time step [s]
    @param hx: float, space step in x direction [m]
    @param hy: float, space step in y direction [m]
    @param nnx: int, number of grid points in x direction
    @param nny: int, number of grid points in y direction
    @param dt: float, time step [s]
    @param tol: float, stop condition (error tolerance)
    @param max_iter: int, maximum number of iterations
    @return method_name: str, method_name used for the solution
    @return T: np.array, temperature at each grid point at the end of the simulation
    @return iter_count: int, number of iterations used for the solution
    @return err: float, error of the solution
    """
    sx = kappa * dt / (hx ** 2)  # parameter for implicit method
    sy = kappa * dt / (hy ** 2)  # parameter for implicit method
    A = operations.matrix_filling_implicit_2d(nnx, nny, sx, sy)  # filling matrix A for implicit method
    err = 10  # initialise error to start loop
    iter_count = 0  # initialise iteration counter
    b = T_init.copy()  # initial RHS with initial temperature profile
    while err > tol and iter_count < max_iter:
        b_old = b.copy()
        b_new = solvers.jacobi(A, b, tol, max_iter)  # solve system of linear equations using jacobi method
        b = b_new.copy()  # update Temperature profile forward in time
        err = np.max(np.abs(b_old - b_new))  # calculate error between b_old and b_new
        iter_count += 1  # iteration counter update
    method_name = 'Implicit'  # storing method used for the solution
    T = b  # storing solution
    return method_name, T, iter_count, err


def fde_advection_2d(T_init, u, v, hx, hy, nnx, nny, dt, max_iter):
    """
    Solve the advection equation in 2D.
    @param T_init:  np.array, initial temperature
    @param u: np.array, velocity[m/s]
    @param v: np.array, velocity[m/s]
    @param hx: float, space step in x direction [m]
    @param hy: float, space step in y direction [m]
    @param nnx: int, number of grid points in x direction
    @param nny: int, number of grid points in y direction
    @param dt: float, time step [s]
    @param tol: float, stop condition (error tolerance)
    @param max_iter: int, maximum number of iterations
    @return method_name: str, method_name used for the solution
    @return T: np.array, temperature at each grid point at the end of the simulation
    @return iter_count: int, number of iterations used for the solution
    @return err: float, error of the solution
    """
    err = 10  # initialise error to start loop
    iter_count = 0  # initialise iteration counter
    b = T_init.copy()  # initial RHS with initial temperature profile
    b_new = np.zeros(nnx * nny)  # initialise b_new
    for n in range(0, max_iter):
        b_old = b.copy()
        for j in range(0, nny):
            for i in range(0, nnx):
                k = j * nnx + i
                if i == 0 or j == 0 or i == nnx - 1 or j == nny - 1:
                    b_new[k] = T_init[k]
                else:
                    b_new[k] = b_old[k] - dt * (
                            u[k] * (b_old[nnx * j + i + 1] - b_old[nnx * j + i - 1]) / (2 * hx) + v[k] * (
                            b_old[nnx * (j + 1) + i] - b_old[nnx * (j - 1) + i]) / (2 * hy))
        b = b_new.copy()  # update Temperature profile forward in time
        err = np.max(np.abs(b_old - b_new))  # calculate error between b_old and b_new
        iter_count += 1  # iteration counter update
    method_name = 'Advection'  # storing method used for the solution
    T = b  # storing solution
    return method_name, T, iter_count, err


def fem_diffusion_1d(T_init, kappa, rho, cp, hx, nnx, dt, tol, max_iter, T_ana):
    """
    Solve the diffusion equation in 1D using the finite element method.
    @param T_init: np.array, initial temperature
    @param kappa: float, thermal conductivity [W/(m*K)]
    @param rho: float, density [kg/m^3]
    @param cp: float, specific heat capacity [J/(kg*K)]
    @param hx: float, space step in x direction [m]
    @param nnx: int, number of grid points in x direction
    @param dt: float, time step [s]
    @param tol: float, stop condition (error tolerance)
    @param max_iter: int, maximum number of iterations
    @param T_ana: np.array, analytical temperature profile
    @return method_name: str, method_name used for the solution
    @return T: np.array, temperature at each grid point at the end of the simulation
    @return iter_count: int, number of iterations used for the solution
    @return err: float, error of the solution
    """
    nelx = nnx - 1 # number of elements in x direction
    err_ana = 10 # initialise error to start loop
    iter_count = 0 # initialise iteration counter
    T = T_init.copy() # initialise RHS (array to store solution)
    while err_ana > tol and iter_count < max_iter:
        A = np.zeros((nnx, nnx))  # initialise matrix A
        b = np.zeros(nnx)  # initialise RHS with initial temperature profile
        for e in range(0, nelx):  # loop over elements
            Me = hx / 3 * rho * cp * np.array([[1, 0.5], [0.5, 1]]) # elemental mass matrix
            Kde = kappa / hx * np.array([[1, -1], [-1, 1]]) # elemental
            Ae = Me + Kde * dt # elemental A matrix
            A[e, e] += Ae[0, 0] # assemble matrix A
            A[e + 1, e] += Ae[1, 0]
            A[e, e + 1] += Ae[0, 1]
            A[e + 1, e + 1] += Ae[1, 1]
            b[e] += (Me[0, 0] + Me[0, 1]) * T[e] # assemble array b
            b[e + 1] += (Me[1, 0] + Me[1, 1]) * T[e + 1]
        A[0, :] = 0 # Boundary conditions
        A[nnx - 1, :] = 0
        A[0, 0] = 1
        A[nnx - 1, nnx - 1] = 1
        b[0] = T_init[0] # boundary condition
        b[nnx - 1] = T_init[nnx - 1] # boundary condition
        T_new = spsolve(csc_matrix(A), b) # solve system of linear equations using scipy for computing misfit
        T = T_new.copy()
        err_ana = np.max(np.abs(T_ana - T_new)) # misfit between calculated and analytical solution
        iter_count += 1 # iteration counter update
    method_name = 'FEM_diffusion' # storing method used for the solution
    return method_name, T, iter_count, err_ana


def fem_advection_1d(T_init, u, kappa, rho, cp, hx, nnx, dt, max_iter, solver):
    nelx = nnx -1 # number of elements in x direction
    iter_count = 0 # initialise iteration counter
    T = T_init.copy() # initialise RHS (array to store solution)
    while iter_count < max_iter:
        A = np.zeros((nnx, nnx))  # initialise matrix A
        b = np.zeros(nnx)  # initialise RHS with initial temperature profile
        for e in range(0, nelx):  # loop over elements
            Me = hx / 3 * rho * cp * np.array([[1, 0.5], [0.5, 1]]) # elemental mass matrix
            # Kde = kappa / hx * np.array([[1, -1], [-1, 1]]) # elemental
            Kde = np.zeros((2, 2))
            Kae = rho * cp * u * np.array([[-0.5, 0.5], [-0.5, 0.5]])
            Fe = np.zeros(2)
            T_current = [T[e], T[e + 1]]
            if solver == 'explicit':
                Ae = Me
                be = Me @ T_current + Fe * dt
            elif solver == 'implicit':
                Ae = Me + (Kde + Kae) * dt
                be = Me @ T_current + Fe * dt
            elif solver == 'crank_nicolson':
                # TODO: implement crank_nicolson
                pass
            A[e, e] += Ae[0, 0] # assemble matrix A
            A[e + 1, e] += Ae[1, 0]
            A[e, e + 1] += Ae[0, 1]
            A[e + 1, e + 1] += Ae[1, 1]
            b[e] += be[0] # assemble array b
            b[e + 1] += be[1]
        A[0, :] = 0 # Boundary conditions
        A[nnx - 1, :] = 0
        A[0, 0] = 1
        A[nnx - 1, nnx - 1] = 1
        b[0] = T_init[0] # boundary condition
        b[nnx - 1] = T_init[nnx - 1] # boundary condition
        T_new = spsolve(csc_matrix(A), b) # solve system of linear equations using scipy for computing misfit
        T = T_new.copy()
        iter_count += 1 # iteration counter update
    method_name = 'FEM_advection' # storing method used for the solution
    return method_name, T, iter_count


def fem_diffusion_2d(T_init, kappa, rho, cp, hx, hy, nnx, nny, dt, tol, max_iter, T_ana):
    """
    Solve the diffusion equation in 1D using the finite element method.
    @param T_init: np.array, initial temperature
    @param kappa: float, thermal conductivity [W/(m*K)]
    @param rho: float, density [kg/m^3]
    @param cp: float, specific heat capacity [J/(kg*K)]
    @param hx: float, space step in x direction [m]
    @param nnx: int, number of grid points in x direction
    @param dt: float, time step [s]
    @param tol: float, stop condition (error tolerance)
    @param max_iter: int, maximum number of iterations
    @param T_ana: np.array, analytical temperature profile
    @return method_name: str, method_name used for the solution
    @return T: np.array, temperature at each grid point at the end of the simulation
    @return iter_count: int, number of iterations used for the solution
    @return err: float, error of the solution
    """
    # TODO: work in progress
    nelx = nnx - 1 # number of elements in x direction
    nely = nny - 1
    err_ana = 10 # initialise error to start loop
    iter_count = 0 # initialise iteration counter
    T = T_init.copy() # initialise RHS (array to store solution)
    while err_ana > tol and iter_count < max_iter:
        A = np.zeros((nnx, nnx))  # initialise matrix A
        b = np.zeros(nnx)  # initialise RHS with initial temperature profile
    #     for ex in range(0, nelx):
    #         for ey in range(0, nely):
    #             # loop over elements
    #             Me = hx / 3 * rho * cp * np.array([[1, 0.5], [0.5, 1]]) # elemental mass matrix
    #             Kde = kappa / hx * np.array([[1, -1], [-1, 1]]) # elemental
    #             Ae = Me + Kde * dt # elemental A matrix
    #             A[e, e] += Ae[0, 0] # assemble matrix A
    #             A[e + 1, e] += Ae[1, 0]
    #             A[e, e + 1] += Ae[0, 1]
    #             A[e + 1, e + 1] += Ae[1, 1]
    #             b[e] += (Me[0, 0] + Me[0, 1]) * T[e] # assemble array b
    #             b[e + 1] += (Me[1, 0] + Me[1, 1]) * T[e + 1]
    #     A[0, :] = 0 # Boundary conditions
    #     A[nnx - 1, :] = 0
    #     A[0, 0] = 1
    #     A[nnx - 1, nnx - 1] = 1
    #     b[0] = T_init[0] # boundary condition
    #     b[nnx - 1] = T_init[nnx - 1] # boundary condition
    #     T_new = spsolve(csc_matrix(A), b) # solve system of linear equations using scipy for computing misfit
    #     T = T_new.copy()
    #     err_ana = np.max(np.abs(T_ana - T_new)) # misfit between calculated and analytical solution
    #     iter_count += 1 # iteration counter update
    # method_name = 'FEM_diffusion_2D' # storing method used for the solution
    # return method_name, T, iter_count, err_ana
