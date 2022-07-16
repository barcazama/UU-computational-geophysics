#! /usr/bin/env python3
"""
This is the main file for the project.
It contains multiple condition that are used to select the method, the solver, the dimension and other parameters.
Here is the following main options:
    - method: diffusion and advection for FDM and FEM and multiple scheme (implicit, explicit, Crank-Nicolson)
    - solver: explicit, scipy, jacobi, gauss_seidel, crank-nicolson, FTCS, Lax-Friedrichs
    - dimension: 1D, 2D
"""
# import
import os
import sys
import timeit

import numpy as np

import methods
import operations


def main_function():
    """
    Main function that will call for all the different functions in the project, see README.md for more information.
    """

    ###### Manual Inputs ######
    # TODO: consider moving this section to yaml file to make it more human readable

    # Method Options (uncomment only one line)
    # ---------------
    # method = 'all_1d_diffusion'  # specify method: both implicit and explicit for 1D diffusion
    # method = '1d_diffusion_explicit' # specify method: explicit for 1D diffusion
    method = '1d_diffusion_implicit' # specify method: implicit for 1D diffusion
    # method = '1d_advection'  # specify method: implicit for 1D advection
    # method = 'all_2d_diffusion'  # specify method: both implicit and explicit for 2D diffusion
    # method = '2d_diffusion_explicit' # specify method: explicit for 2D diffusion
    # method = '2d_diffusion_implicit' # specify method: implicit for 2D diffusion
    # method = '2d_advection'  # specify method: implicit for 2D advection
    # method = '1d_fem_diffusion'  # specify method: finite element method for 1D diffusion
    # method = '1d_fem_advection'  # specify method: finite element method for 1D advection
    # method = '2d_fem_diffusion'  # specify method: finite element method for 2D diffusion

    # Solver Options
    # --------------
    # Solver choice for 1D diffusion - Explicit
    # solver = 'None'
    # --
    # Solver choice for 1D diffusion - Implicit
    # solver = 'all_1d_diffusion'  # all solvers for 1D diffusion
    # solver = '1d_scipy' # SciPy solver for 1D diffusion
    solver = '1d_jacobi' # Jacobi solver for 1D diffusion
    # solver = '1d_gauss_seidel' # Gauss-Seidel solver for 1D diffusion
    # solver = '1d_crank-nicolson' # Crank-Nicolson solver for 1D diffusion
    # --
    # Solver choice for 1D advection
    # solver = 'all_1d_advection'
    # solver = '1d_FTCS' # specify solver 'FTCS'
    # solver = '1d_Lax-Friedrichs'  # specify solver 'Lax-Wendroff'
    # --
    # Solver choice for 2D diffusion
    # solver = None  # there is only one solver implemented for 2D methods
    # --
    # Solver choice for 1D FEM diffusion
    # solver = None # there is only one solver implemented for 1D FEM diffusion
    # --
    # Solver choice for 1D FEM advection
    # solver = 'explicit' # not really a solver but was easier to implement a choice calling them solver (consistency)
    # solver = 'implicit'
    # solver = 'crank-nicolson'
    # --
    # Solver choice for 2D FEM
    # solver = None # there is only one solver implemented for 2D FEM

    # physical variables
    Lx = 1000  # Length of the domain [m] (use this for 1D domain as well)
    Ly = 80e3  # Height of the domain [m]
    start_x = 0  # start of the domain on x-axis [m] (use this for 1D domain as well)
    start_y = 0  # start of the domain on y-axis [m]
    kappa = 1e-6  # thermal conductivity [W/m/K]
    u = 1  # velocity [m/s], relevant when there is advection

    # physical variable 1D specific
    T_left = 100.0  # float, boundary temperature at x=0 [C]
    T_right = 200.0  # flaot, boundary temperature at x=Lx [C]
    T_middle = 123.0  # flaot, initial temperature at t=0 [C]
    T_peak_start_x = 20  # flaot, start of the middle temperature on x-axis[m]
    T_peak_end_x = Lx-20  # flaot, end of the middle temperature on x-axis[m] (rest beside boundary will be zeros)

    # physical variables 2D specific
    T0 = 200  # initial temperature [C]
    T_max = 100  # maximum temperature [C]
    sigma = 1e4  # half-width of the peak [m]
    Q = 0  # heat flux [W/m2]
    xc = 2 / 3
    yc = 2 / 3

    # fem specific
    rho = 1  # density [kg/m3]
    cp = 1  # specific heat capacity [J/kg/K]

    # numerical variables
    nnx = 11  # number of grid points in x direction (use this for 1D domain as well)
    nny = 11  # number of grid points in y direction
    dt = 0.0 # time step [s], if 0 then it will be calculated automatically
    tol = 1e-2  # stop condition (error tolerance)
    max_iter = 500  # maximum number of iterations for the solvers
    CFL = 0.5  # Courant–Friedrichs–Lewy condition, making dt smaller (CFL < 1)
    output_filename = '1d_diffusion_implicit_jacobi'  # string, filename for the output file (e.g. 'plot_temp')

    # script below contain automatic computation, any changes should be made before to avoid hard-codding
    ##############################################################################################################

    ##### script preparation #####
    # setup of path to output files
    path_base = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))  # path to the base directory
    path_output_directory = os.path.join(path_base, 'output/')  # path to output folder
    if not os.path.exists(path_output_directory):  # create output folder if it doesn't exist
        os.makedirs(path_output_directory)
    output_filename_ext = output_filename + '.png'  # add extension to filename
    path_output = os.path.join(path_output_directory, output_filename_ext)  # path to output plot as png

    # creation of variables
    nnp = nnx * nny  # number of nodes total [-]
    hx = Lx / (nnx - 1)  # grid spacing on x-axis [m]
    hy = Ly / (nny - 1)  # grid spacing on y-axis [m]
    x = np.linspace(start_x, Lx + start_x, nnx)  # x-axis (grid points)
    if method == 'all' or method == '1d_diffusion_explicit' or method == '1d_diffusion_implicit' or method == '1d_fem_diffusion':
        dt = hx ** 2 / (2 * kappa) * CFL if dt == 0 else dt  # time step [s] (in case it is not specified)
    if method == '1d_advection' or method == '1d_fem_advection':
        dt = hx / u * CFL if dt == 0 else dt  # time step [s] (in case it is not specified)
    if method == '2d_diffusion_explicit' or method == '2d_diffusion_implicit' or method == '2d_fem_diffusion':
        dt = np.min([hx ** 2, hy ** 2]) / (
                2 * kappa) * CFL if dt == 0 else dt  # time step [s] (in case it is not specified)

    # initial temperature profile array (1D methods)
    T_init = np.zeros(nnx)  # initial temperature
    T_init[round(T_peak_start_x / hx):round(T_peak_end_x / hx)] = T_middle  # middle temperature
    T_init[0:round(T_peak_start_x / hx)] = T_left  # left boundary / left temperature
    T_init[round(T_peak_end_x / hx):-1] = T_right  # right boundary / right temperature
    T_init[-1] = T_right  # right boundary / right temperature

    ###### 1D explicit diffusion ######
    if method == '1d_diffusion_explicit' or method == 'all_1d_diffusion':
        print('Solving diffusion equation in 1d with explicit method...\n')
        print('--------------------------------------------------------')

        # computing solution using explicit method
        start = timeit.default_timer()  # start timer
        method_explicit_1d, T_explicit_1d, iter_count_explicit_1d, err_explicit_1d = methods.fde_explicit(T_init, kappa,
                                                                                                          nnx, hx, dt,
                                                                                                          tol, max_iter)
        stop = timeit.default_timer()  # stop timer
        time_explicit_1d = stop - start  # time needed to compute the solution
        operations.results_print(method_explicit_1d, err_explicit_1d, iter_count_explicit_1d, max_iter,
                                 time_explicit_1d)

    ###### 1D implicit diffusion ######
    if method == '1d_diffusion_implicit' or method == 'all_1d_diffusion':
        print('Solving diffusion equation in 1d with implicit method...\n')
        print('--------------------------------------------------------')

        if solver == 'all_1d_diffusion' or solver == '1d_scipy':
            # compute solution using scipy method
            start = timeit.default_timer()  # start timer
            method_scipy, T_scipy, iter_count_scipy, err_scipy = methods.fde_implicit_scipy_1d(T_init, kappa, dt, hx,
                                                                                               tol, max_iter)
            stop = timeit.default_timer()  # stop timer
            time_scipy = stop - start  # time needed to compute the solution
            operations.results_print(method_scipy, err_scipy, iter_count_scipy, max_iter, time_scipy)
        if solver == 'all_1d_diffusion' or solver == '1d_jacobi':
            # compute solution using jacobi method
            start = timeit.default_timer()  # start timer
            method_jacobi, T_jacobi, iter_count_jacobi, err_jacobi = methods.fde_implicit_jacobi_1d(T_init, kappa, dt,
                                                                                                    hx, tol,
                                                                                                    max_iter)
            stop = timeit.default_timer()  # stop timer
            time_jacobi = stop - start  # time needed to compute the solution
            operations.results_print(method_jacobi, err_jacobi, iter_count_jacobi, max_iter, time_jacobi)
        if solver == 'all_1d_diffusion' or solver == '1d_gauss_seidel':
            # compute solution using gauss-seidel method
            start = timeit.default_timer()  # start timer
            method_gauss, T_gauss, iter_count_gauss, err_gauss = methods.fde_implicit_gauss_1d(T_init, kappa, dt, hx,
                                                                                               tol, max_iter)
            stop = timeit.default_timer()  # stop timer
            time_gauss = stop - start  # time needed to compute the solution
            operations.results_print(method_gauss, err_gauss, iter_count_gauss, max_iter, time_gauss)
        if solver == 'all_1d_diffusion' or solver == '1d_crank-nicolson':
            # compute solution using crank-nicolson method
            start = timeit.default_timer()  # start timer
            method_crank, T_crank, iter_count_crank, err_crank = methods.fde_implicit_crank_1d(T_init, kappa, dt, hx,
                                                                                               tol, max_iter)
            stop = timeit.default_timer()  # stop timer
            time_crank = stop - start  # time needed to compute the solution
            operations.results_print(method_crank, err_crank, iter_count_crank, max_iter, time_crank)

    ###### 1D implicit advection ######
    elif method == '1d_advection':
        print('Solving advection equation in 1d...\n')
        print('--------------------------------------------------------')

        if solver == '1d_FTCS' or solver == 'all_1d_advection':
            # compute solution using FTCS method
            start = timeit.default_timer()  # start timer
            method_ftcs, T_ftcs, iter_count_ftcs, err_ftcs = methods.fde_advection_ftcs_1d(T_init, u, dt, hx, tol,
                                                                                           max_iter)
            stop = timeit.default_timer()  # stop timer
            time_ftcs = stop - start  # time needed to compute the solution
            operations.results_print(method_ftcs, err_ftcs, iter_count_ftcs, max_iter, time_ftcs)

        if solver == '1d_Lax-Friedrichs' or solver == 'all_1d_advection':
            # compute solution using Lax-Friedrichs method
            start = timeit.default_timer()  # start timer
            method_lax, T_lax, iter_count_lax, err_lax = methods.fde_advection_lax_1d(T_init, u, dt, hx, tol, max_iter)
            stop = timeit.default_timer()  # stop timer
            time_lax = stop - start  # time needed to compute the solution
            operations.results_print(method_lax, err_lax, iter_count_lax, max_iter, time_lax)

    ###### 2D ######
    if method == '2d_diffusion_explicit' or method == 'all_2d_diffusion':
        print('Solving diffusion equation in 2d with explicit method...\n')
        print('--------------------------------------------------------')
        # compute coordinates of the grid (aka relation between k index and coordinate x and y)
        print('WARNING, high CFL value may cause instability !!! \n\n') if CFL > 0.5 else None
        xcoords = np.zeros(nnp)  # initialize x coordinates
        ycoords = np.zeros(nnp)  # initialize y coordinates
        for j in range(0, nny):
            for i in range(0, nnx):
                xcoords[j * nnx + i] = i * hx + start_x
                ycoords[j * nnx + i] = j * hy + start_y

        # compute initial profile using analytical solution
        method_name_ana_2d, T_init, iter_count_ana_2d = methods.fde_diffusion_ana_init_2d(T_max, T0, kappa, sigma,
                                                                                          xcoords, ycoords, nnx, nny,
                                                                                          max_iter)
        # compute solution using explicit method
        start = timeit.default_timer()  # start timer
        method_name_explicit_2d, T_explicit_2d, iter_count_explicit_2d, err_explicit_2d = methods.fde_diffusion_explicit_2d(
            T_init, hx, hy, kappa, Q, nnx, nny, dt, tol, max_iter)
        stop = timeit.default_timer()  # stop timer
        time_explicit_2d = stop - start  # time needed to compute the solution
        operations.results_print(method_name_explicit_2d, err_explicit_2d, iter_count_explicit_2d, max_iter,
                                 time_explicit_2d)

    if method == '2d_diffusion_implicit' or method == 'all_2d_diffusion':
        print('Solving diffusion equation in 2d with implicit method...\n')
        print('--------------------------------------------------------')
        # compute coordinates of the grid (aka relation between k index and coordinate x and y)
        print('WARNING, high CFL value may cause instability !!! \n\n') if CFL > 0.5 else None
        xcoords = np.zeros(nnp)  # initialize x coordinates
        ycoords = np.zeros(nnp)  # initialize y coordinates
        for j in range(0, nny):
            for i in range(0, nnx):
                xcoords[j * nnx + i] = i * hx + start_x
                ycoords[j * nnx + i] = j * hy + start_y

        # compute initial profile using analytical solution
        method_name_ana_2d, T_init, iter_count_ana_2d = methods.fde_diffusion_ana_init_2d(T_max, T_init, kappa, sigma,
                                                                                          xcoords, ycoords, nnx, nny,
                                                                                          max_iter)
        # compute solution using explicit method
        start = timeit.default_timer()  # start timer
        method_name_implicit_2d, T_implicit_2d, iter_count_implicit_2d, err_implicit_2d = methods.fde_diffusion_implicit_2d(
            T_init, kappa, dt, hx, hy, nnx, nny, tol, max_iter)
        stop = timeit.default_timer()  # stop timer
        time_implicit_2d = stop - start  # time needed to compute the solution
        operations.results_print(method_name_implicit_2d, err_implicit_2d, iter_count_implicit_2d, max_iter,
                                 time_implicit_2d)

    if method == '2d_advection':
        print('Solving advection equation in 2d...\n')
        print('--------------------------------------------------------')
        # compute coordinates of the grid (aka relation between k index and coordinate x and y)
        print('WARNING, high CFL value may cause instability !!! \n\n') if CFL > 0.5 else None
        # compute initial velocity profile
        xcoords = np.zeros(nnp)  # initialize x coordinates
        ycoords = np.zeros(nnp)  # initialize y coordinates
        T_init = np.zeros(nnp)  # initialize initial temperature profile
        u = np.zeros(nnp)  # initialize velocity profile
        v = np.zeros(nnp)  # initialize velocity profile
        for j in range(0, nny):
            for i in range(0, nnx):
                k = j * nnx + i  # index of current grid point for when array is flattened (i and j together)
                # Grid
                xcoords[k] = i * hx
                ycoords[k] = j * hy
                # Initial temperature
                # if j==14 and i==20:
                if (xcoords[k] - xc) ** 2 + (ycoords[k] - yc) ** 2 <= sigma ** 2:
                    T_init[k] = 1 / 4 * (1 + np.cos(np.pi * (xcoords[k] - xc) / sigma)) * (
                            1 + np.cos(np.pi * (ycoords[k] - yc) / sigma))
                else:
                    T_init[k] = 0
                u[k] = -ycoords[k] + Ly / 2
                v[k] = xcoords[k] - Lx / 2

        # compute solution
        start = timeit.default_timer()  # start timer
        method_name_advection_2d, T_advection_2d, iter_count_advection_2d, err_advection_2d = methods.fde_advection_2d(
            T_init, u, v, hx, hy, nnx, nny, dt, max_iter)
        stop = timeit.default_timer()  # stop timer
        time_advection_2d = stop - start  # time needed to compute the solution
        operations.results_print(method_name_advection_2d, err_advection_2d, iter_count_advection_2d, max_iter,
                                 time_advection_2d)

    ###### FEM ######
    if method == '1d_fem_diffusion':
        # compute analitical solution
        T_ana = T_left + (T_right - T_left) / Lx * x
        # compute solution
        start = timeit.default_timer()  # start timer
        method_name_fem_diffusion_1d, T_fem_diffusion_1d, iter_count_fem_diffusion_1d, err_fem_diffusion_1d = methods.fem_diffusion_1d(
            T_init, kappa, rho, cp, hx, nnx, dt, tol, max_iter, T_ana)
        stop = timeit.default_timer()  # stop timer
        time_fem_diffusion_1d = stop - start  # time needed to compute the solution
        operations.results_print(method_name_fem_diffusion_1d, err_fem_diffusion_1d, iter_count_fem_diffusion_1d,
                                 max_iter,
                                 time_fem_diffusion_1d)

    if method == '1d_fem_advection':
        # compute solution
        err_fem_advection_1d = 0
        start = timeit.default_timer()  # start timer
        method_name_fem_advection_1d, T_fem_advection_1d, iter_count_fem_advection_1d = methods.fem_advection_1d(
            T_init, u, kappa, rho, cp, hx, nnx, dt, max_iter, solver)
        stop = timeit.default_timer()  # stop timer
        time_fem_advection_1d = stop - start  # time needed to compute the solution
        operations.results_print(method_name_fem_advection_1d, err_fem_advection_1d, iter_count_fem_advection_1d,
                                 max_iter, time_fem_advection_1d)

    if method == '2d_fem_diffusion':
        # TODO: work in progress
        pass

    ###### Plot ######
    if method == '1d_diffusion_explicit' or method == 'all_1d_diffusion':
        operations.results_plot_1d(path_output, x, T_init, T_explicit_1d, method_explicit_1d)
    if method == '1d_diffusion_implicit' and solver == '1d_scipy':
        operations.results_plot_1d(path_output, x, T_init, T_scipy, method_scipy)
    if method == '1d_diffusion_implicit' and solver == '1d_jacobi':
        operations.results_plot_1d(path_output, x, T_init, T_jacobi, method_jacobi)
    if method == '1d_diffusion_implicit' and solver == '1d_gauss_seidel':
        operations.results_plot_1d(path_output, x, T_init, T_gauss, method_gauss)
    if method == '1d_diffusion_implicit' and solver == '1d_crank-nicolson':
        operations.results_plot_1d(path_output, x, T_init, T_crank, method_crank)
    if method == '1d_advection' and solver == 'all_1d_advection':
        operations.results_plot_1d(path_output, x, T_init, T_ftcs, method_ftcs, T_lax,
                                   method_lax)  # plot initial and final temperature profiles
    if method == '1d_advection' and solver == 'Lax-Friedrichs':
        operations.results_plot_1d(path_output, x, T_init, T_lax, method_lax)
    if method == '2d_diffusion_explicit' or method == 'all_2d_diffusion':
        operations.results_plot_2d(path_output, T_init, xcoords, ycoords, nnx, nny, title='Initial temperature profile')
        operations.results_plot_2d(path_output, T_explicit_2d, xcoords, ycoords, nnx, nny,
                                   title='Final temperature profile')
    if method == '2d_diffusion_implicit' or method == 'all_2d_diffusion':
        operations.results_plot_2d(path_output, T_init, xcoords, ycoords, nnx, nny, title='Initial temperature profile')
        operations.results_plot_2d(path_output, T_implicit_2d, xcoords, ycoords, nnx, nny,
                                   title='Final temperature profile')
    if method == '2d_advection':
        operations.results_plot_2d(path_output, T_init, xcoords, ycoords, nnx, nny, title='Initial temperature profile')
        operations.results_plot_2d(path_output, T_advection_2d, xcoords, ycoords, nnx, nny,
                                   title='Final temperature profile')
    if method == '1d_fem_diffusion':
        operations.results_plot_1d(path_output, x, T_init, T_fem_diffusion_1d, method_name_fem_diffusion_1d)
    if method == '1d_fem_advection':
        operations.results_plot_1d(path_output, x, T_init, T_fem_advection_1d, method_name_fem_advection_1d)


if __name__ == "__main__":
    print('Starting main script')
    print('========================\n')
    main_function()
    print('\n========================')
    print('Done!')
    sys.exit(0)
