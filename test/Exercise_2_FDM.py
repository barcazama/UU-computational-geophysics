#! /usr/bin/env python3
"""
Exercise 2: Finite Difference Method - Explicit Method - heat equation
====================================================================
"""

# import
import os
import sys
import timeit

import matplotlib.pyplot as plt
import numpy as np


def fde_explicit():
    # Variables
    path_base = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    path_output = os.path.join(path_base, 'output/')
    if not os.path.exists(path_output):
        os.makedirs(path_output)
    ## Physical
    Lx = 1e3  # Length of the domain [m]
    Tleft = 100  # boundary temperature at x=0 [C]
    Tright = 200  # boundary temperature at x=Lx [C]
    Tinit = 123  # initial temperature at t=0 [C]
    kappa = 1e-6  # thermal conductivity [W/m/K]
    ## Numerical
    nnx = 10  # number of grid points in x direction
    tol = 0.1  # stop condition

    # Computing variables
    h = Lx / (nnx - 1)  # grid spacing
    dt = h ** 2 / (2 * kappa)  # time step
    x = np.arange(0, Lx, h)  # grid points
    Told = np.ones(nnx - 1)  # initial temperature
    Told[1:-1] = Tinit  # initial temperature
    Told[0] = Tleft  # left boundary
    Told[-1] = Tright  # right boundary
    Tnew = np.zeros(nnx - 1)  # new temperature (initialized to zero)
    Tnew[0] = Tleft  # left boundary
    Tnew[-1] = Tright  # right boundary
    T_initial = Told.copy()  # initial temperature profile

    # Solving the equation
    error = 1  # initialize error
    n = 1  # initialize iteration counter
    while error > tol:
        for i in range(1, nnx - 2):
            Tnew[i] = Told[i] + dt * kappa * (Told[i - 1] - 2 * Told[i] + Told[i + 1]) / (h ** 2)
        error = np.max(np.abs(Tnew - Told))  # compute error
        print('timestep #', n, ' min/max (boundary excluded) ', np.min(Tnew[1:-1]), np.max(Tnew[1:-1]), ' error = ',
              error)
        Told[:] = Tnew[:]  # update Told
        n += 1  # increment iteration counter

    # Plotting
    plt.figure(figsize=(6, 5))
    plt.title('Reaching steady state with explicit method')
    plt.plot(x, Tnew, color='blue', linewidth=2, label='T_profile last iteration')
    plt.plot(x, T_initial, color='black', linewidth=2, label='T_profile initial situation')
    plt.xlabel('x (m)')
    plt.ylabel('T (K)')
    plt.grid()
    plt.legend()
    path_output_temperature = os.path.join(path_output, 'FDE_implicit_Temperature.png')
    plt.savefig(path_output_temperature)
    plt.show()


if __name__ == "__main__":
    print('Starting FDE implicit solving...')
    start = timeit.default_timer()
    fde_explicit()
    stop = timeit.default_timer()
    print('Done!')
    print('Time: ', stop - start)
    sys.exit(0)
