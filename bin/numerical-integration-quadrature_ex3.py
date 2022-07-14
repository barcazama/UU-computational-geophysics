#! /usr/bin/env python3
'''
Comparing numerical integration: number of quadratic points for gaussian integration. Relate to Exercise 3 of the
numerical integration (see presentation). We will try for 2 and 5 quadratic points approximation and compare them with
the analytical solution calculated by scipy. The three following functions are tested:
f1(x) = sin(x*pi + pi/2)
f2(x) = sqrt(x+1)
f3(x) = x^4 - x^3
Results / Comments:
As expected, the results for a higher number of quadratic points are more accurate (less misfit) and is acceptable for
most application with 5 quadratic points. In contrary, using 2 quadratic points does not give a good approximation
(almost 0.5 mistfit for the f1).
'''

# import
import sys
import timeit

import matplotlib.pyplot as plt
import numpy as np
from scipy import integrate


def main_function():
    # parameters
    for nq in [2,5]: # number of quadrature points
        x_max = 1 # maximum value of x
        x_min = -1 # minimum value of x

        # compute the analytical solution using scipy library
        f1_ana, f1_scipy_error = integrate.quad(f1, x_min, x_max)
        f2_ana, f2_scipy_error = integrate.quad(f2, x_min, x_max)
        f3_ana, f3_scipy_error = integrate.quad(f3, x_min, x_max)

        # computation of quadratic parameters
        if nq == 2:
            xq = [-np.sqrt(1/3), np.sqrt(1/3)] # quadratic points
            wq = [1, 1] # quadratic weights
        elif nq == 5:
            xq = [-1 / 3 * np.sqrt(5 - 2 * np.sqrt(10 / 7)), -1 / 3 * np.sqrt(5 + 2 * np.sqrt(10 / 7)), 0,
                  1 / 3 * np.sqrt(5 - 2 * np.sqrt(10 / 7)), 1 / 3 * np.sqrt(5 + 2 * np.sqrt(10 / 7))] # quadratic points
            wq = [(322 + 13 * np.sqrt(70)) / 900, (322 - 13 * np.sqrt(70)) / 900, 128 / 225,
                  (322 + 13 * np.sqrt(70)) / 900, (322 - 13 * np.sqrt(70)) / 900] # quadratic weights
        else:
            print('Error: number of quadratic points not supported')
            sys.exit()

        # initialize variables
        f1_sol = 0
        f2_sol = 0
        f3_sol = 0

        # computation of the quadratic solution (for nq points)
        for i in range(nq):
            f1_sol += wq[i] * f1(xq[i])
            f2_sol += wq[i] * f2(xq[i])
            f3_sol += wq[i] * f3(xq[i])

        # compute error
        f1_error = np.abs(f1_sol - f1_ana)
        f2_error = np.abs(f2_sol - f2_ana)
        f3_error = np.abs(f3_sol - f3_ana)

        # print the results
        print('Using %d quadratic points:' % nq)
        print("\nError of the approximated solution for f1(x) = %.6f" % f1_error)
        print("\nError of the approximated solution for f2(x) = %.6f" % f2_error)
        print("\nError of the approximated solution for f3(x) = %.6f" % f3_error)
        print('--')

# function to integrate
def f1(x):
    return np.sin(x*np.pi + np.pi/2)
def f2(x):
    return np.sqrt(x+1)
def f3(x):
    return x**4 - x**3


if __name__ == "__main__":
    print('Starting Numerical Integrations\n')
    print('-----------------------------------------------------')
    start = timeit.default_timer()
    main_function()
    stop = timeit.default_timer()
    print('Done!')
    print('Time: ', round(stop - start, 4), 's')
    sys.exit(0)