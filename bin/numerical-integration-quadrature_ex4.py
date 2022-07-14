#! /usr/bin/env python3
'''
Comparing numerical integration: number of quadratic points for gaussian integration. Relate to Exercise 4 of the
numerical integration (see presentation). We will try for 2, 3 and 4 quadratic points approximation and compare them
with the analytical solution calculated by scipy. The following function is tested:
f1(x) = x² + 4y
Results / Comments:
An anomaly is observed for nq = 2, there is no error, and we find the exact solution. Otherwise, the error is
relatively high with butter results for nq = 4 than nq = 3.
'''

# import
import sys
import timeit

import matplotlib.pyplot as plt
import numpy as np
from scipy import integrate


def main_function():
    # parameters
    for nq in [2,3,4]: # number of quadrature points
        x_min, x_max = 7, 10 # min, max value of x
        y_min, y_max = 10, 14 # min, max value of y

        # compute the analytical solution using scipy library
        f1_ana, f1_scipy_error = integrate.dblquad(f1, y_min, y_max, x_min, x_max)

        # computation of quadratic parameters
        if nq == 2:
            xq = [-np.sqrt(1 / 3), np.sqrt(1 / 3)]
            yq = [-np.sqrt(1 / 3), np.sqrt(1 / 3)]
            wq = [1, 1]
        elif nq == 3:
            xq = [-np.sqrt(3 / 5), 0, np.sqrt(3 / 5)]
            yq = [-np.sqrt(3 / 5), 0, np.sqrt(3 / 5)]
            wq = [5 / 9, 8 / 9, 5 / 9]
        elif nq == 4:
            xq = [-np.sqrt(3 / 7 + 2 / 7 * np.sqrt(6 / 5)), -np.sqrt(3 / 7 - 2 / 7 * np.sqrt(6 / 5)),
                  np.sqrt(3 / 7 - 2 / 7 * np.sqrt(6 / 5)), np.sqrt(3 / 7 + 2 / 7 * np.sqrt(6 / 5))]
            yq = [-np.sqrt(3 / 7 + 2 / 7 * np.sqrt(6 / 5)), -np.sqrt(3 / 7 - 2 / 7 * np.sqrt(6 / 5)),
                  np.sqrt(3 / 7 - 2 / 7 * np.sqrt(6 / 5)), np.sqrt(3 / 7 + 2 / 7 * np.sqrt(6 / 5))]
            wq = [(18 - np.sqrt(30)) / 36, (18 + np.sqrt(30)) / 36, (18 + np.sqrt(30)) / 36,
                  (18 - np.sqrt(30)) / 36]
        else:
            print('Error: number of quadratic points not supported')
            sys.exit()

        # initialize variables
        f1_sol = 0

        # computation of the quadratic solution (for nq points)
        for i in range(nq):
            for j in range(nq):
                x = xq[i] * (x_max - x_min)/2 + (x_max + x_min)/2
                y = yq[j] * (y_max - y_min)/2 + (y_max + y_min)/2
                f1_sol += f1(x, y) * wq[0] * wq[1]
        f1_sol *= (x_max - x_min) * (y_max - y_min) / 4

        # compute error
        f1_error = np.abs(f1_sol - f1_ana)

        # print the results
        print('Using %d quadratic points:' % nq)
        print("\nError of the approximated solution for f1(x) = %.6f" % f1_error)
        print('--')

# function to integrate
def f1(x, y):
    return x**2 + 4*y


if __name__ == "__main__":
    print('Starting Numerical Integrations\n')
    print('-----------------------------------------------------')
    start = timeit.default_timer()
    main_function()
    stop = timeit.default_timer()
    print('Done!')
    print('Time: ', round(stop - start, 4), 's')
    sys.exit(0)