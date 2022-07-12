#! /usr/bin/env python3
'''
Comparing numerical integration: middle point rule VS trapezoidal rule

We will try with different function, describe as follow:
f1 = x
f2 = cos(x)

Remarks:
Both rules give a good approximation (less than 1\% error). However, the middle point rule converge faster and is
thefore more accurate ($n>80$ for the trapeze give the same accuracy as $n=10$ for the middle point rule).
'''

# import
import sys
import timeit

import matplotlib.pyplot as plt
import numpy as np


def main_function():
    # input Variables
    a = 0  # lower limit
    b = np.pi / 2  # upper limit
    n_start = 10 # starting number of subintervals
    n_max = 100  # maximum number of subintervals
    n_interval = 10 # interval between two consecutive numbers of subintervals (for the plot)

    # initialize variables for the higher loop
    dn = round((n_max - n_start) / n_interval + 1) # interval between two consecutive numbers of subintervals
    n_axis = np.linspace(n_start, n_max, dn) # numbers of subintervals tested for the higher loop
    mpr_error_axis = np.zeros(dn) # error computed for each number of subintervals
    trap_error_axis = np.zeros(dn) # error computed for each number of subintervals
    iteration_number = 0 # initialise iteration number for main loop

    # loop over various n values to plot error in function of number of subintervals
    for n in range(n_start, n_max, n_interval):
        # computing variables
        nnx = n + 1  # number of points in x-axis (number of interval + 1)
        h = (b - a) / (nnx - 1)  # compute subinterval size
        x = np.linspace(a, b, nnx)  # compute x-axis points

        # initialise empty arrays for storing solution
        mpr_f2_array = np.zeros(nnx)  # middle point rule, array initialization for f(x) = cos(x)
        trap_f2_array = np.zeros(nnx) # trapezoidal rule, array initialization for f(x) = cos(x)

        # compute analytical solution
        ana_f2 = np.sin(b) - np.sin(a) # analytical solution for f(x) = cos(x) (int between 0 and 2pi)

        # compute solution using middle point rule (Exercise 1) and trapezoidal rule (Exercise 2)
        for i in range(1, nnx):
            mpr_f2_array[i] = h * np.cos( (x[i-1]+x[i])/2 )
            trap_f2_array[i] = h * (np.cos(x[i-1]) + np.cos(x[i])) / 2
        mpr_f2 = np.sum(mpr_f2_array) # computing approximate solution of the integral by summing all the subintervals
        trap_f2 = np.sum(trap_f2_array)
        # error_axis[iteration_number] = np.abs(trap_f2 - ana_f2) # computing abs error of the solution
        mpr_error_axis[iteration_number] = np.abs(mpr_f2 - ana_f2)
        trap_error_axis[iteration_number] = np.abs(trap_f2 - ana_f2)

        # print results
        print("\nMiddle point rule error = %.6f" % mpr_error_axis[iteration_number])
        print("\nTrapezoidal rule error = %.6f" % trap_error_axis[iteration_number])
        print('--')

        iteration_number += 1  # increment iteration number

    # plotting
    plt.figure(figsize=(6, 5))
    plt.title('Absolute error as a function of number of subintervals')
    plt.plot(n_axis, mpr_error_axis, color='blue', label='f(x) = cos(x) middle point rule')
    plt.plot(n_axis, trap_error_axis, color='red', label='f(x) = cos(x) trapezoidal rule')
    plt.xlabel('Subinterval size (h)')
    plt.ylabel('Absolute error of approximated solution')
    plt.legend()
    plt.show()



if __name__ == "__main__":
    print('Starting Numerical Integrations')
    start = timeit.default_timer()
    main_function()
    stop = timeit.default_timer()
    print('Done!')
    print('Time: ', round(stop - start, 4), 's')
    sys.exit(0)