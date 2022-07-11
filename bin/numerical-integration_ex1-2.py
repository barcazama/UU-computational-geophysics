#! /usr/bin/env python3
'''
Comparing numerical integration: middle point rule VS trapezoidal rule

While not in used here, a function for each rule has been written. They will be accessible from `modules/methods`

Remarks:
Both rules give a good approximation (less than 1\% error). However, the middle point rule converge faster and is
thefore more accurate ($n>80$ for the trapeze give the same accuracy as $n=10$ for the middle point rule).
'''

# import
import sys
import timeit

import matplotlib.pyplot as plt
import numpy as np


# TODO: add third function for bonus
def main_function():
    # Variables
    a = 0  # lower limit
    b = np.pi / 2  # upper limit
    n_max = 1000  # maximum number of subintervals
    modulo = 10  # control the rate at which the results are printed (higher number = less print)

    # initialize variables
    n = 10  # initialize subinterval variable
    error_array_mpr_f1 = np.empty(0)  # initialize error array for f(x) = x
    error_array_mpr_f2 = np.empty(0)  # initialize error array for f(x) = cos(x) (middle point rule)
    error_array_trap_f2 = np.empty(0)  # initialize error array for f(x) = cos(x) (trapeze rule)
    h_array = np.empty(0)  # initialize h array

    # compute analytical solution
    I_ana_f1 = b ** 2 / 2 - a ** 2 / 2
    I_ana_f2 = np.sin(b) - np.sin(a)  # analytical solution

    # iterations
    while n <= n_max:
        # initilize variables
        I_approx_mpr_f1 = 0
        I_approx_mpr_f2 = 0

        # compute subinterval variables
        h = (b - a) / (n-1)  # compute subinterval size
        x = np.linspace(a, b, n)  # compute subinterval points

        # compute approximated solution using middle point rule (mpr)
        for i in range(n - 1):
            I_approx_mpr_f1 += (x[i] + x[i + 1]) / 2
            I_approx_mpr_f2 += np.cos((x[i] + x[i + 1]) / 2)
        I_approx_mpr_f1 *= h
        I_approx_mpr_f2 *= h

        # compute approximated solution using trapezoidal rule (trap)
        # Note: for f(x)=x, the syntax doesn't change between trap or mpr
        I_approx_trap_f2 = np.cos(x[0]) + np.cos(x[-1])
        for i in range(1, n - 1):
            I_approx_trap_f2 += np.cos(x[i]) + np.cos(x[i + 1])
        I_approx_trap_f2 *= h / 2

        error_mpr_f1 = np.abs(I_approx_mpr_f1 - I_ana_f1)  # compute error
        error_mpr_f2 = np.abs(I_approx_mpr_f2 - I_ana_f2)
        error_trap_f2 = np.abs(I_approx_trap_f2 - I_ana_f2)
        error_array_mpr_f1 = np.append(error_array_mpr_f1, error_mpr_f1)
        error_array_mpr_f2 = np.append(error_array_mpr_f2, error_mpr_f2)
        error_array_trap_f2 = np.append(error_array_trap_f2, error_trap_f2)
        h_array = np.append(h_array, h)  # append h to array

        if n % modulo == 0:
            print('solution for n =', n, ' and h =', h)
            print('I_approx_f1 =', I_approx_mpr_f1, ' I_ana_f1 =', I_ana_f1, ' error_f1 =', error_mpr_f1)
            print('I_approx_mpr_f2 =', I_approx_mpr_f2, ' I_ana_f2 =', I_ana_f2, ' error_mpr_f2 =', error_mpr_f2)
            print('I_approx_trap_f2 =', I_approx_trap_f2, ' I_ana_f2 =', I_ana_f2, ' error_trap_f2 =', error_trap_f2)
            print('\n')
        n += 1  # increment subinterval variable

    # plotting
    plt.figure(figsize=(6, 5))
    plt.title('Absolute error using middle point rule')
    plt.plot(h_array, error_array_mpr_f1, color='black', label='f(x) = x')
    plt.plot(h_array, error_array_mpr_f2, color='blue', label='f(x) = cos(x) (middle point rule)')
    plt.plot(h_array, error_array_trap_f2, color='red', label='f(x) = cos(x) (trapezoidal rule)')
    plt.xlabel('Subinterval size (h)')
    plt.ylabel('Absolute error of approximated solution')
    plt.legend()
    plt.show()


def midpoint_rule(f, a, b, n):
    """
    Midpoint rule for numerical integration
    f: function to be integrated
    a: lower limit
    b: upper limit
    n: number of subintervals
    """
    h = (b - a) / (n-1)  # compute subinterval length
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
    """
    Trapeze rule for numerical integration
    f: function to be integrated
    a: lower limit
    b: upper limit
    n: number of subintervals
    """
    h = (b - a) / (n-1)  # compute subinterval length
    x = np.linspace(a, b, n)  # compute subinterval points
    I = f(x[0]) + f(x[-1])
    for i in range(1, n - 1):
        I += f(x[i]) + f(x[i + 1])
    I *= h / 2
    return I


if __name__ == "__main__":
    print('Starting Numerical Integration - Exercise 1 and 2...')
    start = timeit.default_timer()
    main_function()
    stop = timeit.default_timer()
    print('Done!')
    print('Time: ', stop - start)
    sys.exit(0)
