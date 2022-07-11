#! /usr/bin/env python3
'''
Implementation of the gaussian quadrature method for numerical integration
'''

# import
import sys
import timeit

import matplotlib.pyplot as plt
import numpy as np

def main_function():
    # Variables
    a = 0  # lower limit
    b = np.pi / 2  # upper limit
    n_max = 1000  # maximum number of subintervals

    # initialize variables
    n = 10  # initialize subinterval variable
    x = np.linspace(a, b, n)  # compute subinterval points
    r = 2/(b-a)*(x-a)-1 # normalized x-array between -1 and 1
    dr = 2/(n-1) # normalized subinterval size



if __name__ == "__main__":
    print('Starting Numerical Integration - Exercise 1...')
    start = timeit.default_timer()
    main_function()
    stop = timeit.default_timer()
    print('Done!')
    print('Time: ', stop - start)
    sys.exit(0)
