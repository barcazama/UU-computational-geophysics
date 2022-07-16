# README for Computational Geophysics
## Introduction
This project is a collection of scripts that regroup the deliverable for the course Computational Geophysics at Utrecht University GEO4-1427, directed by [Asst. Prof. Cedric Thieulot](https://cedricthieulot.net/). 

Normally, in any practical situation, we would not rewrite software from scratch, but I believe there is no better way to learn how something works than to completely rebuild it,. Therefore, this project will aim to build all the basic elements of a geomechanical simulation software as well as compare different schemes, solvers and methods to achieve this goal. In the same vein, the focus while coding these scripts will be readability, so the code will not be optimized for speed as it is often the case for python scripts.

## Structure
- `bin/`: contains the executable files
  - `main.py`: answer of the exercises from [Fieldstone](https://cedricthieulot.net/manual.pdf). Parameters can be set in the beginning to fit different way to solve either diffusion or advection equations.
  - `numerical-integration-rules-ex#`: comparison of two methods of numerical approximation of a solution of integral equations: the midpoint rule or the trapezoidal rule. 
  - `numerical-integration-quadrature-ex#`: experimentation of the quadrature rule which is another way to solve integral equations in an approximate way. It consists in replacing the equation by a sum of integrands evaluated at several quadrature points multiplied by quadrature weights.
- `modules/`: contains the functions and classes
  - `methods.py`: highest level of the tree, contain scheme and methodology.
  - `solvers.py`: mid level of the tree, contain the solvers.
  - `operations.py`: lowest level of the tree, contain simple calculations or function used to present results.
- `output/`: contains the output files
- `test/`: contains the test files
- `environment.yml`: file to set environment with conda

All exercises for FDM and FEM are executable from the `bin/main.py` script, with instructions for each exercise present in [Fieldstone](https://cedricthieulot.net/manual.pdf). In addition, three scripts contain the comparison of different numerical integration methods, located in the `bin/` directory. The `modules/` directory contains all the functions that are not called directly and are organized into three different scripts: `methods.py`, `solvers.py` and `operations.py`. `methods. py` is the highest level of the tree, so it is able to call the other scripts in the same folder, and contains the complete schema and methodology. `solvers.py` is one level lower than `methods.py` and contains the different solvers, this script is able to call the other scripts `operations.py`. This last script is at the lowest level of the tree and will not call any other function in this project, it contains simple calculations or formatting functions such as result plotting.
