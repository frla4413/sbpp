This is a set of c++ functions that can be used to solve an intiial-boundary value problem.

The numerical scheme is based uses finite-difference operators on summation-by-parts form.

Requirements:
* gcc
* cmake
* Gtest (for unit testing)
* MKL (for GMRES used in implicit time integration)

To build:
1. mkdir Debug
2. cmake ..
3. make

Schematic overview of the solution process:
1. Specify the geometry - Ex. blocks = Annulus()
2. Form an advection object (boundary conditions are already determined) - advection = (blocks, order)
3. Set initial condition - init  = init_func()
4. Integrate in time - sol = Integration([t0 t1], init, dt, rhs_fun)

See advection/main**.cpp for more details.

Below is the solution of the advection equation (with a forcing function) posed on an annulus.
![](https://github.com/frla4413/sbpp/blob/main/annulus.gif)
