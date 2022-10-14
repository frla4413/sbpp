This is a set of c++ functions that can be used to solve a
intiial-boundary value problem.

The numerical scheme is based uses finite-difference
operators on summation-by-parts form.

Requirements:
  o gcc
  o cmake
  o Gtest (for unit testing)
  o MKL (for GMRES used in implicit time integration)

To build:
  o mkdir Debug
  o cmake ..
  o make

Schematic overwiew of the solution process:
  i) Specify the geometry - Ex. blocks {Annulus()}

  ii) Form an advection object
     (boundary conditions are already determined) -
     advection {blocks, order}

  iii) Set initial condition - init {init_func()}

  iv) Integrate in time -
            sol{Integration([t0 t1], init, dt, rhs_fun)

See advection/main*.cpp for more details.
