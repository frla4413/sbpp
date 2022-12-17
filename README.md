A collection of c++ functions that can be used to solve intiial-boundary value problems.

The numerical scheme uses difference operators on summation-by-parts form. Details are found in my PhD-thesis
http://www.diva-portal.org/smash/record.jsf?pid=diva2%3A1614643&dswid=1629 and references therein.

The script ins/ins_main.cpp shows the typical solution process for the 2D incompressible Navier-Stokes equations.

The lid-driven cavity is a popular benchmark problem. The lower, left and right boundaries are solid walls.
At the upper boundary, (u,v) = (0,1).
![](https://github.com/frla4413/sbpp/blob/main/images/cavity.png)

Another example is shown below. The west side is the inflow boundary, while the top-right one is set to an outflow boundary. The remaining side are walls.
![](https://github.com/frla4413/sbpp/blob/main/gifs/channel.gif)

A final example is shown below, where the flow is coming from the left to the right. The obstacle interferes with the fluid, leading to vertices.
![](https://github.com/frla4413/sbpp/blob/main/gifs/vorticity.gif)

Requirements:
* compiler c++17
* cmake
* Gtest (for unit testing)
* MKL (for GMRES used in implicit time integration)
* Progress bars (https://github.com/mraggi/tqdm-cpp)

Build:
1. mkdir Release 
2. cmake ..
3. make

A simpler setup is found in the advection folder. See advection/main_explicit.cpp for a setup of the advection equation on an annulus, shown below.

![](https://github.com/frla4413/sbpp/blob/main/gifs/annulus.gif)
