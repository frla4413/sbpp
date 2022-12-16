This is a set of c++ functions that can be used to solve an intiial-boundary value problem.

The numerical scheme uses finite-difference operators on summation-by-parts form. Details are found in my PhD-thesis
http://www.diva-portal.org/smash/record.jsf?pid=diva2%3A1614643&dswid=1629 and references therin.

Below is the popular benchmark: The lid-driven cavity. The lower, loeft and rigth boundaries are solid walls. At the upper boundary, we have imposed (u,v) = (0,1). The diffusion coefficient is 0.001.
![](https://github.com/frla4413/sbpp/blob/main/images/cavity.png)

Another example is shown below. The most west boundary is an inflow, while the top-right one is set to be the outflow. The remaining side are walls. The diffusion coefficient is set to 0.001.

A final example is whon below, where the flow is coming from the left to the right. The obstacle in its way is forcing the fluid to go around the sqare, which produces the vertices.
![](https://github.com/frla4413/sbpp/blob/main/gifs/vorticity_video.gif)

The scripts ins/ins_main.cpp shows the typical solution process for the 2D incompressible Navier-Stokes equations.

Requirements:
* compiler c++17
* cmake
* Gtest (for unit testing)
* MKL (for GMRES used in implicit time integration)
* Progress bars (https://github.com/mraggi/tqdm-cpp)

To build:
1. mkdir Debug
2. cmake ..
3. make

A more simple setup is found in the advection folder. See advection/main_explicit.cpp for a setup of the advection equation on an annulus. The solution is shown below.

![](https://github.com/frla4413/sbpp/blob/main/gifs/annulus.gif)
