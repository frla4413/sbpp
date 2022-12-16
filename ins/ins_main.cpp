//==================================================================
// Name        : ins_main.cpp
// Author      : frela05
// Description : main-file for the INS equations.
//               Implicit time integration is used.
//               This folder "sol" must be created before execution.
//==================================================================
#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <cmath>
#include <valarray>
#include "array.hpp"
#include "mbarray.hpp"
#include "mesh.hpp"
#include "mbgrid.hpp"
#include "mbsbp.hpp"
#include "integrate.hpp"
#include "basics.hpp"
#include "ins.hpp"
#include "integrate_ins.hpp"
#include "ins_setup.hpp"

int main() {

//=============== User input ======================================

  // Set time limits and time step
  double t_start = 0, t_end = 10, dt = 0.025;

  // N is a parameter of the grid.
  // The number of grid points depend on N.
  int N = 121;

  // Set Sbp order
  int order = 4;

  // Viscosity parameter
  double mu = 1e-3;

  // Folder to save solution
  auto solution_folder = "ins_sol/sol";

//=================================================================

  auto setup {Channel(N)};

  Tspan tspan {t_start, t_end};
  Ins ins {setup.blocks, order, setup.bd_types, mu};
  InsState init = InitialData(ins);

  InsSolution sol = ImplicitTimeIntegraction(tspan, init, dt, ins,
                                             solution_folder);
  return 0;
}
