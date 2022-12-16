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

//====================== User input ===============================

  // N is a parameter of the grid.
  // The number of grid points depend on N.
  int N = 641;

  // Set Sbp order
  int order = 4;

  // Viscosity parameter
  double mu = 0.01;

  auto setup {ChannelOne(N)};

  Ins ins {setup.blocks, order, setup.bd_types, mu};

  InsState init = InitialData(ins);
  std::cout << "Total size " << init.u.GetTotalSize() << "\n";
  auto lhs = ins.SpatialOperator(0, init);

  return 0;
}
