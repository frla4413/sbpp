//================================================================
// Name        : ins_main.cpp
// Author      : frela05
// Description : main-file for the INS equations.
//               Implicit time integration is used.
//               This folder "sol" must be created before execution.
//================================================================
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

InsState InitialData(Ins& ins) {
  auto initu = [] (double x, double y) {
     return 1.0;
  };
  auto initv = [] (double x, double y) { return 0; };
  auto initp = [] (double x, double y) { return 0; };
  InsState init;
  init.u = ins.Evaluate(initu);
  init.v = ins.Evaluate(initv);
  init.p = ins.Evaluate(initp);
  return init;
}

int main() {

  int Nx = 41;
  int Ny = 81;

  auto block1 {CartesianGrid(Nx, Ny)};
  auto block2 {CartesianGrid(21, 41)};
  block2.x = block2.x + 1;
  std::vector<Block> blocks {block1, block2};
  std::vector<BdType> bd_types(6,BdType::Wall);
  bd_types[2] = BdType::Inflow;
  bd_types[4] = BdType::Outflow;

  int order = 4;
  double mu = 1e-2;
  Ins ins {blocks, order, bd_types, mu};
  InsState init {InitialData(ins)};

  std::string name_base = "plate/sol";
  auto sol {SolveSteadyState(init, ins, name_base)};

  return 0;
}
