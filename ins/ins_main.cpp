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

InsState InitialData(Ins& ins) {
  auto initu = [] (double x, double y) {
     return 1.0 + 0.1*exp(-pow((x - 0.5)/0.1,2))*
                      exp(-pow((y - 0.5)/0.1,2));
  };
  auto initv = [] (double x, double y) { return 0; };
  auto initp = [] (double x, double y) { return 0; };
  InsState init;
  init.u = ins.Evaluate(initu);
  double normalize = abs(init.u.ToValarray()).max();
  init.u = init.u/normalize;
  init.v = ins.Evaluate(initv);
  init.p = ins.Evaluate(initp);
  return init;
}

// TO DO:
// Clean up jacobians
// MKL-operations for Array2D? --
// try on advection and compare operations
int main() {

  double t0 = 0, t1 = 10, dt = 0.05;

  Tspan tspan {t0,t1};
  int N = 51;
  auto block {CartesianGrid(N, N)};
  std::vector<Block> blocks {block};
  std::vector<BdType> bd_types(4,BdType::Wall);

  int order = 4;

  double mu = 0.1;
  Ins ins {blocks, order, bd_types, mu};
  InsState init = InitialData(ins);

  auto name_base = "ins_sol/sol";

  InsSolution sol = ImplicitTimeIntegraction(tspan, init, dt, ins,
                                             name_base);
  return 0;
}
