//============================================================================
// Name        : main_explicit.cpp
// Author      : frela05
// Description : main-file for advection equation u_t + a1 ux + a2 uy = F
//               This test case includes a manufactured solution.
//               Explicit time integration is used.
//               The solution of the finest grid is saved to "adv/sol".
//               This file must be created before execution.
//============================================================================
#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <math.h>
#include <valarray>
#include "array.hpp"
#include "mbarray.hpp"
#include "mesh.hpp"
#include "mbgrid.hpp"
#include "mbsbp.hpp"
#include "advection.hpp"
#include "integrate.hpp"
#include "basics.hpp"

int main() {
  double t0 = 0;
  double t1 = 1;

  std::vector<double> tspan = {t0,t1};

  std::pair<double,double> a{1,1};

  std::valarray<int> N_vec = {21,41};

  std::valarray<double> err_vec(N_vec.size());
  int order = 2;

  int err_vec_pos = 0;

  for (auto& N : N_vec) {
    auto block{CartesianGrid(N, N)};
    std::vector<Block> blocks {block};
    block.x += 1;
    blocks.push_back(block);
    Advection advec{blocks, order, a};
    MbArray initial{advec.AnalyticVec(0.0)};

    double h  = 1.0/(N-1);
    double dt = 0.1*h;

    auto rhs_fun = [&] (double t, const MbArray& y) {
      bool mms;
      return advec.Rhs(t,y, mms = true);
    };

    auto sol = ExplicitIntegration(tspan, initial, dt, rhs_fun);

    err_vec[err_vec_pos] = advec.ComputeError(sol.t,sol.y);
    std::cout << "\nError: " << err_vec[err_vec_pos] << "\n";
    ++err_vec_pos;
  }
  std::cout << "done\n";
  PrintConvergence(N_vec, err_vec);
  return 0;
}
