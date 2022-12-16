//==================================================================
// Name        : ins_test.cpp
// Author      : frela05
// Description : Test openmp
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

#include <chrono>
int main() {

using namespace std::chrono;
//=============== User input ======================================

  // N is a parameter of the grid.
  // The number of grid points depend on N.
  int N = 2000;

  // Set Sbp order
  int order = 4;

  auto setup {Channel(N)};

  MbSbp sbp{setup.blocks, order};
  auto func = [](double x, double y) {

    return x;
  };
  auto f {sbp.Evaluate(func)};

  auto start = high_resolution_clock::now();
#pragma omp parallel for
  for (int i = 0; i < 50; ++i)
    auto df {sbp.Dy(f)};
  auto stop = high_resolution_clock::now();

  auto duration = duration_cast<milliseconds>(stop - start);
  std::cout << duration.count() << std::endl;

  return 0;
}
