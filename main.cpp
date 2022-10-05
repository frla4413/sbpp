//============================================================
// Name        : main.cpp
// Date        : May 2022
// Author      : Fredrik Laurén
// Description : main-file
//============================================================

#include <iostream>
#include <valarray>
#include "array.hpp"
#include "mbarray.hpp"
#include "mesh.hpp"
#include "sbp21.hpp"
#include "sbp42.hpp"
#include "mbgrid.hpp"
#include "mbsbp.hpp"

double F(double x, double y) {
  return std::sin(2*M_PI*x)*std::cos(2*M_PI*y);
}

double dF(double x, double y) {
  return 2*M_PI*std::cos(2*M_PI*x)*std::cos(2*M_PI*y);
}

int main() {

  int N = 81;

  auto blocks{Annulus(N, 0.1, 1)};
  MbSbp sbp {blocks, 4};
  auto f = sbp.Evaluate(F);

  auto error = sbp.Dx(f) - sbp.Evaluate(dF);
  std::cout << sqrt(sbp.Integrate(error*error)) << "\n";
  ExportToTec(sbp.grid(), error, "test");

  return 0;
}
