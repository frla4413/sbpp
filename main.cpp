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
#include "interp.hpp"
#include "interp21.hpp"

//double F(double x, double y) {
//  return std::sin(2*M_PI*x)*std::cos(2*M_PI*y);
//}
//
//double dF(double x, double y) {
//  return 2*M_PI*std::cos(2*M_PI*x)*std::cos(2*M_PI*y);
//}
//
//int main() {
//
//  int N = 81;
//
//  auto blocks{Annulus(N, 0.1, 1)};
//  MbSbp sbp {blocks, 4};
//  auto f = sbp.Evaluate(F);
//
//  auto error = sbp.Dx(f) - sbp.Evaluate(dF);
//  std::cout << sqrt(sbp.Integrate(error*error)) << "\n";
//  ExportToTec(sbp.grid(), error, "test");
//
//  return 0;
//}
//
Array Linspace(double x0, double x1, int N) {
   Array out(1, N);
   double h = (x1 - x0)/(N-1);
   for(int i = 0; i < N; i++)
      out[i] = x0 + i*h;
   return out;
};

Array Ones(int N) {
  Array out(1, N);
  for(int i = 0; i < N; i++)
     out[i] = 1;
  return out;
};


int main() {
  int N_fine = 23;
  int N_coarse = 12;

  Array x_fine = Linspace(0, 1, N_fine);
  Array x_coarse = Linspace(0, 1, N_coarse);

  Interp21 interp{N_coarse,N_fine};
  Array err = x_coarse -
                 interp.Interpolate(x_fine);
  Print(err);
   return 0;
}
