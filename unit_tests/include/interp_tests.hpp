#include <memory> // for unique_ptr
#include "array.hpp"
#include "interp.hpp"
#include "interp21.hpp"
#include "interp42.hpp"

/*
 * Test Interpolation operators:
 *  o TestInterpolation
 */

Array Linspace(double x0, double x1, int N) {
   Array out(1, N);
   double h = (x1 - x0)/(N-1);
   for(int i = 0; i < N; i++)
      out[i] = x0 + i*h;
   return out;
};

int TestInterpolation(int Nc, int Nf, int order) {

  auto x_course = Linspace(0,1,Nc);
  auto x_fine = Linspace(0,1,Nf);

  std::unique_ptr<Interp> interp;

  int stop = 1;
  std::vector<Array> f_coarse{Ones(1,Nc)};
  std::vector<Array> f_fine{Ones(1,Nf)};
  switch (order) {
    case 4:
      interp = std::make_unique<Interp42>(Nc, Nf);
      stop = 2;
      f_coarse.push_back(Linspace(0,1,Nc));
      f_fine.push_back(Linspace(0,1,Nf));
      break;
    default:
      interp = std::make_unique<Interp21>(Nc, Nf);
      break;
  }

  double error = -1;
  for(int i = 0; i < stop; ++i) {
    auto err_vec = f_fine[i] - interp->Interpolate(f_coarse[i]);
    double error_c2f = sqrt(Sum(Abs(err_vec*err_vec)));
    err_vec = f_coarse[i] - interp->Interpolate(f_fine[i]);
    double error_f2c = sqrt(Sum(Abs(err_vec*err_vec)));
    error = std::max(error_c2f, error_f2c);
  }
  return (error < 1e-14);
};

