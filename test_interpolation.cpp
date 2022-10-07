#include <iostream>
#include <chrono>
#include <string>
#include <fstream>
#include <cmath>
#include <valarray>
#include "array.hpp"
#include "mbarray.hpp"
#include "mesh.hpp"
#include "mbgrid.hpp"
#include "mbsbp.hpp"
#include "interpolation.h"
#include "interpolation21.h"
#include "interpolation42.h"

Array2D Linspace(double x0, double x1, int N)
{
   Array2D out(1, N);
   double h = (x1 - x0)/(N-1);
   for(int i = 0; i < N; i++)
      out[i] = x0 + i*h;
   return out;
};

Array2D Ones(int N)
{
   Array2D out(1, N);
   for(int i = 0; i < N; i++)
      out[i] = 1;
   return out;
};


int main()
{
   int N_fine = 23;
   int N_coarse = 12;

   Array2D x_fine = Linspace(0, 1, N_fine);
   Array2D x_coarse = Linspace(0, 1, N_coarse);
   x_fine = x_fine*x_fine;
   x_coarse = x_coarse*x_coarse;
   //x_fine = Ones(N_fine);
   //x_coarse = Ones(N_coarse);

   Interpolation42 interp = Interpolation42(N_coarse,N_fine);
//   Array2D e_fine = x_fine - interp.Interpolate(x_coarse, N_coarse, N_fine);
   interp.Interpolate(x_fine, N_fine, N_coarse).Print();
 //  Array2D e_coarse = x_coarse - interp.Interpolate(x_fine, N_fine, N_coarse);

//    e_coarse.Print();
   //std::cout << "\n Coarse 2 Fine: " << std::endl;
//   e_fine.Print();


   return 0;
}
