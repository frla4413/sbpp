//============================================================================
// Name        : main_explicit.cpp
// Author      : frela05
// Description : main-file for advection equation u_t + a1 ux + a2 uy = 0
//               This file is used to test how discontinuous initial data
//               affect explicit and implicit time integration
//               The solution is saved to "adv/sol".
//               This file must be created before execution.
//============================================================================
#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <math.h>
#include <valarray>
#include "array2D.h"
#include "gridfunction.h"
#include "mesh.h"
#include "multiblockgrid.h"
#include "multiblocksbp.h"
#include "data_types.h"
#include "basics.h"
#include "advection.h"
#include "integrate.h"

Gridfunction InitialData(MultiblockSbp& sbp) 
{ 
   auto init = [] (double x, double y) 
   { 
      double u;
      if (x >= 0){ u = 0.5; }
      else{ u =  1 + exp(-pow((x + 0.5)/0.1,2))*exp(-pow((y - 0.5)/0.1,2));};
      return u;
   };
   return sbp.EvaluateFunction(init);
}

int main()
{
  double t0 = 0;
  double t1 = 10;

  std::vector<double> tspan = {t0,t1};

  std::pair<double,double> a = std::make_pair(1,0);

   int Nx = 40;
   int Ny = 40;
  int order = 4;

   int err_vec_pos = 0;
   Gridfunction computed;

   auto block = Mesh::CartesianGrid(Nx,Ny);
   auto x1 = block.x - 1;

   DataTypes::Block block1 = DataTypes::SetBlock(x1, block.y);
   auto blocks = {block,block1};

   MultiblockSbp sbp = MultiblockSbp(blocks, order);
  
   MultiblockSbp* sbp_ptr = &sbp;
  Advection advec = Advection(sbp_ptr,a);

  Gridfunction initial = InitialData(sbp);

  double h  = 1/(double)(Nx-1);
   double dt = 0.1*h;

   auto rhs_fun =  [&] (double t, const Gridfunction& y)
   {
      bool mms;
      return advec.RhsFunc(t,y, mms = false);
   };

   DataTypes::Solution solution;
   auto ExportToTec = [&] (const Gridfunction& f, std::string name)
   {
      advec.ExportToTec(f, name);
   };

   std::string name_base = "adv/sol";
//   solution = Integrate::ExplicitIntegrationSaveSolution(tspan, 
//                                          initial, dt, rhs_fun, ExportToTec,
//                                          name_base);


   auto lhs_fun =  [&] (double t, const Gridfunction& y)
   {
      bool mms;
      return -1.0*advec.RhsFunc(t,y, mms = false);
   };

   auto Jv_fun =  [&] (double t, const Gridfunction& y)
   {
      return -1.0*advec.Jacobian(t,y);
   };

   solution = Integrate::ImplicitTimeIntegractionSaveSolution(tspan, initial,
                                                  dt, lhs_fun, Jv_fun, 
                                                  ExportToTec, name_base);


   std::cout << "done\n";
   return 0;
}
