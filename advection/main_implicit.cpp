//============================================================================
// Name        : main_implicit.cpp
// Author      : frela05
// Description : main-file for advection equation u_t + a1 ux + a2 uy = F
//               This test case includes a manufactured solution.
//               Implicit time integration is used.
//               The solution of the finest grid is saved to "adv/sol".
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

int main()
{
	double t0 = 0;
	double t1 = 5;

	std::vector<double> tspan = {t0,t1};

	std::pair<double,double> a = std::make_pair(1.0,-1.0);

   std::valarray<int> N_vec = {21,41,61};
   std::valarray<double> err_vec(N_vec.size());
   int err_vec_pos = 0;

	int order = 4;
   for(auto& N : N_vec)
   {
//      auto blocks = Mesh::GetAnnulus(N, 0.1, 1.0);
      auto blocks = Mesh::GetChannelGrid(N, N);

      MultiblockSbp sbp = MultiblockSbp(blocks, order);
	   
      MultiblockSbp* sbp_ptr = &sbp;
	   Advection advec = Advection(sbp_ptr,a);

	   Gridfunction initial = advec.AnalyticVec(0.0);
      
      double h  = 1/(double)(N-1);
      double dt = h;

      auto lhs_fun =  [&] (double t, const Gridfunction& y)
      {
         bool mms;
         return -1.0*advec.RhsFunc(t,y, mms = true);
      };

      auto Jv_fun =  [&] (double t, const Gridfunction& y)
      {
         return -1.0*advec.Jacobian(t,y);
      };

      DataTypes::Solution solution;

      // save finest grid solution only
      if(N == N_vec[N_vec.size()-1])
      {
         auto ExportToTec = [&] (const Gridfunction& f, std::string name)
         {
            advec.ExportToTec(f, name);
         };

         std::string name_base = "adv/sol";
         solution = Integrate::ImplicitTimeIntegractionSaveSolution(tspan, initial,
                                                        dt, lhs_fun, Jv_fun, 
                                                        ExportToTec, name_base);
      }
      else
      {
         solution = Integrate::ImplicitTimeIntegraction(tspan, initial,
                                                        dt, lhs_fun, Jv_fun);
      }
	   err_vec[err_vec_pos] = advec.ComputeError(solution.time,solution.sol);

	   std::cout << "Error: " << err_vec[err_vec_pos] << "\n";
      err_vec_pos ++;
   }
   std::cout << "done\n";
   Basics::PrintConvergence(N_vec, err_vec);
   return 0;
}
