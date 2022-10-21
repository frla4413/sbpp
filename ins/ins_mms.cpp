//============================================================================
// Name        : ins_main.cpp
// Author      : frela05
// Description : main-file for the INS equations.
//               Implicit time integration is used.
//               This folder "sol" must be created before execution.
//============================================================================
#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <cmath>
#include <valarray>
#include "array2D.h"
#include "gridfunction.h"
#include "mesh.h"
#include "multiblockgrid.h"
#include "multiblocksbp.h"
#include "data_types.h"
#include "basics.h"
#include "ins.h"
#include "integrate_ins.h"
#include <cstring>
#include <sstream>
#include <iostream>
#include <string>

DataTypes::InsState InitialData(MultiblockSbp& sbp)
{
   auto initu = [] (double x, double y)
   {
      return 1.0 + exp(-pow((x - 0.5)/0.1,2))*exp(-pow((y - 0.5)/0.1,2));
   };
   auto initv = [] (double x, double y) { return 0; };
   auto initp = [] (double x, double y) { return 0; };
   DataTypes::InsState init;
   init.u = sbp.EvaluateFunction(initu);
   double normalize = abs(init.u.GridfunctionToValarray()).max();
   init.u = init.u/normalize;
   init.v = sbp.EvaluateFunction(initv);
   init.p = sbp.EvaluateFunction(initp);
   return init;
}

// TO DO:
// Clean up jacobians
// MKL-operations for Array2D? -- try on advection and compare operations

int main()
{

   double t0 = 0, t1 = 5, dt = 0.025;

   std::vector<double> tspan = {t0,t1};
   int N = 41;
   auto blocks = Mesh::NonConformingGrid(N);

   std::vector<DataTypes::BdType> bd_types(6,DataTypes::BdType::Pressure);

   for(auto& i : {2})
     bd_types[i] = DataTypes::BdType::Inflow;

   int order = 4;
   MultiblockSbp sbp = MultiblockSbp(blocks, order);

   MultiblockSbp* sbp_ptr = &sbp;
   double mu = 0;
   Ins ins = Ins(sbp_ptr, bd_types, mu);
   DataTypes::InsState init = InitialData(sbp);

   std::string name_base = "nonconforming/sol";
//   DataTypes::InsSolution sol = IntegrateIns::SolveSteadyState(init, ins, name_base);
   DataTypes::InsSolution sol = IntegrateIns::ImplicitTimeIntegractionSaveSolution(tspan, init, dt, ins, name_base);
   return 0;
}
