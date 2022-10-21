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
      return 1.0;
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

int main()
{

   int Nx = 41;
   int Ny = 81;
   auto blocks = Mesh::GetMultiblockCartesianGrid(Nx, Ny);

   std::vector<DataTypes::BdType> bd_types(8);
   bd_types[0] = DataTypes::BdType::SlipWall;
   bd_types[1] = DataTypes::BdType::Pressure;
   bd_types[2] = DataTypes::BdType::Inflow;
   bd_types[3] = DataTypes::BdType::Wall;
   bd_types[4] = DataTypes::BdType::Pressure;
   bd_types[5] = DataTypes::BdType::Wall;
   bd_types[6] = DataTypes::BdType::Pressure;
   bd_types[7] = DataTypes::BdType::Pressure;

   MultiblockGrid grid = MultiblockGrid(blocks);

//   Basics::PrintBoundaries(grid.GetBoundaries());

   int order = 4;

   MultiblockSbp sbp = MultiblockSbp(blocks, order);

   MultiblockSbp* sbp_ptr = &sbp;
   double mu = 1e-4;
   Ins ins = Ins(sbp_ptr, bd_types, mu);
   DataTypes::InsState init = InitialData(sbp);

   std::string name_base = "plate/sol";
   DataTypes::InsSolution sol = IntegrateIns::SolveSteadyState(
                                init, ins, name_base);

   int block_idx = 2;
   auto side = DataTypes::Side::w;
   ins.ExportBoundaryToTec(sol.state, "west", block_idx, side);
   side = DataTypes::Side::e;
   ins.ExportBoundaryToTec(sol.state, "east",block_idx, side);
   side = DataTypes::Side::s;
   ins.ExportBoundaryToTec(sol.state, "south", block_idx, side);
   side = DataTypes::Side::n;
   ins.ExportBoundaryToTec(sol.state, "north", block_idx, side);

   ExportYSlice(sol.state, sbp_ptr, "sol_slice", block_idx, 0);

   return 0;
}
