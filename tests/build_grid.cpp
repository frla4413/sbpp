//============================================================================
// Name        : build_grid.cpp
// Author      : frela05
// Description : Use this routine to build grids with several blocks. 
//               View the result in ParaView.
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
#include "ins_setup.hpp"

int main() {
  int N = 41;

  auto setup = Channel(N);
  MbGrid grid {setup.blocks};

  auto initu = [] (double x, double y) {
    return 1.0;
  };

  MbArray f{grid.Evaluate(initu)};
  ExportToTec(grid, f, "channel");
  std::cout << "Total size " << f.GetTotalSize() << "\n";
  Print(grid.boundaries());
  Print(grid.interfaces());
  return 0;
}
