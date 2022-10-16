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
#include "advection.hpp"
#include "integrate.hpp"
#include "basics.hpp"


int main() {
  int N = 21;
  auto blocks {NonConformingAnnulus(N, 0.2, 1)};

  MbGrid grid {blocks};

  auto initu = [] (double x, double y) {
    return 1.0;
  };

  MbArray f{grid.Evaluate(initu)};
  ExportToTec(grid, f, "annulus_grid");
  return 0;

}
