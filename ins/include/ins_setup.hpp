//==================================================================
// Name        : ins_setp.hpp
// Author      : Fredrik Laurén
// Description : Setup file for ins_main.cpp. Build the grid and
//               specify boundary conditions.
//               The file build_tests/grid.cpp can be used 
//               when building the grid.
//==================================================================


#include "ins.hpp"
#include "mesh.hpp"
#include <algorithm>
#include <cmath>
#include <valarray>

struct InsSetup {
  std::vector<Block> blocks;
  std::vector<BdType> bd_types;
};

/*
 * This is a channel with an obstacle in the middle.
 * For a high Raynolds numbers (small viscisity coefficient, mu),
 * the flow will produce vortices.
 */

InsSetup GetObstacleGrid(int N) {

 /*
  * a4   ---------------------------------------------------------
  *      |                 |    |   |    |
  *      |      0          | 1  | 2 | 3  | 4
  * a3   ---------------------------------------------------------
  *      |                 |    |   |    |
  *      |      5          | 6  | 7 | 8  | 9
  * a2   ---------------------------------------------------------
  *      |                 |    | | |    |
  *      |      10         | 11 | | | 12 | 13
  * -a2  ----------------------------------------------------------
  *      |                 |    |   |    |
  *      |      14         | 15 |16 | 17 | 18
  * -a3  ----------------------------------------------------------
  *      |                 |    |   |    |
  *      |      19         | 20 |21 | 22 | 23
  * -a3  ----------------------------------------------------------
  *      b0                b1   b2   b3   b4                     b5
  */

  std::vector<BdType> bd_types(24,BdType::Pressure);

  for(auto& i : {1,7,10,14,18})
    bd_types[i] = BdType::Inflow;

  for(auto& i : {5,9,13,16,23})
    bd_types[i] = BdType::Outflow;

  for(auto& i : {8,11,12,15})
    bd_types[i] = BdType::Wall;

  double a0 = -1.25, a1 = -0.15, a2 = 0.15, a3 = 0.4, a4 = 1.25;
  double b0 = -2, b1 = -0.5, b2 = -0.1, b3 = 0.1,
         b4 = 0.5, b5 = 3.5;
  int N1 = N, N2 = N1/2+1, N3 = N2/2+1, N4 = N3/2+1, N5 = N4/2+1;

  // Upper row
  Block block0 = CartesianGrid(N4, N4);
  block0.x = block0.x*std::abs(b0-b1) + b0;
  block0.y = block0.y*std::abs(a4-a3) + a3;
  std::vector<Block> blocks = {block0};

  Block block1 = CartesianGrid(N4, N3);
  block1.x = block1.x*std::abs(b1-b2) + b1;
  block1.y = block1.y*std::abs(a4-a3) + a3;
  blocks.push_back(block1);

  Block block2 = CartesianGrid(N4, N3);
  block2.x = block2.x*std::abs(b2-b3) + b2;
  block2.y = block2.y*std::abs(a4-a3) + a3;
  blocks.push_back(block2);

  Block block3 = CartesianGrid(N4, N3);
  block3.x = block3.x*std::abs(b4-b3) + b3;
  block3.y = block3.y*std::abs(a4-a3) + a3;
  blocks.push_back(block3);

  Block block4 = CartesianGrid(N3, N4);
  block4.x = block4.x*std::abs(b5-b4) + b4;
  block4.y = block4.y*std::abs(a4-a3) + a3;
  blocks.push_back(block4);

  //Upper-middle row
  Block block5 = CartesianGrid(N3, N4);
  block5.x = block5.x*std::abs(b0-b1) + b0;
  block5.y = block5.y*std::abs(a3-a2) + a2;
  blocks.push_back(block5);

  Block block6 = CartesianGrid(N3, N3);
  block6.x = block6.x*std::abs(b1-b2) + b1;
  block6.y = block6.y*std::abs(a3-a2) + a2;
  blocks.push_back(block6);

  Block block7 = CartesianGrid(N3, N3);
  block7.x = block7.x*std::abs(b2-b3) + b2;
  block7.y = block7.y*std::abs(a3-a2) + a2;
  blocks.push_back(block7);

  Block block8 = CartesianGrid(N3, N3);
  block8.x = block8.x*std::abs(b4-b3) + b3;
  block8.y = block8.y*std::abs(a3-a2) + a2;
  blocks.push_back(block8);

  Block block9 = CartesianGrid(N2, N4);
  block9.x = block9.x*std::abs(b5-b4) + b4;
  block9.y = block9.y*std::abs(a3-a2) + a2;
  blocks.push_back(block9);

  // Middle row
  Block block10 = CartesianGrid(N2, N4);
  block10.x = block10.x*std::abs(b0-b1) + b0;
  block10.y = block10.y*std::abs(a2-a1) + a1;
  blocks.push_back(block10);

  Block block11 = CartesianGrid(N3, N3);
  block11.x = block11.x*std::abs(b1-b2) + b1;
  block11.y = block11.y*std::abs(a2-a1) + a1;
  blocks.push_back(block11);

  Block block12 = CartesianGrid(N3, N3);
  block12.x = block12.x*std::abs(b4-b3) + b3;
  block12.y = block12.y*std::abs(a2-a1) + a1;
  blocks.push_back(block12);

  Block block13 = CartesianGrid(N1, N4);
  block13.x = block13.x*std::abs(b5-b4) + b4;
  block13.y = block13.y*std::abs(a2-a1) + a1;
  blocks.push_back(block13);

  // bottom middle row
  for(int k = 5; k < 10; ++k) {
    auto block = blocks[k];
    block.y = block.y - a3 - a2;
    blocks.push_back(block);
  }

  // Bottom row
 for(int k = 0; k < 5; ++k) {
   auto block = blocks[k];
   block.y = block.y - a4 - a3;
   blocks.push_back(block);
 }
 return {blocks, bd_types};
}

/*
 * A functiont that stretches the coordinates so that more 
 * points are located at the boundaries.
 */
Array StretchCoordinates(const Array& x) {

  auto f = x.array();
  f = f*2;
  f = f-1;
  f = (std::atan(2*f) + M_PI/2)/M_PI;
  f = f - std::abs(f).min();
  f = f/std::abs(f).max();
  return {Array(x.Nx(), x.Ny(), f)};
}

/*
 * A channel that consists of four blocks.
 */
InsSetup Channel(int N) {

  auto block {CartesianGrid(N, N)};
  block.y = StretchCoordinates(block.y);
  block.x = StretchCoordinates(block.x);
  std::vector<Block> blocks {block};
  block.x = block.x+0.5;
  block.x = block.x*2;
  blocks.push_back(block);
  block.x = block.x-3;
  blocks.push_back(block);
  block = CartesianGrid(N, N);
  block.x = StretchCoordinates(block.x);
  block.y = StretchCoordinates(block.y);
  block.y = block.y-1;
  blocks.push_back(block);

  std::vector<BdType> bd_types(10,BdType::Wall);

  for(auto& i : {6})
    bd_types[i] = BdType::Inflow;

  for(auto& i : {3})
    bd_types[i] = BdType::Pressure;

  return {blocks, bd_types};
}
