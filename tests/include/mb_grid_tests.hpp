#include "mbgrid.hpp"

/*
 * Test the interfaces of an annulus.
 */
bool TestMbGridAnnulusInterfaces() {
  int N = 8;
  auto blocks {Annulus(N, 0.1, 1)};
  MbGrid grid {blocks};
  Interface i0 {0,1, Side::n, Side::s};
  Interface i1 {0,3, Side::s, Side::n};
  Interface i2 {1,2, Side::n, Side::s};
  Interface i3 {2,3, Side::n, Side::s};
  std::vector<Interface> interfaces {i0, i1, i2, i3};

  for(int i = 0; i < grid.interfaces().size(); ++i) {
    auto FindIntf = [&i0](const Interface& intf_)
               { return Equal(i0, intf_);};

    auto it = std::find_if(grid.interfaces().begin(),
                           grid.interfaces().end(),
                           FindIntf);

    auto dist = std::distance(grid.interfaces().begin(), it);
    if (dist >= grid.interfaces().size())
      return false;
  }
  return true;
}

/*
 * Test the interfaces of a square that consists of four 
 * small squares.
 */
bool TestMbGridSquareInterfaces() {
  int N = 3;
  auto block = CartesianGrid(N,N);
  std::vector<Block> blocks{block};
  block.x += 1;
  blocks.push_back(block);
  block.y += 1;
  blocks.push_back(block);
  block.x -= 1;
  blocks.push_back(block);
  MbGrid grid {blocks};

  Interface i0 {0,1, Side::e, Side::w};
  Interface i1 {0,3, Side::n, Side::s};
  Interface i2 {1,2, Side::n, Side::s};
  Interface i3 {2,3, Side::w, Side::e};
  std::vector<Interface> interfaces {i0, i1, i2, i3};

  for(int i = 0; i < grid.interfaces().size(); ++i) {
    auto FindIntf = [&i0](const Interface& intf_)
               { return Equal(i0, intf_);};

    auto it = std::find_if(grid.interfaces().begin(),
                           grid.interfaces().end(),
                           FindIntf);

    auto dist = std::distance(grid.interfaces().begin(), it);
    if (dist >= grid.interfaces().size())
      return false;
  }
  return true;
}
