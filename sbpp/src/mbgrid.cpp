/*
 * multiblockgrid.cpp
 *
 *  Created   : August 2021
 *      Author: Fredrik Laurén
 */

#include "mbgrid.hpp"
#include <assert.h>

//alias to improve readability
using Shapes = std::vector<Shape>;
using Arrays = std::vector<Array>;
using Blocks = std::vector<Block>;
using Corners = std::vector<Corner>;


Corners GetCorners(const Block& block) {
/*
 * Return the corners of a block.
 * Start with x[0], y[0] and continues counter-clockwise
 */
  std::vector<Corner> corners(4);
  int Nx  = block.x.Nx();
  int Ny  = block.x.Ny();
  // south-west
  corners[0] = {block.x[0], block.y[0]};
  // south-east
  corners[1] = {block.x[Nx*Ny-Ny], block.y[Nx*Ny-Ny]};
  // north-west
  corners[2] = {block.x[Nx*Ny-1], block.y[Nx*Ny-1]};
  // north-west
  corners[3] = {block.x[Ny-1], block.y[Ny-1]};
  return corners;
}

bool CompareCornerSize(const Corner& lhs, const Corner& rhs) {
  // is lhs < rhs?
  if(lhs.Size() == rhs.Size()) {
  // same  size, x-component decides
    if (lhs.x < rhs.x)
      return true;
    else
      return false;
  }
  return lhs.Size() < rhs.Size();
}

bool IsSameCorner(const Corner& lhs, const Corner& rhs) {
  return (lhs.x == rhs.x && lhs.y == rhs.y);
};

std::string SideToString(const Side& side) {
  switch (side) {
    case Side::s: {return "s";};
    case Side::e: {return "e";};
    case Side::n: {return "n";};
    case Side::w: {return "w";};
    default: return {"Error SideToString"};
  }
}

int SideToInt(const Side& side) {
  switch (side) {
    case Side::s: {return 0;};
    case Side::e: {return 1;};
    case Side::n: {return 2;};
    case Side::w: {return 3;};
    default: return -1;
  }
}

void Print(const Interface& interface) {
  std::string print_me = "(";
  print_me += std::to_string(interface.block_idx1) + " ";
  print_me += SideToString(interface.side1);
  print_me += " : ";
  print_me += std::to_string(interface.block_idx2) + " ";
  print_me += SideToString(interface.side2);
  print_me += ")";
  std::cout << print_me << std::endl;
}

void Print(const std::vector<Interface>& interfaces) {
  int interface_idx = 0;
  for(auto& interface : interfaces) {
    std::cout << "Interface " << interface_idx << ": ";
    Print(interface);
    ++interface_idx;
  }
}

bool Equal(const Interface& lhs, const Interface& rhs) {

  bool same_blocks = false;

  if ((lhs.block_idx1 == rhs.block_idx1 &&
       lhs.block_idx2 == rhs.block_idx2) ||
      (lhs.block_idx1 == rhs.block_idx2 &&
       lhs.block_idx2 == rhs.block_idx1))
    same_blocks = true;

  bool same_sides = false;

  if ((lhs.side1 == rhs.side1 &&
       lhs.side2 == rhs.side2) ||
      (lhs.side1 == rhs.side2 &&
       lhs.side2 == rhs.side1))
    same_sides = true;

  return (same_blocks && same_sides);
}

MbGrid::MbGrid(const Blocks& blocks) {

  for(auto& block : blocks) {
    if(GetShape(block.x) != GetShape(block.y))
      throw std::invalid_argument("@MbGrid invalid shape");
  }

  num_blocks_ = blocks.size();
  blocks_ = blocks;

  for (auto& block : blocks)
     shapes_.push_back(GetShape(block.x));

  SetGrid();
}

int MbGrid::num_blocks() {
  return num_blocks_;
};

Blocks& MbGrid::blocks() {
  return blocks_;
};

const Blocks& MbGrid::blocks() const {
  return blocks_;
};

Block& MbGrid::blocks(int block_idx) {
  return blocks_[block_idx];
};

Shapes& MbGrid::shapes() {
  return shapes_;
};

const Shapes& MbGrid::shapes() const {
   return shapes_;
};

Shape& MbGrid::shapes(int block_idx) {
  return shapes_[block_idx];
};

int MbGrid::GetTotalSize() const {
  int tot_size = 0;
  for(auto& shape : shapes_)
     tot_size += shape.Nx * shape.Ny;
  return tot_size;
}

const std::vector<Interface>& MbGrid::interfaces() const {
  return interfaces_;
};

std::vector<std::vector<BlockInterface>>MbGrid::GetBlockInterfaces
() {
  return block_interfaces_;
};

std::vector<BoundaryInfo> MbGrid::boundaries() {
  return boundaries_;
};
/* Check if an interface has flipped orientation compared to its neighbour.
 * Ex. an east-to-south interface is flipped.
 *
 * Input: 
 *    o interface_idx: The index of the interface
 * Output:
 *    o True if flipped, False otherwise
 */
bool MbGrid::IsFilppedInterface(int interface_idx) {

  bool is_flipped = false;
  Side side1 = interfaces_[interface_idx].side1;
  Side side2 = interfaces_[interface_idx].side2;

  std::pair<Side,Side> sides = std::make_pair(side1,side2);

  auto FindInterface = [&sides](std::pair<Side,Side>& sides_) {
    return sides_.first == sides.first &&
           sides_.second == sides.second;
  };
  auto it = std::find_if(flipped_interfaces_.begin(),
                         flipped_interfaces_.end(),
                         FindInterface);

  auto dist = std::distance(flipped_interfaces_.begin(), it);

  if(dist < flipped_interfaces_.size())
     is_flipped = true;
  return is_flipped;
}

Block MbGrid::GetBoundary(int block_idx, Side side) {

   auto f{blocks_[block_idx].x};
   auto shapes {GetShape(f)};
   Block boundary;
   boundary.x = ToBlockBoundary({f}, block_idx, side);
   f = blocks_[block_idx].y;
   boundary.y = ToBlockBoundary({f}, block_idx, side);
   return boundary;
}

Array MbGrid::ToBlockBoundary(const MbArray& f,
                              int block_idx, Side side) {
  std::pair<int, std::slice> bd =
  GetBlockBoundarySliceAndSize(block_idx,side);
  return Array(1,bd.first,f[block_idx][bd.second]);
}

MbArray MbGrid::Evaluate(std::function<double(double,double)> f) {

  std::vector<Array> z;
  int block_idx = 0;

  for(auto& block : blocks_) {
    Shape shape = shapes(block_idx);
    Array zi = Array(shape.Nx, shape.Ny);

    for(int j = 0; j < block.x.size(); ++j)
       zi[j] = f(block.x[j], block.y[j]);
    z.push_back(zi);
    ++block_idx;
  }
  return MbArray(z);
}

/* 
 * Retrurn a slice-object and the size to be used in ToBlockBoundary.
 * Input:
 *    o block_idx - the block index
 *    o side - the side, one of {s,e,n,w}
 *  Returns: std::pair - the size of the boundary (i.e. Nx or Ny) and the slice object.
 */
std::pair<int,std::slice> MbGrid::
GetBlockBoundarySliceAndSize(int block_idx, Side side) {

  int start  = 0, stride = 0, size = 0;

  Shape shape{shapes(block_idx)};
  int Nx = shape.Nx;
  int Ny = shape.Ny;

  switch (side) {
    case Side::w: {start = 0;  stride = 1; size = Ny; break;}
    case Side::e: {start = Nx*Ny-Ny;stride = 1; size = Ny; break;}
    case Side::s: {start = 0; stride = Ny; size = Nx; break;}
    case Side::n: {start = Ny-1; stride = Ny; size = Nx; break;}
  }
 return std::make_pair(size,std::slice(start,size,stride));
}

Array MbGrid::ToBlockBoundary(const Array& f,
                              int block_idx,
                              Side side) {
  std::pair<int,std::slice> bd =
    GetBlockBoundarySliceAndSize(block_idx,side);
   return Array(1,bd.first,f[bd.second]);
}

// ----------- private function, used in the constructor ---------

void MbGrid::SetGrid() {

  CollocateCorners();
  SaveUniqueCorners();
  InterfacesAsCorners();
  SaveUniqueEdges();
  SaveFaceEdges();
  FindInterfaces();
  FindPhysicalBoundaries();
  SetFlippedInterfaces();
}

void MbGrid::SaveUniqueCorners() {
  /*
   * Extract all corners from each block.
   * Sort according to CompareCornerSize
   * Save unique corners only
  */

  std::vector<Corner> block_corners;

  for(auto& block : blocks_) {
    block_corners = GetCorners(block);
    for (auto& corner : block_corners)
      corners_.push_back(corner);
  }

  std::sort(corners_.begin(),
            corners_.end(),
            CompareCornerSize);
  auto it = std::unique(corners_.begin(),
                        corners_.end(),
                        IsSameCorner);
  corners_.resize(std::distance(corners_.begin(), it) );
}

void MbGrid::InterfacesAsCorners() {
 /*
  * Save faces in terms of unique corners
  * faces_[i] contains the faces of the i:th block.
  * Ex. faces_[0] = {0,2,3,1} means that element {0,2,3,1}
  * in corners_
  * are the corners of block 0.
  */

  std::vector<int> indices;
  for(auto& block : blocks_) {
    auto block_corners = GetCorners(block);
    for (auto& corner : block_corners) {
      auto FindCorner = [corner](Corner corner_) {
        return corner_.x == corner.x &&
               corner_.y == corner.y; };

      auto it = std::find_if(corners_.begin(),
                             corners_.end(),
                             FindCorner);
      indices.push_back(std::distance(corners_.begin(), it));
     }
     faces_.push_back(indices);
     indices.clear();
  }
}

void MbGrid::SaveUniqueEdges() {
  std::vector<int> index_vec(2);
  for(auto& face : faces_) {
    for(int k = 0; k < 4; k++) {
      index_vec[0] = face[k];
      index_vec[1] = face[(k+1)%4];
      std::sort(index_vec.begin(), index_vec.begin()+2);
      edges_.push_back(std::make_pair(index_vec[0],index_vec[1]));
    }
    index_vec.clear();
  }
  std::sort(edges_.begin(), edges_.end(), CompareEdgeSize);
  auto it = std::unique(edges_.begin(), edges_.end(), IsSameEdge);
  edges_.resize(std::distance(edges_.begin(), it));
}

// to be used in SaveUniqueEdges
bool MbGrid::CompareEdgeSize(const std::pair<int,int>& lhs,
                             const std::pair<int,int>& rhs) {
  // is lhs < rhs?
  int lhs_size = lhs.first*lhs.first + lhs.second*lhs.second;
  int rhs_size = rhs.first*rhs.first + rhs.second*rhs.second;
  if(lhs_size == rhs_size) {
  // same  size, x-component decides
    if (lhs.first < rhs.first)
      return true;
     else
      return false;
  }
  return lhs_size < rhs_size;
}

// to be used in SaveUniqueEdges
bool MbGrid::IsSameEdge(const std::pair<int,int>& lhs,
                        const std::pair<int,int>& rhs) {
  return lhs.first == rhs.first && lhs.second == rhs.second;
}

void MbGrid::SaveFaceEdges() {
  /*
   * faces_edges_[i] contains the faces of the i:th block 
   * in the specific order {s, e, n ,w}.
   * The difference to faces_ is that they are not
   * in the specific order {s, e, n ,w}.
   * face_edges_[k] is a vector<pair<Side, int>>.
   * Ex. face_edges_[0] = {[s,1], [e,3], [n,2], [w,0]},
   * where the indices {1,3,2,0} corresponds to the
   * elements in edges_.
   */
  std::pair<Side, int> face_edge;
  std::vector<std::pair<Side, int>> block_face_edges;
  std::vector<Side> sides{Side::s, Side::e, Side::n, Side::w};

  std::vector<int> index_vec(2);
  std::pair<int, int> edge;

  for(auto& face : faces_) {
    for(int k = 0; k < 4; k++) {
      index_vec[0] = face[k];
      index_vec[1] = face[(k+1)%4];
      std::sort(index_vec.begin(), index_vec.begin()+2);
      edge = std::make_pair(index_vec[0], index_vec[1]);
      auto FindEdge = [edge](std::pair<int,int> edge_) {
        return edge_.first == edge.first &&
               edge_.second == edge.second;
      };
      auto it = std::find_if(edges_.begin(),
                             edges_.end(),
                             FindEdge);
      auto dist = std::distance(edges_.begin(), it);
      block_face_edges.push_back(std::make_pair(sides[k], dist));
    }
    face_edges_.push_back(block_face_edges);
    block_face_edges.clear();
  }
}

void MbGrid::FindPhysicalBoundaries() {
  std::vector<Side> sides{Side::s, Side::e, Side::n, Side::w};

  BoundaryInfo boundary;
  for(int block_idx = 0; block_idx < num_blocks_; block_idx++) {
    for(auto& side : sides) {
      if(!IsInterface(block_idx, side)) {
         boundary.block_idx = block_idx;
         boundary.side = side;
         boundaries_.push_back(boundary);
      }
    }
  }
}

bool MbGrid::IsInterface(int block_idx, Side side) {
  bool is_interface = false;
  for(auto& block_interface : block_interfaces_[block_idx]) {
    if(block_interface.my_side == side) {
      is_interface = true;
      break;
    }
  }
  return is_interface;
}

void MbGrid::FindInterfaces() {
  std::vector<Side> sides = {Side::s, Side::e, Side::n, Side::w};

  BlockInterface block_interface1;
  BlockInterface block_interface2;

  std::vector<std::vector<BlockInterface>>
    block_interfaces(num_blocks_);

  Interface interface;

  for(int i = 0; i < num_blocks_; i++) {
    auto faces_block1 = face_edges_[i];
    for(int j = i+1; j < num_blocks_; j++) {
      auto faces_block2 = face_edges_[j];
      for(int m = 0; m < 4; m++) {
        auto edge1 = faces_block1[m].second;
        for(int n = 0; n < 4; n++) {
          auto edge2 = faces_block2[n].second;
          if(edge1 == edge2) {
            interface.block_idx1 = i;
            interface.side1 = faces_block1[m].first;
            interface.block_idx2 = j;
            interface.side2 = faces_block2[n].first;
            interfaces_.push_back(interface);

            block_interface1.my_side = interface.side1;
            block_interface1.adjacent_block_idx = j;
            block_interface1.adjacent_block_side = interface.side2;

            block_interface2.my_side = interface.side2;
            block_interface2.adjacent_block_idx = i;
            block_interface2.adjacent_block_side = interface.side1;

            block_interfaces[i].push_back(block_interface1);
            block_interfaces[j].push_back(block_interface2);
          }
        }
      }
    }
  }
  block_interfaces_ = block_interfaces;
}

/*
 * Collocate corners of blocks if they are equal up to some tolerance
 */
void MbGrid::CollocateCorners() {
  double tol = 1e-12;
  Block block1, block2;
  std::vector<Corner> corners_block1, corners_block2;

  std::valarray<int> ind_vec1, ind_vec2;
  Shape shape_block1, shape_block2;
  int Nx, Ny, ind1, ind2;

  for(int i = 0; i < blocks_.size(); i++) {
    block1 = blocks_[i];
    shape_block1 = shapes(i);
    Nx = shape_block1.Nx;
    Ny = shape_block1.Ny;
    ind_vec1 = {0, Nx*Ny-Ny, Nx*Ny-1,Ny-1};
    corners_block1 = GetCorners(block1);

    for(int j = i; j < blocks_.size(); j++) {
      block2 = blocks_[j];
      shape_block2 = shapes(j);
      Nx = shape_block2.Nx;
      Ny = shape_block2.Ny;
      ind_vec2 = {0, Nx*Ny-Ny, Nx*Ny-1,Ny-1};
      corners_block2 = GetCorners(block2);

      ind1 = 0;
      for(auto& corner_block1 : corners_block1) {
        ind2 = 0;
        for(auto& corner_block2 : corners_block2) {
          if(abs(corner_block1.x  - corner_block2.x) < tol &&
             abs(corner_block1.y  - corner_block2.y) < tol) {
             blocks_[i].x[ind_vec1[ind1]] =
               blocks_[j].x[ind_vec2[ind2]];
             blocks_[i].y[ind_vec1[ind1]] =
               blocks_[j].y[ind_vec2[ind2]];
          }
          ind2 += 1;
        }
        ind1 += 1;
      }
    }
  }
}

// std::vector of pair<Side,Side> of all interfaces that 
// have flipped orientation.
void MbGrid::SetFlippedInterfaces() {
  //('s','e'), ('s','s'),
  //('e','s'), ('e','e'),
  //('n','w'), ('n','n'),
  //('w','n'), ('w','w')]:
  flipped_interfaces_.push_back(std::make_pair(Side::s,Side::e));
  flipped_interfaces_.push_back(std::make_pair(Side::s,Side::s));
  flipped_interfaces_.push_back(std::make_pair(Side::e,Side::s));
  flipped_interfaces_.push_back(std::make_pair(Side::e,Side::e));
  flipped_interfaces_.push_back(std::make_pair(Side::n,Side::w));
  flipped_interfaces_.push_back(std::make_pair(Side::n,Side::n));
  flipped_interfaces_.push_back(std::make_pair(Side::w,Side::n));
  flipped_interfaces_.push_back(std::make_pair(Side::w,Side::w));
}

void ExportToTec(const MbGrid& grid, const MbArray& vec,
                 const std::string name) {

  std::vector<Block> blocks = grid.blocks();

  std::ofstream outfile;
  outfile.open(name + ".tec");
  outfile << "TITLE = 'ins_solution.tec' \n";
  outfile << "VARIABLES = x,y,u \n";

  Array X,Y;
  Shape shape;

  std::string str;
  int Nx, Ny, index;

  for(int k = 0; k < blocks.size(); ++k) {
    X = blocks[k].x;
    Y = blocks[k].y;
    shape = GetShape(X);
    Nx = shape.Nx;
    Ny = shape.Ny;

    str = "ZONE I = " + std::to_string(Ny) + ", J = "
                      + std::to_string(Nx) + ", F = POINT\n";
    outfile <<  str;

    for(int i = 0; i < X.size(); i++)
      outfile << std::to_string(X[i]) + " " +
                 std::to_string(Y[i]) +  " " +
                 std::to_string(vec[k][i]) + "\n";
  }
  outfile.close();
}
