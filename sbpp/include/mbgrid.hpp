//============================================================
// Name        : MbGrid.hpp
// Date        : September 2022
// Author      : Fredrik Laurén
// Description : hpp file for the class multiblockgrid.
//============================================================

/* MbGrid (MultiblockGrid).
 *
 * This class takes as input a std::vector<Block> and forms
 * a multi-block grid.
 *
 * Each block has four Corners.
 *
 * The block has the following structure:
 *
 *     --------------
 *    |      n       |
 *    |              |
 *    |w           e |
 *    |              |
 *    |      s       |
 *     --------------
 *   (0,0)
 *
 * The first grid point is in the s/w corner and thereafter 
 * orginized columnwise.
 *
 * Member functions:
 *  o blocks(), blocks(idx)
 *  o shapes(), shapes(idx)
 *  o int num_blocks()
 *  o GetTotalSize()
 *  o EvaluateFunction(function<double(double,double)>)
 *  o <int, slice> GetBlockBoundarySliceAndSize(idx, side)
 *  o ToBlockBoundary(MbArray& f, idx, side)
 *  o ArrayPair(Nx, Ny) GetNormals(idx, side)
 *  o interfaces()
 *  o vector<BoundaryInfo> boundaries()
 *  o Block GetBoundary(idx, side)
 *  o IsInterface(idx, side)
 *  o IsFilppedInterface(idx)
 *
 * Private functions:
 *  o void SetGrid() - find adjacent blocks, interfaces,
 *                       boundaries, ...
 *  o CollocateCorners()
 *  o SetNormals()
 *  o SaveUniqueCorners()
 *  o InterfacesAsCorners()
 *  o SaveUniqueEdges()
 *  o SaveFaceEdges()
 *  o FindInterfaces()
 *  o FindPhysicalBoundaries()
 *  o SetFlippedInterfaces()
 *  o static CompareEdgeSize(lhs, rhs)
 *  o static IsSameEdge(lhs, rhs)
 *
 * Member variables:
 *  o blocks - vector<Block>
 *  o num_blocks - int
 *  o shapes - vector<Shape>
 *  o corners - vector<Corner>
 *  o faces - vector<vector<int>>
 *  o edges - vector<pair<int,int>>
 *  o face_edges - vector<vector<pair<Side,int>>>
 *  o block_interfaces - vector<vector<BlockInterface>>
 *  o intrefaces vector<Interface>
 *  o boundaries - vector<BoundaryInfo>
 *  o n_w, n_e, n_s, n_n - vector<ArrayPair>
 *  o flipped_interfaces - vector<pair<Side,Side>>
 *
 * Non-member functions:
 *  o ExportToTec(mbgrid, mbarray, name)
 *
 * Comments:
 *    o Normals are not included in this class
 *      since they require differentiation operators.
 *      Normals are included in MbSbp.
 *
 * Other data structes defined here:
 *  o Corner - x,y
 *    - Size() --> x*x + y*y
 *    - GetCorners(block)
 *    - CompareCornerSize(lhs, rhs)
 *    - IsSameCorner(lhs, rhs)
 *  o enum class Side {s, e, n, w}
 *  o BlockInterface - my_side, adjacent_block_idx,
 *                     adjacent_block_size
 *  o Interface - block_idx1, block_idx2, side1, side2
 *    - Equal(lhs, rhs)
 *    - Print(Interface)
 *    - Print(vector<Interface>)
 *  o BoundaryInfo - block_idx, side
 *  o ArrayPair - a1, a2
 */

#pragma once
#include <iostream>
#include <vector>
#include <fstream>
#include <math.h>
#include "mbarray.hpp"
#include "mesh.hpp"

struct Corner {
   double x, y;
   double Size() const { return x*x + y*y; };
};

std::vector<Corner> GetCorners(const Block& block);
bool CompareCornerSize(const Corner& lhs, const Corner& rhs);
bool IsSameCorner(const Corner& lhs, const Corner& rhs);

/*
 * Calss to determine the side of a quadrilateral
 */
enum class Side {s, e, n, w};

/*
 * f[block_idx][slice] returns f at the boundary
 */

struct BdSlice {
  Side side;
  int block_idx;
  std::slice slice; //slice.size() to get size
};

/*
 * Contains information on the intreface to an adjacent block.
 */
struct BlockInterface {
  Side my_side;
  int adjacent_block_idx;
  Side adjacent_block_side;
};

/*
 * Contains information on an interface between two blocks.
 */
struct Interface {
  int block_idx1, block_idx2;
  Side side1, side2;
};

/*
 * True if intf1 and inft2 hold the same information,
 * not necessarily in the same order.
 */
bool Equal(const Interface& lhs, const Interface& rhs);

/*
 * Print an interface to the terminal.
 */
void Print(const Interface&);

/*
 * Print a vector of interfaces to the terminal.
 */
void Print(const std::vector<Interface>&);

/*
 * Information about a specific boundary.
 */
struct BoundaryInfo {
  int block_idx;
  Side side;
};

/*
 * A pair of Arrays.
 */
struct ArrayPair {
  Array a1, a2;
};


class MbGrid {
  //alias to improve readability
  using Shapes = std::vector<Shape>;
  using Arrays = std::vector<Array>;
  using Blocks = std::vector<Block>;
  using Corners = std::vector<Corner>;

  public:
    MbGrid(const Blocks& blocks);

    Blocks& blocks();
    const Blocks& blocks() const;
    Block& blocks(int block_idx);

    Shapes& shapes();
    const Shapes& shapes() const;
    Shape& shapes(int block_idx);
    int num_blocks();
    int GetTotalSize() const;

    // Evaluate function z = f(x,y)
    MbArray Evaluate(std::function<double(double,double)>f);

    BdSlice GetBdSlice(int block_idx, Side side);

    Array ToBlockBoundary(const MbArray& f,
                          int block_idx, Side side);

    // f is an Array, that lives on block block_idx
    Array ToBlockBoundary(const Array& f, int block_idx,
                          Side side);

    ArrayPair GetNormals(int block_idx, Side side);
    const std::vector<Interface>& interfaces() const;
    std::vector<std::vector<BlockInterface>> GetBlockInterfaces();
    std::vector<BoundaryInfo> boundaries();
    Block GetBoundary(int block_idx, Side);

    bool IsInterface(int block_idx, Side side);
    bool IsFilppedInterface(int interface_idx);

  private:
    // ------------------- variables ---------------------------
    int num_blocks_ = -1;
    Shapes shapes_;

    Blocks blocks_;
    Corners corners_;
    std::vector<std::vector<int>> faces_;
    std::vector<std::pair<int,int>> edges_;
    std::vector<std::vector<std::pair<Side,int>>> face_edges_;

    // vector of block-interfaces.
    // Ex: block_interfaces[1] = (w: 0 e), (e: 2w)
    // means that block 1 west side is adjacent to block 0 
    // east side and
    // that block 1 east side is adjacent to block 2 west side.
    std::vector<std::vector<BlockInterface>> block_interfaces_;

    // vector of interfaces.
    // Ex: interfaces[0] = (0, e: 1 w)
    // means that interface 0 is between the east side of block 0
    // and west side of block 1.
    std::vector<Interface> interfaces_;

    // vector containing all physical boundaries.
    // Ex boundaries_[0] = (0, Side::w) means that west side
    // of block 0
    // is a physical boundary.
    std::vector<BoundaryInfo> boundaries_;

    // vectors containing the normals of each block.
    // Ex. n_w_[i] is the Array2D-pair [n1, n2] containing 
    // the normal at each
    // point of the west normal in the i:th block.
    std::vector<ArrayPair> n_w_;
    std::vector<ArrayPair> n_e_;
    std::vector<ArrayPair> n_s_;
    std::vector<ArrayPair> n_n_;

    std::vector<std::pair<Side,Side>> flipped_interfaces_;

    // ------------------- functions ---------------------------
    // these functions are used in the constuctor to find,
    // corners, edges,
    // faces, interfaces, physical boundaries.
    void SetGrid();
    void CollocateCorners();
    void SetNormals();
    void SaveUniqueCorners();
    void InterfacesAsCorners();
    void SaveUniqueEdges();
    void SaveFaceEdges();
    void FindInterfaces();
    void FindPhysicalBoundaries();
    void SetFlippedInterfaces();

    bool static CompareEdgeSize(const std::pair<int,int>& lhs,
                                const std::pair<int,int>& rhs);

    bool static IsSameEdge(const std::pair<int,int>& lhs,
                           const std::pair<int,int>& rhs);

};

/*
 * Export the grid and array to name.tec.
 * Can be used in paraview to view the solution.
 */
void ExportToTec(const MbGrid& grid, const MbArray& array,
                 const std::string name);

