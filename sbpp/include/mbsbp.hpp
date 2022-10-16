//============================================================
// Name        : MbSbp.hpp
// Date        : October 2022
// Author      : Fredrik Laurén
// Description : hpp file for the class MbSbp.
//============================================================

/* MbSbp (MbSbp).
 *
 * This class combines MbGrid (MultiblockGrid) and Sbp.
 *
 * The class takes as input a std::vector<Block>
 * and forms a MultiblockGrid.
 * The class also creates an Sbp-object on each Block
 * in the MbGrid. Interfaces are handled automatically during
 * differentiation (Dx, Dy).
 *
 * Member functions:
 *  o ArrayPair GetNormals(block_idx, side)
 *  o MbArray Dx(f)
 *  o MbArray Dy(f)
 *  o MbArray DxAndInterface(f)
 *  o MbArray DyAndInterface(f)
 *  o Integrate(f)
 *  o Array GetBoundaryQuadrature(block_idx, side)
 *  o Array GetPinvAtBoundary(block_idx, side)
 *  o MbGrid grid()
 * functions to be used for SAT:
 *  o Array Dn(f, block_idx, side)
 *  o Array Dx(f, block_idx)
 *  o Array Dy(f, block_idx)
 *  o Array DnT(f, block_idx, side)
 *  o Array DxT(f, block_idx)
 *  o Array DyT(f, block_idx)
 * MbGrid-functions:
 *  o MbArray Evaluate(function<double(double,double)> f)
 *  o std::vector<Block> blocks()
 *  o Block blocks(int block_idx)
 *  o int num_blocks()
 *  o std::vector<Shape> shapes()
 *  o Shape shapes(block_idx)
 *  o std::vector<Interface> interfaces()
 *  o std::vector<BoundaryInfo> boundaries()
 *  o Array ToBlockBoundary(f, block_idx, side)
 *  o Array ToBlockBoundary(f, block_idx, side)
 *  o std::pair<int,std::slice> GetBlockBoundarySliceAndSize(block_idx, Side side)
 *  o bool IsFilppedInterface(interface_idx)
 *
 * Private functions:
 *  o SetNormalsAndBoundaryQuadratures(block_idx)
 *  o SetInterpolationOperators()
 *  o InterfaceTreatment(f, df, direction)
 *
 * Member variables:
 *  o MbGrid grid
 *  o int order
 *  o vector<std::unique_ptr<Sbp>> sbp
 *  o vector<std::unique_ptr<Interp>> interp
 *  o vector<Array> x_xi, x_eta, y_xi, y_eta, J, invJ
 *  o p_inv_west_, p_inv_east, p_inv_south, p_inv_north
 *  o bd_quad_west_, bd_quad_east_, bd_quad_south_, bd_quad_north
 *  o nw, ne, ns, nn
 *
 * Other data structes defined here:
 *  o enum class Direction {x,y}
 */

#pragma once
#include <iostream>
#include <vector>
#include <fstream>
#include <math.h>
#include <memory> // for unique_ptr
#include "mbgrid.hpp"
#include "sbp.hpp"
#include "sbp21.hpp"
#include "sbp42.hpp"
#include "interp.hpp"
#include "interp21.hpp"
#include "interp42.hpp"

// This is used in InterfaceTreatment
enum class Direction {x, y};

class MbSbp {
  public:
    MbSbp(const std::vector<Block>& blocks, int order);
    ~MbSbp(){};

    ArrayPair GetNormals(int block_idx, Side side);

    MbArray Dx(const MbArray& f);
    MbArray Dy(const MbArray& f);

    MbArray DxAndInterface(const MbArray& f);
    MbArray DyAndInterface(const MbArray& f);

    double Integrate(const MbArray& f);

    Array GetBoundaryQuadrature(int block_idx, Side side);
    Array GetPinvAtBoundary(int block_idx, Side side);

    MbGrid grid();

    // functions to be used for SAT
    Array Dn(const MbArray& f, int block_idx, Side side);
    Array Dx(const MbArray& f, int block_idx);
    Array Dy(const MbArray& f, int block_idx);

    Array DnT(const MbArray& f, int block_idx, Side side);
    Array DxT(const MbArray& f, int block_idx);
    Array DyT(const MbArray& f, int block_idx);

    // ------------ MbGrid-functions --------------------
    MbArray Evaluate(std::function<double(double,double)> f);
    std::vector<Block> blocks();
    Block blocks(int block_idx);
    int num_blocks();

    std::vector<Shape> shapes();

    Shape shapes(int block_idx);

    std::vector<Interface> interfaces();
    std::vector<BoundaryInfo> boundaries();

    Array ToBlockBoundary(const MbArray& f,
                          int block_idx, Side side);

    Array ToBlockBoundary(const Array& f, int block_idx, Side side);

   BdSlice GetBdSlice(int block_idx, Side side);

    bool IsFilppedInterface(int interface_idx);


  private:
    // ------------------- functions ---------------------------
    void SetNormalsAndBoundaryQuadratures(int block_idx);
    void SetInterpolationOperators();
    void InterfaceTreatment(const MbArray& f, MbArray& df,
                            const Direction& direction);

    // ------------------- variables ---------------------------
    MbGrid grid_;
    int order_ = 0;
    std::vector<std::unique_ptr<Sbp>> sbp_;
    std::vector<std::unique_ptr<Interp>> interp_;

    // metric terms
    std::vector<Array> x_xi_;
    std::vector<Array> x_eta_;
    std::vector<Array> y_xi_;
    std::vector<Array> y_eta_;
    std::vector<Array> J_;
    std::vector<Array> invJ_;

    // p_inv at boundaries for all blocks (for SATs)
    // Ex. p_inv_west_[i] is p_inv at west boundary of the
    // i:th block.
    std::vector<Array> p_inv_west_;
    std::vector<Array> p_inv_east_;
    std::vector<Array> p_inv_south_;
    std::vector<Array> p_inv_north_;

    // boundary quadratures for each block
    // Ex. bd_quad_west_[i] is the west boundary
    // quadratue of the i:th block.
    std::vector<Array> bd_quad_west_;
    std::vector<Array> bd_quad_east_;
    std::vector<Array> bd_quad_south_;
    std::vector<Array> bd_quad_north_;

    // std::vectors containing the normals of each block.
    // Ex. nw_[k] is the Array-pair [n1, n2]
    // containing the normal at each
    // point of the west normal in the k:th block.
    std::vector<ArrayPair> nw_;
    std::vector<ArrayPair> ne_;
    std::vector<ArrayPair> ns_;
    std::vector<ArrayPair> nn_;

};
