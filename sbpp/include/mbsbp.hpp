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
 * the differentiation (Dx, Dy).
 *
 * Member functions:
 *  o
 *  o
 *  o
 *
 * Private functions:
 *  o
 *
 * Member variables:
 *  o
 *
 * Non-member functions:
 *  o
 *
 * Other data structes defined here:
 *  o
 *
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
//#include "interpolation.hpp"
//#include "interpolation21.hpp"
//#include "interpolation42.hpp"

// This is used in InterfaceTreatment
enum class Direction {x, y};

class MbSbp {
  public:
    MbSbp(const std::vector<Block>& blocks, int order);
    ~MbSbp(){};
    //MbSbp(){};

    ArrayPair GetNormals(int block_idx, Side side);

    void InterfaceTreatment(const MbArray& f,MbArray& df,
                            Direction direction);

    MbArray Dx(const MbArray& f);
    MbArray Dy(const MbArray& f);

    MbArray DxAndInterface(const MbArray& f);
    MbArray DyAndInterface(const MbArray& f);

    double Integrate(const MbArray&f);

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

    std::pair<int,std::slice>
      GetBlockBoundarySliceAndSize(int block_idx,
                             Side side);

    bool IsFilppedInterface(int interface_idx);

    // ------------------ MbGrid-functions -------------------

  private:
    // ------------------- functions ---------------------------
    void SetNormalsAndBoundaryQuadratures(int block_idx);
    void SetInterpolationOperators();

    // ------------------- variables ---------------------------
    MbGrid grid_;
    int order_ = 0;
    std::vector<std::unique_ptr<Sbp>> sbp_;
//    std::vector<std::unique_ptr<Interpolation>> interp_;

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
