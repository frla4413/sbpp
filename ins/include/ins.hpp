/*
 * Ins.h
 *
 *  Created on: 11 jan. 2021
 *      Author: frela05
 *
 *  Class for the problem I_tilde w_t + SpatialOperator(w) + Sat(w) = 0
 *  A "split" is performed on the spatial_operator for stability reasons.
 *  w - u,v,p. Struct of MbArrays (denoted InsState)
 *
 */

#pragma once
#include <iostream>
#include <valarray>
#include "array.hpp"
#include "mesh.hpp"
#include "mbarray.hpp"
#include "mbgrid.hpp"
#include "mbsbp.hpp"

struct InsState {
  MbArray u, v, p;
};

enum class BdType {Wall, Inflow, Outflow, Pressure};//, SlipWall};

class Ins {
  using Blocks = std::vector<Block>;
  using BdTypes = std::vector<BdType>;
  using valarray = std::valarray<double>;
  public:
    Ins(const Blocks& blocks, int order,
        const BdTypes& bd_types, double mu);
    valarray SpatialOperator(double t, const InsState& w);

    MbArray Evaluate(std::function<double(double,double)> f);

// -------------------- Sat -------------------------------------
    void ApplySat(double t, const InsState& w,
                  MbArray& l1, MbArray& l2, MbArray& l3);

    void WallSat(double t, const InsState& w,
                 MbArray& l1, MbArray& l2, MbArray& l3,
                 int block_idx, Side side);

    //void SlipWallSat(double t, const InsState& w,
    //                 MbArray& l1, MbArray& l2, MbArray& l3,
    //                 int block_idx, Side side);

    void InflowSat(double t, const InsState& w,
                   MbArray& l1, MbArray& l2, MbArray& l3,
                   int block_idx, Side side,
                   double wn_data, double wt_data);

    void PressureSat(double t, const InsState& w,
                     MbArray& l1, MbArray& l2, MbArray& l3,
                     int block_idx, Side side, double p_data);

    void OutflowSat(double t, const InsState& w,
                    MbArray& l1, MbArray& l2, MbArray& l3,
                    int block_idx, Side side);

// -------------------------------------------------------------

    valarray InsStateToValArray(const InsState& state);

    InsState ValArrayToInsState(const valarray& state);

    valarray MbArraysToValArray(const MbArray& u,
                                const MbArray& v,
                                const MbArray& p);

    void ExportToTec(const InsState& w, const std::string name);
    void ExportBoundaryToTec(const InsState& w,
                             const std::string& name,
                             int block_idx, const Side side);

// ---------------------- Jacobian ------------------------------
    valarray Jacobian(const InsState& f);

    void JacobianSat(const InsState& f, MbArray& J_l1_f, 
                     MbArray& J_l2_f, MbArray& J_l3_f);

    void JacobianWallSat(const InsState& f, MbArray& J_l1_f, 
                MbArray& J_l2_f, MbArray& J_l3_f, 
                int block_idx, Side side);

    //void JacobianSlipWallSat(const InsState& f,
    //                        MbArray& J_l1_f,
    //                        MbArray& J_l2_f,
    //                        MbArray& J_l3_f,
    //                        int block_idx, Side side);

    void JacobianInflowSat(const InsState& f, MbArray& J_l1_f,
                           MbArray& J_l2_f, MbArray& J_l3_f,
                           int block_idx, Side side,
                           double wn_data, double wt_data);

    void JacobianPressureSat(const InsState& f, MbArray& J_l1_f,
                             MbArray& J_l2_f, MbArray& J_l3_f,
                             int block_idx, Side side);

    void JacobianOutflowSat(const InsState& f, MbArray& J_l1_f,
                           MbArray& J_l2_f, MbArray& J_l3_f,
                           int block_idx, Side side);

  private:
   double mu_ = 0; // diffusion coefficient
   std::unique_ptr<MbSbp> sbp_;
   std::vector<BdType> bd_types_;

// Current state_ and gradients
// (used when evaluating the Jacobian)
   InsState state_;
   MbArray dudx_, dudy_, dvdx_, dvdy_;
};

/*
 * Set initial data. Edit in the cpp file
 */
InsState InitialData(Ins& ins);

/*
 *
 */
