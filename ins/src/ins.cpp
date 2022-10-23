/*
 * Ins.cpp
 *
 *  Created on: 16 jan. 2021
 *      Author: frela05
 */

#include "ins.hpp"
#include <stdexcept>
#include <iostream>
#include <fstream>
#include <string>

using Blocks = std::vector<Block>;
using BdTypes = std::vector<BdType>;
using valarray = std::valarray<double>;

Ins::Ins(const Blocks& blocks, int order,
         const BdTypes& bd_types, double mu) {

  sbp_ = std::make_unique<MbSbp>(blocks, order);
  if(sbp_->boundaries().size() != bd_types.size()) {
    std::string msg = "bd_types and number \
                       of boundaries do not agree!";
     throw std::invalid_argument(msg);
  }
  bd_types_ = bd_types;
  mu_ = mu;
}


MbArray Ins::Evaluate(std::function<double(double,double)> f) {

  return sbp_->Evaluate(f);
}

/*
 * Build the spatial_operator
 * I_tilde w_t + SpatialOperator(w) + Sat(w) = 0
 */
valarray Ins::SpatialOperator (double t, const InsState& state) {
  auto dudx {sbp_->DxAndInterface(state.u)};
  auto dudy {sbp_->DyAndInterface(state.u)};

  auto duudx {sbp_->DxAndInterface(state.u*state.u)};

  auto duvdy {sbp_->DyAndInterface(state.u*state.v)};

  auto dvdx {sbp_->DxAndInterface(state.v)};
  auto dvdy {sbp_->DyAndInterface(state.v)};

  auto dvvdy {sbp_->DyAndInterface(state.v*state.v)};
  auto duvdx {sbp_->DxAndInterface(state.u*state.v)};

  auto dpdx {sbp_->DxAndInterface(state.p)};
  auto dpdy {sbp_->DyAndInterface(state.p)};

  auto l1 {0.5*(state.u*dudx + state.v*dudy +
                duudx + duvdy) + dpdx};
  auto l2 {0.5*(state.u*dvdx + state.v*dvdy +
                dvvdy + duvdx) + dpdy};
  auto l3 {dudx + dvdy};
  ApplySat(t, state, l1, l2, l3);

  // save for evaluating Jacobian
  state_ = state;
  dudx_ = dudx;
  dudy_ = dudy;
  dvdx_ = dvdx;
  dvdy_ = dvdy;

  if (mu_ != 0) {
    l1 -= mu_*(sbp_->DxAndInterface(dudx) +
               sbp_->DyAndInterface(dudy));
    l2 -= mu_*(sbp_->DxAndInterface(dvdx) +
               sbp_->DyAndInterface(dvdy));
  }
  return MbArraysToValArray(l1, l2, l3);
}

/*
 * Add SATs to l1, l2, l3.
 */
void Ins::ApplySat(double t, const InsState& w,
                   MbArray& l1, MbArray& l2, MbArray& l3) {
   int idx = 0;
  for(auto& bd : sbp_->boundaries()) {
    switch(bd_types_[idx]) {
      case BdType::Wall:
        WallSat(t, w, l1, l2, l3, bd.block_idx, bd.side);
        break;
//      case BdType::SlipWall:
//        SlipWallSat(t, w, l1, l2, l3, bd.block_idx, bd.side);
//        break;
      case BdType::Inflow:
        InflowSat(t, w, l1, l2, l3, 
                  bd.block_idx, bd.side, -1, 0);
        break;
//      case BdType::Pressure:
//        PressureSat(t, w, l1, l2, l3, bd.block_idx, bd.side, 0);
//        break;
      case BdType::Outflow:
        OutflowSat(t, w, l1, l2, l3, bd.block_idx, bd.side);
        break;
    }
    ++idx;
  }
}

void Ins::WallSat(double t, const InsState& w,
                  MbArray& l1, MbArray& l2, MbArray& l3,
                  int block_idx, Side side) {

  auto normals {sbp_->GetNormals(block_idx, side)};
  auto nx {normals.a1};
  auto ny {normals.a2};

  auto bd_slice {sbp_->GetBdSlice(block_idx, side)};

  auto Pinv {sbp_->GetPinvAtBoundary(block_idx, side)};
  auto Pbd {sbp_->GetBoundaryQuadrature(block_idx, side)};
  auto lift {Pinv*Pbd};

  int size = bd_slice.slice.size();
  auto u_bd {Array(1,size, w.u[block_idx][bd_slice.slice])};
  auto v_bd {Array(1,size, w.v[block_idx][bd_slice.slice])};

  auto wn {nx*u_bd + ny*v_bd};

  auto pen {-0.5*lift*u_bd*wn};
  l1[block_idx][bd_slice.slice] += pen.array();
  pen = -0.5*lift*v_bd*wn;
  l2[block_idx][bd_slice.slice] += pen.array();
  pen = -1.0*lift*wn;
  l3[block_idx][bd_slice.slice] += pen.array();

  if(mu_ != 0) {
    int Nx = sbp_->shapes(block_idx).Nx;
    int Ny = sbp_->shapes(block_idx).Ny;
    size = Pbd.size();
    auto u_Pbd {Array(1,size,
                w.u[block_idx][bd_slice.slice]*Pbd.array())};
    pen = mu_*Pinv*sbp_->DnT(u_Pbd, block_idx, side);
    l1[block_idx][bd_slice.slice] += pen.array();

    auto v_Pbd {Array(1,size,
                w.v[block_idx][bd_slice.slice]*Pbd.array())};

    pen = mu_*Pinv*sbp_->DnT(v_Pbd, block_idx, side);
    l2[block_idx][bd_slice.slice] += pen.array();
  }
}

//void Ins::SlipWallSat(double t, const InsState& w,
//                    MbArray& l1, MbArray& l2,
//                    MbArray& l3,
//                    int block_idx, Side side) {
//  auto normals = sbp_->GetNormals(block_idx, side);
//  auto nx = normals.a1, ny = normals.a2;
//
//  auto bd_slice = sbp_->GetBdSlice(block_idx, side);
//
//  auto Pinv {sbp_->GetPinvAtBoundary(block_idx, side)};
//  auto Pbd {sbp_->GetBoundaryQuadrature(block_idx, side)};
//  auto lift {Pinv*Pbd};
//
//  auto u_bd = Array2D(1,size, w.u[block_idx][slice]);
//  Array2D v_bd = Array2D(1,size, w.v[block_idx][slice]);
//  
//  Array2D wn = nx*u_bd + ny*v_bd;
//  
//  Array2D pen = -0.5*lift*u_bd*wn;
//  l1[block_idx][slice] += pen.array();
//  pen = -0.5*lift*v_bd*wn;
//  l2[block_idx][slice] += pen.array();
//  pen = -1.0*lift*wn;
//  l3[block_idx][slice] += pen.array();
//
//  if(mu_ != 0)
//  {
//     MbArray n_du = MbArray(sbp_->GetNumBlocks());
//     n_du[block_idx] = sbp_->Dn(w.u, block_idx, side);
//     MbArray n_dv = MbArray(sbp_->GetNumBlocks());
//     n_dv[block_idx] = sbp_->Dn(w.v, block_idx, side);
//
//     MbArray tau_t = MbArray(sbp_->GetNumBlocks());
//
//     tau_t[block_idx] = nx*n_dv[block_idx] - ny*n_du[block_idx];
//
//     MbArray wn_Pbd_nx = MbArray(sbp_->GetNumBlocks());
//     wn_Pbd_nx[block_idx] = Array2D(sbp_->GetShape(block_idx).first,
//                                sbp_->GetShape(block_idx).second);
//     wn_Pbd_nx[block_idx][slice] = wn.array()
//                                             *Pbd.array()
//                                             *nx.array();
//
//     pen = mu_*Pinv*(sbp_->DnT(wn_Pbd_nx, block_idx, side) -
//                     ny*Pbd*tau_t[block_idx]);
//
//     l1[block_idx][slice] += pen.array();
//
//     MbArray wn_Pbd_ny = MbArray(sbp_->GetNumBlocks());
//     wn_Pbd_ny[block_idx] = Array2D(sbp_->GetShape(block_idx).first,
//                                sbp_->GetShape(block_idx).second);
//     wn_Pbd_nx[block_idx][slice] = wn.array()
//                                             *Pbd.array()
//                                             *ny.array();
//     pen = mu_*Pinv*(sbp_->DnT(wn_Pbd_ny, block_idx, side) +
//                    nx*Pbd*tau_t[block_idx]);
//     l2[block_idx][slice] += pen.array();
//  }
//}

void Ins::InflowSat(double t, const InsState& w,
                     MbArray& l1, MbArray& l2, MbArray& l3,
                     int block_idx, Side side,
                     double wn_data, double wt_data) {

  auto normals {sbp_->GetNormals(block_idx, side)};
  auto nx {normals.a1}, ny {normals.a2};

  auto bd_slice {sbp_->GetBdSlice(block_idx, side)};

  auto Pinv {sbp_->GetPinvAtBoundary(block_idx, side)};
  auto Pbd {sbp_->GetBoundaryQuadrature(block_idx, side)};

  int size = bd_slice.slice.size();
  auto slice {bd_slice.slice};
  auto u_bd {Array(1,size, w.u[block_idx][slice])};
  auto v_bd {Array(1,size, w.v[block_idx][slice])};

  auto wn {nx*u_bd + ny*v_bd};
  auto wt {-1.0*ny*u_bd + nx*v_bd};
  auto lift {Pinv*Pbd};

  auto ybd {sbp_->ToBlockBoundary(sbp_->blocks(block_idx).y,
                                   block_idx,side)};

  //Array2D u_data = Array2D(1,ybd.Size(),-1 + std::tanh(20*(ybd.array()+1))
  //                                         - std::tanh(20*(ybd.array()-1)));
  auto u_data {nx*wn_data - ny*wt_data};

  auto pen {-1.0*lift*wn*(u_bd - u_data)};
  l1[block_idx][slice] +=  pen.array();

  auto v_data {ny*wn_data + nx*wt_data};
  pen = -1.0*lift*wn*(v_bd - v_data);
  l2[block_idx][slice] += pen.array();

  pen = -1.0*lift*(wn-wn_data);
  l3[block_idx][slice] += pen.array();

  if(mu_ != 0) {
    auto u_Pbd {w.u[block_idx]};
    u_Pbd[slice] -= u_data.array();
    u_Pbd[slice] *= Pbd.array();
    pen = mu_*Pinv*sbp_->DnT(u_Pbd, block_idx, side);
    l1[block_idx][slice] += pen.array();

    auto v_Pbd {w.v[block_idx]};
    v_Pbd[slice] -= v_data.array();
    v_Pbd[slice] *= Pbd.array();
    pen = mu_*Pinv*sbp_->DnT(v_Pbd, block_idx, side);
    l2[block_idx][slice] += pen.array();
  }
}

//void Ins::PressureSat(double t, const InsState& w, 
//                        MbArray& l1, MbArray& l2, MbArray& l3, 
//                        int block_idx, Side side, double p_data)
//{
//
//   std::pair<Array2D, Array2D> normals = sbp_->GetNormals(block_idx, side);
//   Array2D nx = normals.first, ny = normals.second;
//   
//   auto bd_slice = sbp_->GetBlockBoundarySliceAndSize(block_idx, side);
//   
//   Array2D Pinv = sbp_->GetPinvAtBoundary(block_idx, side);
//   Array2D Pbd = sbp_->GetBoundaryQuadrature(block_idx, side);
//   Array2D lift = Pinv*Pbd;
//   
//   Array2D p_bd = Array2D(1,size, w.p[block_idx][slice]);
//   
//   Array2D pen = -1.0*lift*nx*(p_bd - p_data);
//   l1[block_idx][slice] +=  pen.array();
//   pen = -1.0*lift*ny*(p_bd-p_data);
//   l2[block_idx][slice] +=  pen.array();
//
//   if(mu_ != 0)
//   {
//      pen = mu_*lift*sbp_->Dn(w.u, block_idx, side);
//      l1[block_idx][slice] += pen.array();
//
//      pen = mu_*lift*sbp_->Dn(w.v, block_idx, side);
//      l2[block_idx][slice] += pen.array();
//   }
//}

void Ins::OutflowSat(double t, const InsState& w,
                    MbArray& l1, MbArray& l2, MbArray& l3,
                    int block_idx, Side side) {

  auto normals {sbp_->GetNormals(block_idx, side)};
  auto nx {normals.a1}, ny {normals.a2};

  auto bd_slice = sbp_->GetBdSlice(block_idx, side);

  auto Pinv = sbp_->GetPinvAtBoundary(block_idx, side);
  auto Pbd = sbp_->GetBoundaryQuadrature(block_idx, side);
  auto lift = Pinv*Pbd;

  int size = bd_slice.slice.size();
  auto slice {bd_slice.slice};
  auto u_bd {Array(1,size, w.u[block_idx][slice])};
  auto v_bd {Array(1,size, w.v[block_idx][slice])};
  auto p_bd {Array(1,size, w.p[block_idx][slice])};

  auto wn {nx*u_bd + ny*v_bd};

  double a = 30, c = 0;
  double mag = 2;
  auto h {mag*Array(1, size,
      exp(c)/(std::exp(a*wn.array()) + exp(c)))};

  auto w2 {u_bd*u_bd + v_bd*v_bd};

  auto pen {-1.0*lift*(nx*(h*w2+p_bd))};
  l1[block_idx][slice] +=  pen.array();
  pen =  -1.0*lift*(ny*(h*w2+p_bd));
  l2[block_idx][slice] +=  pen.array();

  if(mu_ != 0) {
    pen = mu_*lift*sbp_->Dn(w.u, block_idx, side);
    l1[block_idx][slice] += pen.array();

    pen = mu_*lift*sbp_->Dn(w.v, block_idx, side);
    l2[block_idx][slice] += pen.array();
  }
}

valarray Ins::InsStateToValArray(const InsState& w) {
  valarray u = w.u.ToValarray();
  valarray v = w.v.ToValarray();
  valarray p = w.p.ToValarray();

  valarray array(3*u.size());
  int start = 0, size = u.size();
  array[std::slice(start,size,1)] = u;
  start = u.size();
  array[std::slice(start,size,1)] = v;
  start = 2*u.size();
  array[std::slice(start,size,1)] = p;
  return array;
}

valarray Ins::MbArraysToValArray(const MbArray& u,
                                 const MbArray& v,
                                 const MbArray& p) {
   int start = 0, size = u.GetTotalSize();
   valarray array(3*size);

   // Concatenate L to one long std::valarray
   array[std::slice(start,size,1)] = u.ToValarray();
   start = size;
   array[std::slice(start,size,1)] = v.ToValarray();
   start = 2*size;
   array[std::slice(start,size,1)] = p.ToValarray();
   return array;
}

InsState Ins::ValArrayToInsState(const valarray& w) {

  int start = 0, size = w.size()/3;
  auto shapes {sbp_->shapes()};
  auto u = MbArray(shapes, w[std::slice(start,size,1)]);
  start += size;
  auto v = MbArray(shapes, w[std::slice(start,size,1)]);
  start += size;
  auto p = MbArray(shapes,w[std::slice(start,size,1)]);
  return {u,v,p};
}

void Ins::ExportToTec(const InsState& state,
                      const std::string name) {

  Blocks blocks = sbp_->blocks();

  std::ofstream outfile;
  outfile.open(name + ".tec");
  outfile << "TITLE = 'solution.tec' \n";
  outfile << "VARIABLES = x,y,u,v,p,vort \n";

  Array X,Y;

  std::string str;
  int Nx, Ny, index;

  MbArray vx = sbp_->DxAndInterface(state.v);
  MbArray uy = sbp_->DyAndInterface(state.u);

  MbArray vorticity = vx - uy;

  for(int k = 0; k < blocks.size(); k++) {
    X = blocks[k].x;
    Y = blocks[k].y;
    auto shape = GetShape(X);
    Nx = shape.Nx;
    Ny = shape.Ny;

    str = "ZONE I = " + std::to_string(Ny) + ", J = "
                      + std::to_string(Nx) + ", F = POINT\n";
    outfile <<  str;

    for(int i = 0; i < X.size(); ++i)
      outfile << std::to_string(X[i]) + " " +
                 std::to_string(Y[i]) +  " " +
                 std::to_string(state.u[k][i]) + " " +
                 std::to_string(state.v[k][i]) + " " +
                 std::to_string(state.p[k][i]) + " " +
                 std::to_string(vorticity[k][i]) + "\n";
  }
  outfile.close();
}

//void Ins::ExportBoundaryToTec(const InsState& state,
//                              const std::string name,
//                              int block_idx,
//                              const Side side) {
//
//   std::vector<DataTypes::Block> blocks = sbp_->GetBlocks();
//
//   std::ofstream outfile;
//   outfile.open(name + ".tec");
//   outfile << "Incompressible_navier_stokes_solution.txt\n";
//   outfile << "VARIABLES = x,y,u,v,p,uy\n";
//
//   auto uy = sbp_->DyAndInterface(state.u);
//
//   int Nx = state.u.GetShape(1).first;
//   int Ny = state.u.GetShape(1).second;
//   int start, stride, size;
//
//   auto block_x = blocks[block_idx].x;
//   auto block_y = blocks[block_idx].y;
//
//   auto x = sbp_->GetGrid().ToBlockBoundary(block_x,
//                                            block_idx, side);
//   auto y = sbp_->GetGrid().ToBlockBoundary(block_y,
//                                            block_idx, side);
//   auto u = sbp_->GetGrid().ToBlockBoundary(state.u,
//                                            block_idx, side);
//   auto v = sbp_->GetGrid().ToBlockBoundary(state.v,
//                                            block_idx, side);
//   auto p = sbp_->GetGrid().ToBlockBoundary(state.p,
//                                            block_idx, side);
//
//   auto uy_bd = sbp_->GetGrid().ToBlockBoundary(uy, block_idx,
//                                                side);
//
//   for(int i = 0; i < x.Size(); i++)
//      outfile << std::fixed << std::setprecision(12) <<
//                 x[i] << " " <<
//                 y[i] << " " <<
//                 u[i] << " " <<
//                 v[i] << " " <<
//                 p[i] << " " <<
//                 uy_bd[i] << "\n";
//   outfile.close();
//}
//
//void ExportYSlice(const DataTypes::InsState& state,
//                  MultiblockSbp* sbp,
//                  const std::string name, int block_idx, int idx) {
//
//   std::vector<DataTypes::Block> blocks = sbp->GetBlocks();
//
//   std::ofstream outfile;
//   outfile.open(name + ".tec");
//   outfile << "Incompressible_navier_stokes_solution.txt\n";
//   outfile << "VARIABLES = x,y,u,v,p,uy\n";
//
//   auto dudy = sbp->DyAndInterface(state.u);
//
//   int Nx = state.u.GetShape(1).first;
//   int Ny = state.u.GetShape(1).second;
//
//   auto block_x = blocks[block_idx].x;
//   auto block_y = blocks[block_idx].y;
//
//   auto x = sbp->GetBlock(block_idx).x;
//   auto y = sbp->GetBlock(block_idx).y;
//   auto u = state.u[block_idx];
//   auto v = state.u[block_idx];
//   auto p = state.u[block_idx];
//
//   auto uy = dudy[block_idx];
//
//   int start = Ny*40;
//   int stop = Ny*41;
//
//   for(int i = start; i < stop; i++)
//      outfile << std::fixed << std::setprecision(12) <<
//                 x[i] << " " <<
//                 y[i] << " " <<
//                 u[i] << " " <<
//                 v[i] << " " <<
//                 p[i] << " " <<
//                 uy[i] << "\n";
//   outfile.close();
//}

//---------------------------------- Jacobian ---------------------------------------
/*
 *  Jacobian for L(w) + S(w).
 *  The SAT-term is evaluated for g = 0 (no data).
 *  If a force is used in L(w),
 *  it is not included since it is independent of the
 *  InsState.
 *  The function returns J(state_)*f
 */
valarray Ins::Jacobian(const InsState& f) {
  auto dfudx = sbp_->DxAndInterface(f.u);
  auto dfudy = sbp_->DyAndInterface(f.u);

  auto J_l1_f = 0.5*(state_.u*dfudx + dudx_*f.u +
                     state_.v*dfudy +
                     sbp_->DyAndInterface(state_.v*f.u)) +
                     sbp_->DxAndInterface(state_.u*f.u);

  J_l1_f += 0.5*(dudy_*f.v + sbp_->DyAndInterface(state_.u*f.v));
  J_l1_f += sbp_->DxAndInterface(f.p);

  auto dfvdx = sbp_->DxAndInterface(f.v);
  auto dfvdy = sbp_->DyAndInterface(f.v);

  auto J_l2_f = 0.5*(dvdx_*f.u +
                     sbp_->DxAndInterface(state_.v*f.u));

  J_l2_f += 0.5*(state_.v*dfvdy + dvdy_*f.v + state_.u*dfvdx +
                 sbp_->DxAndInterface(state_.u*f.v)) +
                 sbp_->DyAndInterface(state_.v*f.v);
  J_l2_f += sbp_->DyAndInterface(f.p);

  auto J_l3_f = sbp_->DxAndInterface(f.u) +
                sbp_->DyAndInterface(f.v);

  if (mu_ != 0) {
    J_l1_f -= mu_*(sbp_->DxAndInterface(dfudx) +
                   sbp_->DyAndInterface(dfudy));
     J_l2_f -= mu_*(sbp_->DxAndInterface(dfvdx) +
                    sbp_->DyAndInterface(dfvdy));
  }
  JacobianSat(f, J_l1_f, J_l2_f, J_l3_f);
  return MbArraysToValArray(J_l1_f, J_l2_f, J_l3_f);
}

void Ins::JacobianSat(const InsState& f,  MbArray& J_l1_f,
                      MbArray& J_l2_f, MbArray& J_l3_f) {

  int bd_idx = 0;
  for(auto& bd : sbp_->boundaries()) {
    switch(bd_types_[bd_idx]) {
      case BdType::Wall:
        JacobianWallSat(f, J_l1_f, J_l2_f, J_l3_f,
                        bd.block_idx, bd.side);
        break;
       //case BdType::SlipWall:
       //      JacobianSlipWallSat(f, J_l1_f, J_l2_f, J_l3_f,
       //                          bd.block_idx, bd.side);
       //      break;
      case BdType::Inflow:
        JacobianInflowSat(f, J_l1_f, J_l2_f, J_l3_f,
                          bd.block_idx, bd.side, -1, 0);
        break;
    //   case BdType::Pressure:
    //         JacobianPressureSat(f, J_l1_f, J_l2_f, J_l3_f,
    //                             bd.block_idx, bd.side);
    //         break;
      case BdType::Outflow:
       JacobianOutflowSat(f, J_l1_f, J_l2_f, J_l3_f,
                          bd.block_idx, bd.side);
       break;
    }
    ++bd_idx ;
  }
}

void Ins::JacobianWallSat(const InsState& f, MbArray& J_l1_f,
                          MbArray& J_l2_f, MbArray& J_l3_f,
                          int block_idx, Side side) {

  auto normals {sbp_->GetNormals(block_idx, side)};
  auto nx {normals.a1}, ny {normals.a2};

  auto bd_slice {sbp_->GetBdSlice(block_idx, side)};

  auto Pinv {sbp_->GetPinvAtBoundary(block_idx, side)};
  auto Pbd {sbp_->GetBoundaryQuadrature(block_idx, side)};
  auto lift {Pinv*Pbd};

  int size = bd_slice.slice.size();
  auto slice {bd_slice.slice};
  auto wu_bd {Array(1,size,
              state_.u[block_idx][slice])};
  auto wv_bd {Array(1,size,
                       state_.v[block_idx][bd_slice.slice])};

  auto wwn = nx*wu_bd + ny*wv_bd;
  auto fu_bd {Array(1,size, f.u[block_idx][slice])};
  auto fv_bd {Array(1,size, f.v[block_idx][slice])};

  auto pen {-0.5*lift*((wu_bd*nx + wwn)*fu_bd +
                        wu_bd*ny*fv_bd)};
  J_l1_f[block_idx][slice] +=  pen.array();

  pen = -0.5*lift*(wv_bd*nx*fu_bd +
                  (wv_bd*ny + wwn)*fv_bd);
  J_l2_f[block_idx][slice] +=  pen.array();
  pen = -1.0*lift*(nx*fu_bd + ny*fv_bd);
  J_l3_f[block_idx][slice] +=  pen.array();

  if(mu_ != 0) {
    auto u_Pbd = f.u;
    u_Pbd[block_idx][slice] *= Pbd.array();
    pen = mu_*Pinv*sbp_->DnT(u_Pbd[block_idx], block_idx, side);
    J_l1_f[block_idx][slice] += pen.array();

    auto v_Pbd = f.v;
    v_Pbd[block_idx][slice] *= Pbd.array();
    pen = mu_*Pinv*sbp_->DnT(v_Pbd[block_idx], block_idx, side);
    J_l2_f[block_idx][slice] += pen.array();
  }
}

//void Ins::JacobianSlipWallSat(const InsState& f,
//                              MbArray& J_l1_f,
//                              MbArray& J_l2_f,
//                              MbArray& J_l3_f,
//                              int block_idx, Side side) {
//
//  std::pair<Array2D,Array2D> normals = sbp_->GetNormals(block_idx, side);
//  Array2D nx = normals.first, ny = normals.second;
//
//  auto bd_slice = sbp_->GetBlockBoundarySliceAndSize(block_idx, side);
//
//  Array2D Pinv = sbp_->GetPinvAtBoundary(block_idx, side);
//  Array2D Pbd = sbp_->GetBoundaryQuadrature(block_idx, side);
//  Array2D lift = Pinv*Pbd;
//
//  Array2D wu_bd = Array2D(1,size, state_.u[block_idx][slice]);
//  Array2D wv_bd = Array2D(1,size, state_.v[block_idx][slice]);
//
//  Array2D wwn = nx*wu_bd + ny*wv_bd;
//  Array2D fu_bd = Array2D(1,size, f.u[block_idx][slice]);
//  Array2D fv_bd = Array2D(1,size, f.v[block_idx][slice]);
//
//  Array2D pen = -0.5*lift*((wu_bd*nx + wwn)*fu_bd + 
//                                 wu_bd*ny*fv_bd);
//  J_l1_f[block_idx][slice] +=  pen.array();
//
//  pen = -0.5*lift*(wv_bd*nx*fu_bd + 
//                         (wv_bd*ny + wwn)*fv_bd);
//  J_l2_f[block_idx][slice] +=  pen.array();
//  pen = -1.0*lift*(nx*fu_bd + ny*fv_bd);
//  J_l3_f[block_idx][slice] +=  pen.array();
//
//  if(mu_ != 0)
//  {
//     MbArray n_du = MbArray(sbp_->GetNumBlocks());
//     n_du[block_idx] = sbp_->Dn(f.u, block_idx, side);
//     MbArray n_dv = MbArray(sbp_->GetNumBlocks());
//     n_dv[block_idx] = sbp_->Dn(f.v, block_idx, side);
//
//     Array2D f_wn = nx*fu_bd + ny*fv_bd;
//
//     MbArray tau_t = MbArray(sbp_->GetNumBlocks());
//
//     tau_t[block_idx] = nx*n_dv[block_idx] - ny*n_du[block_idx];
//
//     MbArray fwn_Pbd_nx = MbArray(sbp_->GetNumBlocks());
//     fwn_Pbd_nx[block_idx] = Array2D(sbp_->GetShape(block_idx).first,
//                                sbp_->GetShape(block_idx).second);
//     fwn_Pbd_nx[block_idx][slice] = f_wn.array()
//                                              *Pbd.array()
//                                              *nx.array();
//
//     pen = mu_*Pinv*(sbp_->DnT(fwn_Pbd_nx, block_idx, side) -
//                     ny*Pbd*tau_t[block_idx]);
//
//     J_l1_f[block_idx][slice] += pen.array();
//
//     MbArray fwn_Pbd_ny = MbArray(sbp_->GetNumBlocks());
//     fwn_Pbd_ny[block_idx] = Array2D(sbp_->GetShape(block_idx).first,
//                                sbp_->GetShape(block_idx).second);
//     fwn_Pbd_nx[block_idx][slice] = f_wn.array()
//                                              *Pbd.array()
//                                              *ny.array();
//     pen = mu_*Pinv*(sbp_->DnT(fwn_Pbd_ny, block_idx, side) +
//                    nx*Pbd*tau_t[block_idx]);
//     J_l2_f[block_idx][slice] += pen.array();
//  }
//}

void Ins::JacobianInflowSat(const InsState& f, MbArray& J_l1_f,
      MbArray& J_l2_f, MbArray& J_l3_f, int block_idx, Side side,
      double wn_data, double wt_data) {

  auto normals {sbp_->GetNormals(block_idx, side)};
  auto nx {normals.a1}, ny {normals.a2};

  auto bd_slice {sbp_->GetBdSlice(block_idx, side)};

  auto Pinv {sbp_->GetPinvAtBoundary(block_idx, side)};
  auto Pbd {sbp_->GetBoundaryQuadrature(block_idx, side)};
  auto lift {Pinv*Pbd};

  int size = bd_slice.slice.size();
  auto slice {bd_slice.slice};

  auto wu_bd {Array(1,size, state_.u[block_idx][slice])};
  auto wv_bd {Array(1,size, state_.v[block_idx][slice])};

  auto wwn {nx*wu_bd + ny*wv_bd};
  auto wwt {-1.0*ny*wu_bd + nx*wv_bd};
  auto fu_bd {Array(1,size, f.u[block_idx][slice])};
  auto fv_bd {Array(1,size, f.v[block_idx][slice])};

  // u
  auto u_data {nx*wn_data - ny*wt_data};
  auto pen {-1.0*lift*((wwn + nx*(wu_bd - u_data))*fu_bd +
            ny*(wu_bd-u_data)*fv_bd)};

  J_l1_f[block_idx][slice] +=  pen.array();
  Array v_data {ny*wn_data + nx*wt_data};

  // v
  pen = -1.0*lift*(nx*(wv_bd - v_data)*fu_bd +
                  (wwn + ny*(wv_bd - v_data))*fv_bd);
  J_l2_f[block_idx][slice] +=  pen.array();

  // p
  pen = -1.0*lift*(nx*fu_bd + ny*fv_bd);
  J_l3_f[block_idx][slice] +=  pen.array();

  if(mu_ != 0) {
    auto u_Pbd {f.u[block_idx]};
    u_Pbd[slice] *= Pbd.array();
    pen = mu_*Pinv*sbp_->DnT(u_Pbd, block_idx, side);
    J_l1_f[block_idx][slice] += pen.array();

    auto v_Pbd {f.v[block_idx]};
    v_Pbd[slice] *= Pbd.array();
    pen = mu_*Pinv*sbp_->DnT(v_Pbd, block_idx, side);
    J_l2_f[block_idx][slice] += pen.array();
   }
}

//void Ins::JacobianPressureSat(const InsState& f, MbArray& J_l1_f, 
//      MbArray& J_l2_f, MbArray& J_l3_f, int block_idx, Side side)
//{
//
//   std::pair<Array2D,Array2D> normals = sbp_->GetNormals(block_idx, side);
//   Array2D nx = normals.first, ny = normals.second;
//
//   auto bd_slice = sbp_->GetBlockBoundarySliceAndSize(block_idx, side);
//
//   Array2D Pinv = sbp_->GetPinvAtBoundary(block_idx, side);
//   Array2D Pbd = sbp_->GetBoundaryQuadrature(block_idx, side);
//   Array2D lift = Pinv*Pbd;
//
//   Array2D fp_bd = Array2D(1,size, f.p[block_idx][slice]);
//
//   Array2D pen = -1.0*lift*nx*fp_bd;
//   J_l1_f[block_idx][slice] +=  pen.array();
//
//   pen = -1.0*lift*ny*fp_bd;
//   J_l2_f[block_idx][slice] +=  pen.array();
//   if(mu_ != 0)
//   {
//      pen = mu_*lift*sbp_->Dn(f.u, block_idx, side);
//      J_l1_f[block_idx][slice] += pen.array();
//
//      pen = mu_*lift*sbp_->Dn(f.v, block_idx, side);
//      J_l2_f[block_idx][slice] += pen.array();
//   }
//}
//
void Ins::JacobianOutflowSat(const InsState& f, MbArray& J_l1_f,
      MbArray& J_l2_f, MbArray& J_l3_f,
      int block_idx, Side side) {

  auto normals {sbp_->GetNormals(block_idx, side)};
  auto nx {normals.a1}, ny {normals.a2};

  auto bd_slice {sbp_->GetBdSlice(block_idx, side)};

  auto Pinv {sbp_->GetPinvAtBoundary(block_idx, side)};
  auto Pbd {sbp_->GetBoundaryQuadrature(block_idx, side)};
  auto lift {Pinv*Pbd};

  int size = bd_slice.slice.size();
  auto slice {bd_slice.slice};

  auto u_bd {Array(1,size, state_.u[block_idx][slice])};
  auto v_bd {Array(1,size, state_.v[block_idx][slice])};
  auto p_bd {Array(1,size, state_.p[block_idx][slice])};

  auto wn {nx*u_bd + ny*v_bd};
  double a = 30, c = 0, mag = 2;;
  auto h {mag*Array(1, size,
              exp(c)/(std::exp(a*wn.array()) + exp(c)))};

  Array w2 {u_bd*u_bd + v_bd*v_bd};

  valarray tmp = -a*exp(a*wn.array()+c)/
                       (pow((exp(a*wn.array())+ exp(c)),2));

  auto dhdwn {mag*Array(1,size,tmp)};
  auto dhdu  {dhdwn*nx};
  auto dhdv  {dhdwn*ny};

  auto dw2du {2*u_bd};
  auto dw2dv {2*v_bd};

  auto fu_bd {Array(1,size, f.u[block_idx][slice])};
  auto fv_bd {Array(1,size, f.v[block_idx][slice])};
  auto fp_bd {Array(1,size, f.p[block_idx][slice])};

  auto pen {-1.0*lift*nx*((dhdu*w2+h*dw2du)*fu_bd +
                           (dhdv*w2+h*dw2dv)*fv_bd + fp_bd)};
  J_l1_f[block_idx][slice] +=  pen.array();

  pen = -1.0*lift*ny*((dhdu*w2+h*dw2du)*fu_bd +
                      (dhdv*w2+h*dw2dv)*fv_bd +
                      fp_bd);

  J_l2_f[block_idx][slice] +=  pen.array();
  if(mu_ != 0) {
    pen = mu_*lift*sbp_->Dn(f.u, block_idx, side);
    J_l1_f[block_idx][slice] += pen.array();

    pen = mu_*lift*sbp_->Dn(f.v, block_idx, side);
    J_l2_f[block_idx][slice] += pen.array();
  }
}
