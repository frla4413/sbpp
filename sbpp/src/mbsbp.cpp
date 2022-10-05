/*
 * multiblockgrid.cpp
 *
 *  Created   : August 2021
 *      Author: Fredrik Laurén
 */

#include "mbsbp.hpp"


MbSbp::MbSbp(const std::vector<Block>& blocks, int order) :
  grid_(MbGrid(blocks)), order_ (order) {

  int block_idx = 0, idx = 0;

  for (auto& block : blocks) {
    Shape shape = grid_.shapes(block_idx);
    int Nx = shape.Nx;
    int Ny = shape.Ny;

    double d_xi = 1.0/static_cast<double>(Nx-1);
    double d_eta = 1.0/static_cast<double>(Ny-1);

    switch (order_) {
      case 4:
        sbp_.push_back(std::make_unique<Sbp42>(Nx,Ny,d_xi,d_eta));
        break;
      default:
      sbp_.push_back(std::make_unique<Sbp21>(Nx,Ny,d_xi,d_eta));
    }
    x_xi_.push_back(sbp_[block_idx]->DXi(block.x));
    x_eta_.push_back(sbp_[block_idx]->DEta(block.x));
    y_xi_.push_back(sbp_[block_idx]->DXi(block.y));
    y_eta_.push_back(sbp_[block_idx]->DEta(block.y));
    J_.push_back(x_xi_[block_idx]*y_eta_[block_idx] -
                 x_eta_[block_idx]*y_xi_[block_idx]);
    invJ_.push_back(1.0/J_[block_idx]);
    SetNormalsAndBoundaryQuadratures(block_idx);
    ++block_idx;
    ++idx;
  }
//  SetInterpolationOperators();
}



int MbSbp::num_blocks() {
  return grid_.num_blocks();
}

MbGrid MbSbp::grid() {
  return grid_;
};

std::vector<Interface> MbSbp::interfaces() {
  return grid_.interfaces();
};

std::vector<BoundaryInfo> MbSbp::boundaries() {
   return grid_.boundaries();
};

Array MbSbp::GetPinvAtBoundary(int block_idx, Side side) {
  switch (side) {
    case Side::w: {return p_inv_west_[block_idx];}
    case Side::e: {return p_inv_east_[block_idx];}
    case Side::s: {return p_inv_south_[block_idx];}
    case Side::n: {return p_inv_north_[block_idx];}
  }
  return Array(1,1,{-1});
}

Array MbSbp::GetBoundaryQuadrature(int block_idx, Side side) {
  switch (side) {
    case Side::w: {return bd_quad_west_[block_idx];}
    case Side::e: {return bd_quad_east_[block_idx];}
    case Side::s: {return bd_quad_south_[block_idx];}
    case Side::n: {return bd_quad_north_[block_idx];}
  }
  return Array(1,1,{-1});
}

MbArray MbSbp::Dx(const MbArray& f) {

   std::vector<Array> df_dx(num_blocks());
   for (int i = 0; i < num_blocks(); ++i) {
     Array temp = y_eta_[i]*f[i];
     Array df = y_eta_[i]*sbp_[i]->DXi(f[i]) +
                sbp_[i]->DXi(temp);
     temp = y_xi_[i]*f[i];
     df = df - y_xi_[i]*sbp_[i]->DEta(f[i]) -
               sbp_[i]->DEta(temp);
     df_dx[i] = 0.5*invJ_[i]*df;
   }
   return {df_dx};
}


MbArray MbSbp::Dy(const MbArray& f) {

   std::vector<Array> df_dy(num_blocks());;

   for (int i = 0; i < num_blocks(); ++i) {
     Array temp = x_xi_[i]*f[i];
     Array df = x_xi_[i]*sbp_[i]->DEta(f[i]) +
                sbp_[i]->DEta(temp);
     temp = x_eta_[i]*f[i];
     df = df - x_eta_[i]*sbp_[i]->DXi(f[i]) -
          sbp_[i]->DXi(temp);
     df_dy[i] = 0.5*invJ_[i]*df;
   }
   return {df_dy};
}

double MbSbp::Integrate(const MbArray& f) {
  // integrate a multiblock function
  double integral = 0;
  for (int i = 0; i < num_blocks(); ++i) {
    integral += sbp_[i]->Integrate(f[i]*J_[i]);
  }
  return integral;
}

ArrayPair MbSbp::GetNormals(int block_idx, Side side) {

  ArrayPair normals;
  switch (side) {
    case Side::w: { normals = nw_[block_idx]; break; }
    case Side::e: { normals = ne_[block_idx]; break; }
    case Side::s: { normals = ns_[block_idx]; break; }
    case Side::n: { normals = nn_[block_idx]; break; }
  }
  return normals;
}

Array MbSbp::Dn(const MbArray& f, int block_idx, Side side) {
  auto normals = GetNormals(block_idx, side);
  auto bd_slice = GetBlockBoundarySliceAndSize(block_idx, side);

  std::valarray<double> dfdx = Dx(f, block_idx)[bd_slice.second];
  std::valarray<double> dfdy = Dy(f, block_idx)[bd_slice.second];

  return Array(1, bd_slice.first, normals.a1.array()*dfdx +
                                  normals.a2.array()*dfdy);
}

Array MbSbp::Dx(const MbArray& f, int block_idx) {

  Array temp = y_eta_[block_idx]*f[block_idx];
  Array df = y_eta_[block_idx]*
             sbp_[block_idx]->DXi(f[block_idx]) +
             sbp_[block_idx]->DXi(temp);
  temp = y_xi_[block_idx]*f[block_idx];
  df = df - y_xi_[block_idx]*
       sbp_[block_idx]->DEta(f[block_idx]) -
       sbp_[block_idx]->DEta(temp);
  return 0.5*invJ_[block_idx]*df;
}

Array MbSbp::Dy(const MbArray& f, int block_idx) {

  Array temp = x_xi_[block_idx]*f[block_idx];
  Array df = x_xi_[block_idx]*
             sbp_[block_idx]->DEta(f[block_idx]) +
             sbp_[block_idx]->DEta(temp);
  temp = x_eta_[block_idx]*f[block_idx];
  df = df - x_eta_[block_idx]*sbp_[block_idx]->DXi(f[block_idx]) -
       sbp_[block_idx]->DXi(temp);
  return 0.5*invJ_[block_idx]*df;
}

/*
 * To be used in SAT:s where DnT is needed.
 * Dn = Nx*Dx + Ny*Dy at block_idx, boundary side
 * Uses DxT and DyT
 */
//Array MbSbp::DnT(const MbArray& f, int block_idx, Side side) {
//   auto bd_slice = GetBlockBoundarySliceAndSize(block_idx, side);
//   auto normals = GetNormals(block_idx, side);
//
//   std::vector<Array> nx_f (num_blocks());
//   nx_f[block_idx] = Array(shapes(block_idx).Nx,
//                           shapes(block_idx).Ny);
//   nx_f[block_idx][bd_slice.second] = normals.a1.array()*
//                                  f[block_idx][bd_slice.second];
//
//   std::valarray DxT_nx_f = DxT(nx_f, block_idx)[bd_slice.second];
//
//   std::vector<Array> ny_f (num_blocks());
//   ny_f[block_idx] = Array(shapes(block_idx).Nx,
//                           shapes(block_idx).Ny);
//   ny_f[block_idx][bd_slice.second] = normals.a2.array()*
//                                  f[block_idx][bd_slice.second];
//
//   std::valarray DyT_ny_f = DyT(ny_f, block_idx)[bd_slice.second];
//   return Array(1, bd_slice.first, DxT_nx_f + DyT_ny_f);
//}

/*
 * Transpose(Dx) * f on block block_idx. 
 * To be used in SAT:s where DnT is needed.
 */
//Array MbSbp::DxT(const MbArray& f, int block_idx) {
//
//  Array Jinv_f = 0.5*invJ_[block_idx]*f[block_idx];
//  Array y_eta = y_eta_[block_idx];
//  Array y_xi = y_xi_[block_idx];
//
//  return sbp_[block_idx]->DXiT(y_eta*Jinv_f) +
//         y_eta*sbp_[block_idx]->DXiT(Jinv_f) -
//         sbp_[block_idx]->DEtaT(y_xi*Jinv_f) -
//         y_xi*sbp_[block_idx]->DEtaT(Jinv_f);
//}

/*
 * Transpose(Dy) * f on block block_idx.
 * To be used in SAT:s where DnT is needed.
 */
//Array MbSbp::DyT(const MbArray& f, int block_idx) {
//
//  Array Jinv_f = 0.5*invJ_[block_idx]*f[block_idx];
//  Array x_eta = x_eta_[block_idx];
//  Array x_xi = x_xi_[block_idx];
//
//  return sbp_[block_idx]->DEtaT(x_xi*Jinv_f) +
//         x_xi*sbp_[block_idx]->DEtaT(Jinv_f) -
//         sbp_[block_idx]->DXiT(x_eta*Jinv_f) -
//         x_eta*sbp_[block_idx]->DXiT(Jinv_f);
//}


// This function is used in the constructor
void MbSbp::SetNormalsAndBoundaryQuadratures(int block_idx) {
  /*
   * This function sets all normals: nx, ny at each
   * boundary for each block
   * ns_, ne_, nn_, nw_ are vectors.
   * Ex ns_[k] is a std::pair<Array,Array> = nx_s, ny_s 
   * are the normals of the south boundary of block k.
   */

  Array p = sbp_[block_idx]->GetP()*J_[block_idx];
  Array p_inv = 1.0/p;

  double d_eta = sbp_[block_idx]->d_eta();
  double d_xi  = sbp_[block_idx]->d_xi();
  int Nx = shapes(block_idx).Nx;
  int Ny = shapes(block_idx).Ny;

  // west and east
  Array n1(Nx,Ny);
  Array n2(Nx,Ny);
  Array nrm(Nx,Ny);

  for(int i = 0; i < Nx*Ny; ++i) {
    nrm[i]= sqrt(pow(x_eta_[block_idx][i],2) +
            pow(y_eta_[block_idx][i],2));
    n1[i] = y_eta_[block_idx][i]/nrm[i];
    n2[i] = x_eta_[block_idx][i]/nrm[i];
  }

  Array bd_quad = sbp_[block_idx]->GetQuadratureWeights(Ny)*d_eta;
  Side side = Side::w;
  nw_.push_back({ToBlockBoundary(-1.0*n1, block_idx, side),
                 ToBlockBoundary(n2, block_idx, side)});
  bd_quad_west_.push_back(
      ToBlockBoundary(nrm, block_idx, side)*bd_quad);
  p_inv_west_.push_back(ToBlockBoundary(p_inv, block_idx, side));

  side = Side::e;
  ne_.push_back({ToBlockBoundary(n1,block_idx,side),
                 ToBlockBoundary(-1.0*n2, block_idx, side)});

  bd_quad_east_.push_back(
      ToBlockBoundary(nrm, block_idx, side)*bd_quad);
  p_inv_east_.push_back(
      ToBlockBoundary(p_inv, block_idx, side));

  // south and north
  for(int i = 0; i < Nx*Ny; ++i) {
    nrm[i]= sqrt(pow(x_xi_[block_idx][i],2) +
                 pow(y_xi_[block_idx][i],2));
    n1[i] = y_xi_[block_idx][i]/nrm[i];
    n2[i] = x_xi_[block_idx][i]/nrm[i];
  }

  bd_quad = sbp_[block_idx]->GetQuadratureWeights(Nx)*d_xi;

  side = Side::s;
  ns_.push_back({ToBlockBoundary(n1,block_idx, side),
                 ToBlockBoundary(-1.0*n2,block_idx,side)});
  bd_quad_south_.push_back(
      ToBlockBoundary(nrm,block_idx,side)*bd_quad);
  p_inv_south_.push_back(
      ToBlockBoundary(p_inv, block_idx,side));

  side = Side::n;
  nn_.push_back({ToBlockBoundary(-1.0*n1,block_idx,side),
                 ToBlockBoundary(n2,block_idx,side)});
  bd_quad_north_.push_back(
      ToBlockBoundary(nrm,block_idx,side)*bd_quad);
  p_inv_north_.push_back(
      ToBlockBoundary(p_inv, block_idx,side));
}

// -------------------- MbGrid-functions --------------------------
std::pair<int,std::slice> MbSbp::GetBlockBoundarySliceAndSize(
    int block_idx, Side side) {
  return grid_.GetBlockBoundarySliceAndSize(block_idx, side);
}

std::vector< Block > MbSbp::blocks() {
  return grid_.blocks();
}

Block MbSbp::blocks(int block_idx) {
  return grid_.blocks(block_idx);
}

std::vector<Shape> MbSbp::shapes() {
  return grid_.shapes();
}

Shape MbSbp::shapes(int block_idx) {
  return grid_.shapes(block_idx);
}

MbArray MbSbp::Evaluate(std::function<double(double,double)> f) {
  return MbArray(grid_.Evaluate(f));
}

Array MbSbp::ToBlockBoundary(const MbArray& f,
                             int block_idx, Side side) {
  return grid_.ToBlockBoundary(f, block_idx, side);
}

Array MbSbp::ToBlockBoundary(const Array& f,
                             int block_idx, Side side) {
  return grid_.ToBlockBoundary(f, block_idx, side);
}

bool MbSbp::IsFilppedInterface(int interface_idx) {
  return grid_.IsFilppedInterface(interface_idx);
}

//void MbSbp::SetInterpolationOperators() {
//
//  for(auto& interface : grid_.interfaces()) {
//    int N1 = GetPinvAtBoundary(interface.block_idx1,
//                               interface.side1).size();
//    int N2 = GetPinvAtBoundary(interface.block_idx1,
//                               interface.side2).size();
//    switch (order_) {
//      case 4:
//        interp_.push_back(std::make_unique<Interpolation42>(N1,N2));
//         break;
//      default:
//         interp_.push_back(std::make_unique<Interpolation21>(N1,N2));
//     }
//  }
//}
//
/*
 * InterfaceTreatment. If used together with Dx and Dy,
 * this interface treatment gives a skew-symmetric
 * "D-matrix" on a multi-block domain.
 * Use it to encapsulate all interfaces on a multi-block domain.
 * Values from f are added onto df.
 *
 * Use:
 *  void InterfaceTreatment(df, f) -->
 *  values from f at each interface are added onto df.
 */
//void MbSbp::InterfaceTreatment(const MbArray& f,
//                                     MbArray& df,
//                               const Direction& direction)
//{
//
//   int block_idx1, block_idx2, N1, N2;
//   std::pair<Array,Array> normals1, normals2;
//   Side side1, side2;
//   Array Pinv1, Pinv2, Pgamma1, Pgamma2, f1, f2, interface1, interface2, n1, n2;
//
//   int idx = 0; //interface index
//
//   for(auto& interface : grid_.GetInterfaces())
//   {
//      block_idx1 = interface.block_idx1;
//      block_idx2 = interface.block_idx2;
//      side1 = interface.side1;
//      side2 = interface.side2;
//
//      std::pair<int, std::slice> size_and_slice1 = grid_.GetBlockBoundarySliceAndSize(
//                                          block_idx1,side1);
//      std::pair<int, std::slice> size_and_slice2 = grid_.GetBlockBoundarySliceAndSize(
//                                       block_idx2,side2);
//      Pinv1 = GetPinvAtBoundary(block_idx1, side1);
//      Pinv2 = GetPinvAtBoundary(block_idx2, side2);
//      Pgamma1 = GetBoundaryQuadrature(block_idx1, side1);
//      Pgamma2 = GetBoundaryQuadrature(block_idx2, side2);
//      N1 = Pinv1.Size();
//      N2 = Pinv2.Size();
//
//      f1 = ToBlockBoundary(f, block_idx1, side1);
//      f2 = ToBlockBoundary(f, block_idx2, side2);
//
//      // Interfiace X-direction
//      if(direction == Direction::x)
//      {
//
//         n1 = GetNormals(block_idx1, side1).first;
//         n2 = GetNormals(block_idx2, side2).first;
//      }
//      else // Interfiace Y-direction
//      {
//         n1 = GetNormals(block_idx1, side1).second;
//         n2 = GetNormals(block_idx2, side2).second;
//      }
//
//      if(grid_.IsFilppedInterface(idx))
//         f2.Reverse();
//
//      //interface1 = -0.5*Pinv1*Pgamma1*(n1*f1 - 0.5*(n1 - n2)*f2);
//      interface1 = -0.5*Pinv1*Pgamma1*(n1*f1 - 
//                    0.5*(n1*interp_[idx]->Interpolate(f2,N2,N1) - 
//                         interp_[idx]->Interpolate(n2*f2,N2,N1)));
//
//      if(grid_.IsFilppedInterface(idx))
//      {
//         f2.Reverse();
//         f1.Reverse();
//      }
//      
//      //interface2 = -0.5*Pinv2*Pgamma2*(1.0*n2*f2 - 0.5*(n2 - n1)*f1);
//      interface2 = -0.5*Pinv2*Pgamma2*(n2*f2 - 
//                    0.5*(n2*interp_[idx]->Interpolate(f1,N1,N2) - 
//                         interp_[idx]->Interpolate(n1*f1,N1,N2)));
//
//      df[block_idx1][size_and_slice1.second] += interface1.GetArray();
//      df[block_idx2][size_and_slice2.second] += interface2.GetArray();
//      ++idx;
//   }
//}

//MbArray MbSbp::DyAndInterface(const MbArray& f) {
//  MbArray dfdy = Dy(f);
//  InterfaceTreatment(f, dfdy, Direction::y);
//  return dfdy;
//}
