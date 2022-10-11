/*
 * Advection.cpp
 *
 *  Created on: 15 jan. 2020
 *      Author: frela05
 */

#include "advection.hpp"

Advection::Advection(std::vector<Block> blocks, int order,
                     std::pair<double,double> a) {
  sbp_ = std::make_unique<MbSbp>(blocks, order);
  a_   = a;
}

double Advection::Analytic(double t, double x, double y) {
  return sin(2*M_PI*x-t)*sin(2*M_PI*y-t);
}

MbArray Advection::AnalyticVec(double t) {
  auto func = [&t](double x, double y) { 
    return Analytic(t,x,y);
  };
  return sbp_->Evaluate(func);
}

double Advection::Force(double t, double x, double y,
                        double a1, double a2) {
  double ut = - cos(2*M_PI*x-t)*sin(2*M_PI*y-t) -
                sin(2*M_PI*x-t)*cos(2*M_PI*y-t);
  double ux = 2*M_PI*cos(2*M_PI*x-t)*sin(2*M_PI*y-t);
  double uy = 2*M_PI*sin(2*M_PI*x-t)*cos(2*M_PI*y-t);
  return ut + a1*ux + a2*uy;
}

MbArray Advection::Force(double t) {
  auto func = [a = a_, tau = t](double x, double y) {
     return Force(tau,x,y, a.first, a.second);
   };
  return sbp_->Evaluate(func);
}

/*
 * Build the right-hand side.
 * ut = rhs
 * where
 * rhs = -(a1 ux + a2 uy) + SAT + Force
 */
MbArray Advection::Rhs(double t, const MbArray& sol, bool mms) {
  MbArray d_sol_dx = sbp_->DxAndInterface(sol);
  MbArray d_sol_dy = sbp_->DyAndInterface(sol);

  MbArray rhs = -1.0*(a_.first*d_sol_dx + a_.second*d_sol_dy);

  auto func_analytic = [tau = t](double x, double y) {
    return Analytic(tau,x,y);
  };
  if(mms) {
    auto g_func = [&](double t) {
      return sbp_->Evaluate(func_analytic);
    };
    ApplySat(t,sol,rhs,g_func);
  }
  else {
    auto g_func_helper = [](double x, double y) {
      return  1;
    };
    auto g_func = [&](double t){
      return sbp_->Evaluate(g_func_helper);
    };
    ApplySat(t,sol,rhs,g_func);
  };

  if (mms)
    return rhs + Force(t);
  return rhs;
}

/*
 * Add SAT to rhs
 * The boundary condition is u = g,
 * where g is the analytic solution
 * The boundary condition is applied to inflow boundaries only:
 * SAT = (a*n - |a*n|)*Pinv*Pgamma*(u - g)
 */
void Advection::ApplySat(double t, const MbArray& sol,
                         MbArray& rhs,
      std::function<MbArray(double t)> bd_data) {

  Array bd_quad, Pinv, Pgamma;
  ArrayPair normals;

  double direction;
//   MbArray analytic = AnalyticVec(t);
  MbArray analytic = bd_data(t);

  for(auto& bd : sbp_->boundaries()) {
    bd_quad = sbp_->GetBoundaryQuadrature(bd.block_idx, bd.side);
    normals = sbp_->GetNormals(bd.block_idx, bd.side);

    auto bd_slice = sbp_->GetBlockBoundarySliceAndSize(
                                bd.block_idx, bd.side);

    Pinv = sbp_->GetPinvAtBoundary(bd.block_idx, bd.side);
    Pgamma = sbp_->GetBoundaryQuadrature(bd.block_idx, bd.side);

    Shape shape = sbp_->shapes(bd.block_idx);

    Array direction = a_.first*normals.a1 + a_.second*normals.a2;
    direction = 0.5*(direction - Abs(direction));
    Array sol_bd = Array(1,bd_slice.first,
                         sol[bd.block_idx][bd_slice.second]);
    Array analytic_bd = Array(1,bd_slice.first,
                       analytic[bd.block_idx][bd_slice.second]);

    Array pen = direction*Pinv*Pgamma*(sol_bd - analytic_bd);

    rhs[bd.block_idx][bd_slice.second] +=  pen.array();
  }
}

/*
 * Jacobian for rhs multiplied by a vector f
 *  J f  = (a1 Dxf + a2 Dy f) + SAT(f) 
 *  The SAT-term is evaluated for g = 0 (no data), 
 *  since it is not included in the Jacobian
 */
MbArray Advection::Jacobian(double t, const MbArray& f) {
  MbArray d_f_dx = sbp_->DxAndInterface(f);
  MbArray d_f_dy = sbp_->DyAndInterface(f);
  MbArray rhs = -1.0*(a_.first*d_f_dx + a_.second*d_f_dy);

  auto data_func = [tau = t](double x, double y){
    return 0;
  };
  auto g_func = [&](double t){
    return sbp_->Evaluate(data_func);
  };
  ApplySat(t,f,rhs, g_func);
  return rhs;
}

double Advection::ComputeError(const double t,
                               const MbArray& computed) {

  MbArray err, analytic = AnalyticVec(t);
  err = computed - analytic;
  err = err*err;
 return sqrt(sbp_->Integrate(err));
}

//void Advection::ExportToTec(const MbArray& arr,
//                            const std::string name) {
//  ExportToTec(sbp_->grid(), arr, name);
//}

void Advection::ExportToTec(const std::vector<MbArray>& gf_vec,
                            const std::string name_base) {
  std::string name;
  for (int i = 0; i < gf_vec.size(); ++i) {
    name = name_base + std::to_string(i);
    ::ExportToTec(sbp_->grid(), gf_vec[i], name);
  }
}
