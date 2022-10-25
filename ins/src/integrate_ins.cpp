#include <math.h>
#include <mkl.h>
#include <algorithm>    // std::max
#include "integrate_ins.hpp"
#include "tqdm.hpp"
#include "ins.hpp"

//---------------------------------------------------------------
//------------------- Implicit time integration ------------------
//----------------------------------------------------------------
/*
 * Implicit integration.
 * Integrates state' + spatial_op(t,state) = 
 * from tspan[0] to tspan[1].
 * Input:
 *  o tspan[t0, t1] - limits of the integration
 *  o state0 - initial condition
 *  o dt - step size
 *  o spatial_op - The spatial operator.
 *  o Jv - std::function.
 *    The action of the Jacobian to the spatial operator.
 *
 * Output:
 *  o InsSolution containing the approximation u,v,p at t1.
 */

using valarray = std::valarray<double>;

InsSolution ImplicitTimeIntegraction(const Tspan& tspan,
                               const InsState& init, double dt,
                               Ins& ins, std::string name_base) {

  int NoS  = static_cast<int>(ceil((tspan.t1 - tspan.t0)/dt));
  int steps_per_time = static_cast<int>(ceil(NoS/(tspan.t1 -
                                        tspan.t0)));
  int interval = static_cast<int>(ceil(steps_per_time/15));
  interval = std::max(1, interval);

  auto state {init}, prev_state {init}, prev_prev_state {init};
  double t = tspan.t0;

  int tec_pos = 0;
  double tol = 1e-12;
  for(int i : tq::trange(NoS)) {
    if(i % interval == 0) {
      std::string file_name = name_base + std::to_string(tec_pos);
      ins.ExportToTec(state,file_name);
      ++tec_pos;
    }

    t = tspan.t0 + (i+1)*dt;
    if(i == 0)
      state = InsBdf1Step(state, dt, t, ins, tol);
    else {
       prev_prev_state = prev_state;
       prev_state = state;
       state = InsBdf2Step(prev_prev_state,
                           prev_state, dt, t, ins, tol);
    }
  }
  return {t,state};;
}

//InsSolution SolveSteadyState(const InsState& init,
//                             Ins& ins, std::string name_base) {
//
//  InsState new_state = init, prev_state = init, prev_prev_state = init;
//
//  int tec_pos = 0;
//  double tol = 1e-9;
//  int N = init.u.GetTotalSize();
//  double t = 0;
//
//  auto F = [&ins](double t, const InsState& state)
//  {
//     return ins.SpatialOperator(t,state);
//  };
//
//  auto J_Fv = [&ins](const InsState& state)
//  {
//     return ins.Jacobian(state);
//  };
//
//
//  std::valarray<double> delta;
//  std::valarray<double> b = F(t,new_state);
//  double linf_err = abs(b).max();
//  for(int k = 0; k < 25; k++)
//  {
//     prev_state = new_state;
//
//     std::string file_name = name_base + std::to_string(k);
//     ins.ExportToTec(new_state,file_name);
//
//     //auto J_Fv = [&ins,&b,&F,prev_state](const InsState& state)
//     //{
//     //   double mu = 1e-5;
//     //   DataTypes::InsState tmp = state;
//     //   tmp.u = tmp.u*mu + state.u;
//     //   tmp.v = tmp.v*mu + state.v;
//     //   tmp.p = tmp.p*mu + state.p;
//     //   return (F(0,tmp) - b)/mu;
//     //};
//
//     delta = InsGmresMkl(J_Fv, b, new_state, ins);
//
//     new_state.u -= MbArray(delta[std::slice(0,N,1)], 
//                                 new_state.u.GetShapes());
//
//     new_state.v -= MbArray(delta[std::slice(N,N,1)], 
//                                 new_state.u.GetShapes());
//     new_state.p -= MbArray(delta[std::slice(2*N,N,1)], 
//                                 new_state.u.GetShapes());
//     b = F(t,new_state);
//     linf_err = abs(b).max();
//     std::cout << "err: " << linf_err << std::endl;
//     if(linf_err < tol){ break; }
//  }
//  if (linf_err > tol) {std::cout << "BDF1: GMRES NOT CONVERGED, L_inf err"
//                    << linf_err << std::endl;}
//
//  DataTypes::InsSolution solution;
//  solution.state = new_state;
//  return solution;
//}


InsState InsBdf1Step(const InsState& prev_state,
                     double dt, double t, Ins& ins, double tol) {
  auto new_state = prev_state;
  int N = prev_state.u.GetTotalSize();

  auto F = [&prev_state,&ins,&dt,&N]
            (double t, const InsState& state) {

    valarray T = dt*ins.SpatialOperator(t,state);

    int size = N, start = 0, stride = 1;
    T[std::slice(start,size,stride)] += (state.u.ToValarray() -
                                       prev_state.u.ToValarray());
    start += N;
    T[std::slice(start,size,stride)] +=
      (state.v.ToValarray() - prev_state.v.ToValarray());
    return T;
  };

  auto J_Fv = [&N,&dt,&ins](const InsState& state) {
    valarray T = dt*ins.Jacobian(state);
    int size = N, start = 0, stride = 1;
    T[std::slice(start,size,stride)] += state.u.ToValarray();
    start += N;
    T[std::slice(start,size,stride)] += state.v.ToValarray();
    return T;
  };
  valarray delta;
  double linf_err;

  auto b {F(t,new_state)};
  linf_err = abs(b).max();

  for(int k = 0; k < 25; ++k) {
    delta = InsGmresMkl(J_Fv, b, new_state, ins);
    auto shapes {new_state.u.shapes()};

    new_state.u -= MbArray(shapes, delta[std::slice(0,N,1)]);
    new_state.v -= MbArray(shapes, delta[std::slice(N,N,1)]);
    new_state.p -= MbArray(shapes, delta[std::slice(2*N,N,1)]);

    b = F(t,new_state);
    linf_err = abs(b).max();
    std::cout << linf_err << "\n";
    if(linf_err < tol){ break; }
  }
  if (linf_err > tol) {
    std::cout << "BDF1: GMRES NOT CONVERGED, L_inf err"
             << linf_err << std::endl;
  }
  return new_state;
}

InsState InsBdf2Step(const InsState& prev_prev_state,
                     const InsState& prev_state,
                     double dt, double t, Ins& ins, double tol) {
  InsState new_state = prev_state;
  int N = prev_state.u.GetTotalSize();

  auto F = [&](double t, const InsState& state) {
    valarray T = (2.0/3.0)*dt*ins.SpatialOperator(t,state);

    int size = N, start = 0, stride = 1;
    T[std::slice(start,size,stride)] += state.u.ToValarray() -
                         (4.0/3.0)*prev_state.u.ToValarray() +
                         (1.0/3.0)*prev_prev_state.u.ToValarray();
    start += N;
    T[std::slice(start,size,stride)] += state.v.ToValarray() -
                         (4.0/3.0)*prev_state.v.ToValarray() +
                         (1.0/3.0)*prev_prev_state.v.ToValarray();
    return T;
  };

  auto J_Fv = [&](const InsState& state) {
     valarray T = (2.0/3.0)*dt*ins.Jacobian(state);
     int size = N, start = 0, stride = 1;
     T[std::slice(start,size,stride)] += state.u.ToValarray();
     start += N;
     T[std::slice(start,size,stride)] += state.v.ToValarray();
     return T;
  };

  valarray delta;
  double linf_err;

  valarray b = F(t,new_state);
  linf_err = abs(b).max();
  for(int k = 0; k < 20; ++k) {
    auto shapes = new_state.u.shapes();
    delta = InsGmresMkl(J_Fv, b, new_state, ins);

    new_state.u -= MbArray(shapes,delta[std::slice(0,N,1)]);

    new_state.v -= MbArray(shapes,delta[std::slice(N,N,1)]);
    new_state.p -= MbArray(shapes, delta[std::slice(2*N,N,1)]);

    b = F(t,new_state);
    linf_err = abs(b).max();
    if(linf_err < tol){ break;}
  }
  if (linf_err > tol) {
    std::cout << "BDF2: GMRES NOT CONVERGED, L_inf err = "
              << linf_err << std::endl;
  }
  return new_state;
}


/*
 * Wrapper of Intel MKL dfgmres. Solves the system Ax = b. 
 * Input: 
 *   Ax - A function returning the action of A on a MbArray
 *   b - right-hand side
 *   x0 - initial guess
 */
valarray InsGmresMkl(
      std::function<valarray(const InsState&)> Jv,
      const valarray& b, const InsState& cur_state,
      Ins& ins) {

  int N = b.size();
  int iter;
  int max_iter = 1500, subspace_size = 150;
  double rel_tol = 1e-6, abs_tol = 1e-6;

  valarray sol(N), b_arr = b;
  int RCI_request;
  std::valarray<int> ipar(128);
  std::vector<double> dpar(128);

  int tmp_size = (2*subspace_size + 1)*N +
                 subspace_size*(subspace_size + 9)/2 + 1;
  valarray tmp(tmp_size);

  dfgmres_init(&N, &sol[0], &b_arr[0], &RCI_request,
               &ipar[0], &dpar[0], &tmp[0]);

  ipar[0]  = N;
  ipar[1]  = 6; // error messages
  ipar[2]  = 1;
  ipar[3]  = 0; // current iteration number
  ipar[4]  = max_iter; // maximum iterations
  ipar[5]  = 1; // error messages
  ipar[6]  = 1; // error messages
  ipar[7]  = 1; // 1--> check ipar[3] < ipar[4], defalut: 1
  ipar[8]  = 0;
  ipar[9]  = 0;
  ipar[10] = 0; // 0 --> not use preconditioner
  ipar[11] = 0;
  ipar[12] = 0;
  ipar[13] = 0;
  ipar[14] = subspace_size;

  dpar[0] = rel_tol; //relative tolerance
  dpar[1] = abs_tol; //absolute tolerance

  dfgmres_check(&N, &sol[0], &b_arr[0], &RCI_request,
                &ipar[0], &dpar[0], &tmp[0]);

GMRES: {
  dfgmres(&N, &sol[0], &b_arr[0], &RCI_request,
          &ipar[0], &dpar[0], &tmp[0]);

    switch (RCI_request) {
       case 0: {
         goto COMPLETE;
       }
       case 1: {
         auto shapes = cur_state.u.shapes();
         int start = ipar[21]-1;
         InsState state =
           ins.ValArrayToInsState(tmp[std::slice(start,N,1)]);
         start = ipar[22]-1;
         tmp[std::slice(start,N,1)] = Jv(state);
         goto GMRES;
       }
       case 4: {
         if(dpar[6] < 1e-12)
           goto COMPLETE;
         else
           goto GMRES;
       }
    }
  }

COMPLETE: {
  ipar[12] = 0;
  dfgmres_get(&N, &sol[0], &b_arr[0], &RCI_request,
              &ipar[0], &dpar[0], &tmp[0], &iter);
  return sol;
  }
}
