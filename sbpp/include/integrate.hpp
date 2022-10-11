#pragma once

#include <iostream> 
#include <vector>
#include "array.hpp"
#include "mbarray.hpp"

struct Solution {
  double t;
  MbArray y;
};

// --------------------------------------------------------------
// ----------------- Explicit time integration ------------------
// --------------------------------------------------------------
Solution ExplicitIntegration(std::vector<double>& tspan,
                             const MbArray& y0, double dt,
      std::function<MbArray(double t, const MbArray&)> odefun);

//Solution ExplicitIntegrationSaveSolution(std::vector< double >& tspan, 
//                                                    const MbArray& y0, 
//                                                    double dt,
//                 std::function<MbArray(double t, const MbArray&)> odefun,
//        std::function<void(const MbArray&, const std::string&)> export_to_tec,
//        std::string name_base);

MbArray RK4Step(double t, const MbArray& y, double dt,
           std::function<MbArray(double, const MbArray&)> odefun);

//// -------------------------------------------------------------------------
//// ---------------------------- Implicit time integration ------------------
//// -------------------------------------------------------------------------
//Solution ImplicitTimeIntegraction(std::vector< double >& tspan,
//                                     const MbArray& y0, double dt,
//        std::function<MbArray(double, const MbArray&)> odefun,
//        std::function<MbArray(double, const MbArray&)> Jv);
//
//Solution ImplicitTimeIntegractionSaveSolution(std::vector< double >& tspan,
//                                     const MbArray& y0, double dt,
//        std::function<MbArray(double, const MbArray&)> odefun,
//        std::function<MbArray(double, const MbArray&)> Jv,
//        std::function<void(const MbArray&, const std::string&)> export_to_tec,
//        std::string name_base);
//
//
//MbArray BDF1Step(const MbArray& y_prev,
//                      double dt, double t,
//        std::function<MbArray(double, const MbArray&)> odefun,
//        std::function<MbArray(double, const MbArray&)> Jv);
//
//MbArray BDF2Step(const MbArray& y_prev_prev,
//                      const MbArray& y_prev,
//                      double dt, double t,
//        std::function<MbArray(double, const MbArray&)> odefun,
//        std::function<MbArray(double, const MbArray&)> Jv);
//
//MbArray BDF3Step(const MbArray& y_prev_prev_prev,
//                      const MbArray& y_prev_prev,
//                      const MbArray& y_prev, 
//                      double dt, double t,
//        std::function<MbArray(double, const MbArray&)> odefun,
//        std::function<MbArray(double, const MbArray&)> Jv);
//
//// documentation: 
//// https://software.intel.com/content/www/us/en/develop/documentation/onemkl-developer-reference-c/top.html
////
////The algorithm for MKL's GMRES is inspired by 
//// http://sep.stanford.edu/sep/claudio/Research/Prst_ExpRefl/ShtPSPI/intel/mkl/10.0.3.020/examples/solver/source/dcsrilu0_exampl1.c
//MbArray GMRESMKL(std::function<MbArray
//                   (const MbArray&)> Ax,
//                    const MbArray& b, const MbArray& x0);
