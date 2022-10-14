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

/*
 * Explicit integration.
 * Optional to save every XX solution step
 * so that 15 snapshots are saved per time unit.
 *
 * Integrates y' + odefun(t,y) = 0 from tspan[0] to tspan[1].
 * Input:
 *   o tspan[t0, t1] - limits of the integration
 *   o y0 - initial condition
 *   o dt - step size
 *   o odefun - std::function so that
 *              MbArray y' = - odefun(double t, MbArray y)
 *   o write_to_file - true/false
 *
 *   o export_to_tec - std::function to write y to .tec format.
 *                  void export_to_tec(MbArray y, std::string name)
 *
 *   o std::string name_base - name_base of the files to be saved.
 *                     Each saved file will get an index extension.
 *                     Ex: name_base = "save/sol"
 *                     Then each file will be saved as 
 *                     "save/sol1, "save/sol2", ...
 *
 * Output:
 *   o Solution containing the approximation to y(t1) and t1.
 */

Solution ExplicitIntegration(std::vector<double>& tspan,
                             const MbArray& y0, double dt,
      std::function<MbArray(double t, const MbArray&)> odefun,
      bool write_to_file = false,
      std::function<void(const MbArray&,const std::string&)>
      export_to_tec = NULL,
      std::string name_base = "none");

/*
 * 4th order RK-step. Used in ExplicitIntegration.
 */
MbArray RK4Step(double t, const MbArray& y, double dt,
           std::function<MbArray(double, const MbArray&)> odefun);

//// --------------------------------------------------------------
//// ------------------ Implicit time integration ------------------
//// ---------------------------------------------------------------
Solution ImplicitTimeIntegraction(std::vector< double >& tspan,
                                     const MbArray& y0, double dt,
        std::function<MbArray(double, const MbArray&)> odefun,
        std::function<MbArray(double, const MbArray&)> Jv);

//Solution ImplicitTimeIntegractionSaveSolution(std::vector< double >& tspan,
//                                     const MbArray& y0, double dt,
//        std::function<MbArray(double, const MbArray&)> odefun,
//        std::function<MbArray(double, const MbArray&)> Jv,
//        std::function<void(const MbArray&, const std::string&)> export_to_tec,
//        std::string name_base);
//
//
MbArray BDF1Step(const MbArray& y_prev,
                      double dt, double t,
        std::function<MbArray(double, const MbArray&)> odefun,
        std::function<MbArray(double, const MbArray&)> Jv);

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

// documentation:
// https://software.intel.com/content/www/us/en/develop/documentation/onemkl-developer-reference-c/top.html
//
//The algorithm for MKL's GMRES is inspired by 
// http://sep.stanford.edu/sep/claudio/Research/Prst_ExpRefl/ShtPSPI/intel/mkl/10.0.3.020/examples/solver/source/dcsrilu0_exampl1.c
MbArray GMRESMKL(std::function<MbArray(const MbArray&)> Ax,
                 const MbArray& b, const MbArray& x0);
