/*
 * Advection.h
 *
 *  Created on: 15 jan. 2020
 *      Author: frela05
 *
 *  Class for the problem u_t + a1*u_x + a2*u_y = 0
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


class Advection {
  public:
   Advection(MbSbp* sbp, std::pair<double,double> a);

   Advection(std::vector<Block> blocks, int order,
             std::pair<double,double> a);
   MbArray AnalyticVec(double t);

   MbArray Force(double t);

   MbArray Rhs(double t, const MbArray& f, bool mms);

   void ApplySat(double t, const MbArray& sol,
                 MbArray& rhs, std::function<MbArray(double t)>);

   void WriteToFile(const MbArray& sol,
                    const std::string name_base);

   void SolutionToFile(const std::vector<Array>& sol,
                       const std::string name_base,
                       const int stride);

   double ComputeError(const double t, const MbArray& computed);

//   void ExportToTec(const MbArray& gf,const std::string name);

   void ExportToTec(const std::vector<MbArray>& gf_vec,
                    const std::string name_base);

   MbArray Jacobian(double t, const MbArray& f);

  private:

   double static Analytic(double t, double x, double y);
   double static Force(double t, double x, double y,
                       double a1, double a2);

   std::pair<double,double> a_;
   std::unique_ptr<MbSbp> sbp_;
  // MbSbp* sbp_;
};
