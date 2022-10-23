#pragma once

#include <iostream>
#include <vector>
#include "array.hpp"
#include "mbarray.hpp"
#include "ins.hpp"

struct InsSolution {
  double time;
  InsState state;
};

struct Tspan {
  double t0, t1;
};

//--------------------------------------------------------------
//----------------- Implicit time integration ------------------
//---------------------------------------------------------------
InsSolution ImplicitTimeIntegraction(const Tspan& tspan,
                    const InsState& init, double dt,
                    Ins& ins, std::string name_base);

InsSolution SolveSteadyState(const InsState& init,
                             Ins& ins, std::string name_base);

InsState InsBdf1Step(const InsState& prev_state,
                     double dt, double t, Ins& eul, double tol);

InsState InsBdf2Step(const InsState& prev_prev_state,
                     const InsState& prev_state,
                     double dt, double t, Ins& eul, double tol);

std::valarray<double> InsGmresMkl(
      std::function<std::valarray<double>(const InsState&)> Jv,
      const std::valarray<double>& b, const InsState& cur_state,
      Ins& eul);
