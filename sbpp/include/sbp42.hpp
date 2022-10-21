//============================================================
// Name        : sbp42.hpp
// Date        : June 2022
// Author      : Fredrik Laurén
// Description : hpp file for the class Sbp21.
//============================================================

/*
 * Class for Sbp42 operators.
 * Inherit from base class Sbp.
 */

#pragma once
#include <iostream>
#include "sbp.hpp"

class Sbp42:public Sbp {
  public:
    Sbp42(){};
    Sbp42(int Nx, int Ny, double d_xi, double d_eta);
    ~Sbp42(){};

    double Integrate(const Array&) const;
    Array GetQuadratureWeights(int size) const;

  private:
    void Diff(const double* f, double* df, int, int, double) const;
    void DT(const double* f, double* dtf, int, int, double);

    // weights
    double h1_ = 17.0/48.0;
    double h2_ = 59.0/48.0;
    double h3_ = 43.0/48.0;
    double h4_ = 49.0/48.0;
};
