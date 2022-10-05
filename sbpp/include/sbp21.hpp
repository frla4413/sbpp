//============================================================
// Name        : sbp21.hpp
// Date        : June 2022
// Author      : Fredrik Laurén
// Description : hpp file for the class Sbp21.
//============================================================

/*
 * Class for Sbp21 operators.
 * Inherit from base class Sbp.
 * See sbp.hpp for further details.
 */

#pragma once
#include <iostream>
#include "sbp.hpp"

class Sbp21:public Sbp {
  public:
    Sbp21(){};
    Sbp21(int Nx, int Ny, double d_xi, double d_eta);
    ~Sbp21(){};

    double Integrate(const Array&) const;
    Array GetQuadratureWeights(int N) const;

  private:
    void Diff(const double* f,double* df, int, int, double) const;
//    void DT(const double* f,double* dtf, int, int, double);
    // weights
    //double h1_ = 0.5;
};
