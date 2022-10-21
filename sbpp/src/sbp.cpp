/*
 * sbp.cpp
 *
 *  Created on: 13 dec. 2019
 *      Author: Fredrik Laurén
 */

#include "sbp.hpp"
Sbp::Sbp(int Nx, int Ny, double d_xi, double d_eta) :
  Nx_(Nx), Ny_(Ny), d_xi_(d_xi), d_eta_(d_eta)
{}

double Sbp::d_xi() const {
  return d_xi_;
}

double Sbp::d_eta() const {
  return d_eta_;
}

int Sbp::Nx() const {
  return Nx_;
}

int Sbp::Ny() const {
  return Ny_;
}

Array Sbp::DEta(const Array& f) const {
  int stride = 1;
  int length = Ny_;

  Array dfEta{Nx_,Ny_};

  for(int i = 0; i < Nx_; ++i) {
    Diff(&f[i*Ny_], &dfEta[i*Ny_], stride, length, d_eta_);
  }
  return dfEta;
}

Array Sbp::DEtaT(const Array& f) {
// DyT on a rectangular grid
   int stride = 1;
   int length = Ny_;

   auto dfEtaT {Array(Nx_,Ny_)};

   const double* f_pointer;

   for(int i = 0; i < Nx_; ++i) {
     DT(&f[i*Ny_], &dfEtaT[i*Ny_], stride, length, d_eta_);
   }
   return dfEtaT;
}

Array Sbp::DXi(const Array& f) const {
  int stride = Ny_;
  int length = Nx_;

  Array dfXi{Nx_,Ny_};

  for(int i = 0; i < Ny_; ++i) {
    Diff(&f[i], &dfXi[i], stride, length, d_xi_);
  }
  return dfXi;
}

Array Sbp::DXiT(const Array& f) {
// DxT on a rectangular grid
  int stride = Ny_;
  int length = Nx_;

  auto dfXiT {Array(Nx_,Ny_)};

  for(int i = 0; i < Ny_; ++i) {
     DT(&f[i], &dfXiT[i], stride, length, d_xi_);
  }
  return dfXiT;
}

Array Sbp::GetP() const {

  auto weights_x {GetQuadratureWeights(Nx_)};
  auto weights_y {GetQuadratureWeights(Ny_)};

  std::valarray<double> p(Nx_*Ny_);

  // iterate along columns

  int pos = 0;
  for(int i = 0; i < Nx_; ++i) {
    for(int j = 0; j < Ny_; ++j) {
      p[pos] = weights_x[i]*weights_y[j];
      ++pos;
    }
  }
  return {Nx_,Ny_, p*(d_xi_*d_eta_)};
}
