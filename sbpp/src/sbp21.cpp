/*
 * sbp21.h
 *
 *  Created on: 13 dec. 2019
 *      Author: Fredrik Laurén
 */

#include "sbp21.hpp"

Sbp21::Sbp21(int Nx, int Ny, double dXi, double dEta) :
  Sbp(Nx,Ny,dXi,dEta){}

void Sbp21::Diff(const double* f, double* df,
                 int stride, int length, double h) const {
  df[0] = (f[stride]-f[0])/h;
  df[stride*(length-1)] = (f[stride*(length-1)]-
                           f[stride*(length-2)])/h;

  for(int i = 1; i < length-1; i++)
    df[stride*i] = (f[stride*(i+1)] - f[stride*(i-1)])/(2*h);
}

void Sbp21::DT(const double* f, double* dtf,
               int stride, int length, double h) {
  dtf[0] = (-f[stride]/2-f[0])/h;
  dtf[1] = (-f[stride]/2+f[0])/h;
  dtf[stride*(length-2)] = (-f[stride*(length-1)] +
                            f[stride*(length-2)]/2)/h;
  dtf[stride*(length-1)] = (f[stride*(length-1)] +
                            f[stride*(length-2)]/2)/h;
  for(int i = 2; i < length-2; ++i) {
    dtf[stride*i] = (-f[stride*(i+1)] + f[stride*(i-1)])/(2*h);
  }
}

Array Sbp21::GetQuadratureWeights(int N) const {
  std::valarray <double> weights(1,N);
  weights[0] = 0.5;
  weights[N-1] = 0.5;
  return {1, N, weights};
}

double Sbp21::Integrate(const Array& f) const {

  double sum_tot = 0;

  // integrate along columns
  for(int i = 0; i < Nx_; ++i) {
    double sum  = 0;
    int start = i*Ny_; //start at column i

    sum += 0.5*(f[start] + f[start+Ny_-1]);

    for(int j = 1; j < Ny_-1; ++j)
      sum += f[start+j];

    sum *= d_eta_;

    if(i == 0 || i == Nx_-1)
      sum_tot += 0.5*sum*d_xi_;
    else
      sum_tot += sum*d_xi_;
  }
  return sum_tot;
}
