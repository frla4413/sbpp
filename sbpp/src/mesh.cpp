//============================================================
// Name        : mesh.cpp
// Date        : May 2022
// Author      : Fredrik Laurén
// Description : source file for meshes.
//============================================================

#include "mesh.hpp"

Block CartesianGrid(int Nx, int Ny) {

  if (Nx < 1 || Ny < 1)
    throw std::invalid_argument("Nx, Ny > 1 required");

  double d_xi = 1.0/static_cast<double>(Nx-1);
  double d_eta = 1.0/static_cast<double>(Ny-1);

  Array x{Nx, Ny}, y{Nx, Ny};

  int pos = 0;
  for(int i = 0; i < Nx; ++i) {
    for(int j = 0; j < Ny; ++j) {
      x[pos] = i*d_xi;
      y[pos] = j*d_eta;
      ++pos;
    }
  }
  return {x, y};
}

Block CircleSector(int Nr, int Nth, double r_in, double r_out,
                   double angle_0, double angle_1) {
  double d_r  = (r_out-r_in)/(double)(Nr-1);
  double d_th = (angle_1-angle_0)/(double)(Nth-1);
  Array x{Nr,Nth};
  Array y{Nr,Nth};

  int pos    = 0;
  double th  = angle_0;
  double r   = r_in;

  for(int i = 0; i < Nr; i++) {
    th = angle_0;
      for(int j = 0; j < Nth; j++) {
        x[pos] = r*cos(th);
        y[pos] = r*sin(th);
        ++pos;
        th += d_th;
      }
      r += d_r;
  }
  return {x,y};
}

std::vector<Block> Annulus(int N, double r_in, double r_out) {
  std::vector<Block> blocks;
  blocks.push_back(CircleSector(N,N,r_in,r_out,0.0,0.5*M_PI));
  blocks.push_back(CircleSector(N,N,r_in,r_out,0.5*M_PI, M_PI));
  blocks.push_back(CircleSector(N,N,r_in,r_out,M_PI,1.5*M_PI));
  blocks.push_back(CircleSector(N,N,r_in,r_out,1.5*M_PI,2*M_PI));
  return blocks;
}
