//============================================================
// Name        : compare_mkl_valarray.cpp
// Date        : May 2022
// Author      : Fredrik Laurén
// Description : Performance test of MKL:s vector operations and
//               std::valarray.
//============================================================
#include <iostream>
#include <valarray>
#include <vector>
#include <chrono>
#include "mkl.h"

/*
 * Minimalistic script to comapre
 *  i) MKL VM Mathematical Vector Functions
 *  ii) std::valarray
 *
 *  Results shows that std::valarray is at least twice as fast
 *  compared to MKL for large test cases.
 */

int main() {

  int Nx = 100000, Ny = 1, N = 100000;
  std::vector<double> a(Nx,3.14);
  std::vector<double> b(Nx,2.46);

  std::vector<double> c(Nx);

  // time MKL.
  auto end = std::chrono::steady_clock::now();
  auto start = std::chrono::steady_clock::now();

  for(int i = 0; i < N; i++)
    vdMul(a.size(), &a[0], &b[0], &c[0]);

  end = std::chrono::steady_clock::now();
  std::chrono::duration<double> elapsed_seconds = end-start;
  std::cout << "time MKL: " << elapsed_seconds.count() << "s\n";

  std::valarray<double> a2(3.14,Nx);
  std::valarray<double> b2(2.456,Nx);
  std::valarray<double> c2(1,Nx);

  // time std::valarray
  start = std::chrono::steady_clock::now();

  for(int i = 0; i < N; i++)
    c2 = a2 * c2;
  end = std::chrono::steady_clock::now();
  elapsed_seconds = end-start;
  std::cout << "time array: " << elapsed_seconds.count() << "s\n";
  return 0;

}
