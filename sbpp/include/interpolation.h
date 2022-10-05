/*
 * interpolation.h
 * Interpolation operators from 
 * "Mattsson and Karpenter, Stable and Accurate Interpolation Operators 
 *                          for High-Order Multi-Block Finite-Difference Methods"
 *  Created on: 16 Jan. 2021
 *      Author: Fredrik Laurén
 */

#pragma once
#include <iostream>
#include "array2D.h"

class Interpolation
{
  public:
     Interpolation(){};
     Interpolation(int Nc, int Nf);
     virtual ~Interpolation(){};

     // interpolate array from grid of size N1 to grid of size N2, using either
     // Coarse2Fine or Fine2Coarse
     Array2D Interpolate(const Array2D& array, int N1, int N2);

  protected:
     virtual Array2D Coarse2Fine(const Array2D& coarse) = 0;
     virtual Array2D Fine2Coarse(const Array2D& fine) = 0;
     int Nc_, Nf_;
};
