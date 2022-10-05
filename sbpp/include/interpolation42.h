/*
 * interpolation42.h
 *
 *  Created on: 18 Jan. 2021
 *      Author: Fredrik Laurén
 */

#pragma once
#include "array2D.h"
#include <valarray>
#include "interpolation.h"


class Interpolation42:public Interpolation
{
  public:
     Interpolation42(){}; Interpolation42(int Nc, int Nf);
     ~Interpolation42(){};

  private:
     Array2D Coarse2Fine(const Array2D& coarse);
     Array2D Fine2Coarse(const Array2D& fine);
     void SetInterpolationWeights();

     // weights for Coarse2Fine()
     std::valarray<double> C2F0_, C2F1_,C2F2_,C2F3_,C2F4_,
                           C2F5_,C2F6_,C2F7_,C2F8_,C2F9_,
                           C2F_stencil_;

     std::valarray<double> C2F0_rev_, C2F1_rev_,C2F2_rev_,C2F3_rev_,C2F4_rev_,
                           C2F5_rev_,C2F6_rev_,C2F7_rev_,C2F8_rev_,C2F9_rev_;

     // weights for Fine2Coarse()
     std::valarray<double> F2C0_, F2C1_,F2C2_,F2C3_,F2C_stencil_;
     std::valarray<double> F2C0_rev_, F2C1_rev_,F2C2_rev_,F2C3_rev_;
};
