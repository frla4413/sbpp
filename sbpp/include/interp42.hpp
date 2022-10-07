//============================================================
// Name        : interp42.hpp
// Date        : October 2022
// Author      : Fredrik Laurén
// Description : hpp file for the class interp.
//============================================================
/*
 * Description
 * (4,2) order interpolation operator. Extends Interp.
 *
 * Member functions:
 *  o fine = Coarse2Fine(coarse)
 *  o coarse = Fine2Coarse(fine)
 *
 * Member variables:
 *  o Nc_ - N course
 *  o Nf_ - N fine
*/

#pragma once
#include "array.hpp"
#include "interp.hpp"

class Interp42:public Interp {
  public:
    Interp42(){}; Interp42(int Nc, int Nf);
    ~Interp42(){};

  private:
    Array Coarse2Fine(const Array& coarse);
    Array Fine2Coarse(const Array& fine);
    void SetInterpWeights();

    // weights for Coarse2Fine()
    std::valarray<double> C2F0_, C2F1_,C2F2_,C2F3_,C2F4_,
                          C2F5_,C2F6_,C2F7_,C2F8_,C2F9_,
                          C2F_stencil_;

    std::valarray<double> C2F0_rev_, C2F1_rev_,C2F2_rev_,
                          C2F3_rev_,C2F4_rev_,
                          C2F5_rev_,C2F6_rev_,C2F7_rev_,
                          C2F8_rev_,C2F9_rev_;

    // weights for Fine2Coarse()
    std::valarray<double> F2C0_, F2C1_,F2C2_,F2C3_,F2C_stencil_;
    std::valarray<double> F2C0_rev_, F2C1_rev_,F2C2_rev_,F2C3_rev_;
};
