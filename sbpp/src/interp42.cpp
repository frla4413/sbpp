//============================================================
// Name        : interp42.cpp
// Date        : October 2022
// Author      : Fredrik Laurén
// Description : source file for the class interp42.
//============================================================


#include "interp42.hpp"
#include <stdexcept>

Interp42::Interp42(int Nc, int Nf) : Interp(Nc, Nf) {
  SetInterpWeights();
}

Array Interp42::Fine2Coarse(const Array& fine) {
  Array coarse(1, (fine.size() + 1)/2);
  int fine_ind = 6;
  int N = fine.size();

  for(int i = 0; i < coarse.size(); i++) {
    if( i == 0) {
      coarse[i] = (F2C0_*fine[std::slice(0,10,1)]).sum();
    }
    else if (i == 1) {
      coarse[i] = (F2C1_*fine[std::slice(0,10,1)]).sum();
    }
    else if (i == 2) {
      coarse[i] = (F2C2_*fine[std::slice(0,10,1)]).sum();
    }
    else if (i == coarse.size() - 3) {
      coarse[i] = (F2C2_rev_*fine[std::slice(N-10,10,1)]).sum();
    }
    else if (i == coarse.size() - 2) {
      coarse[i] = (F2C1_rev_*fine[std::slice(N-10,10,1)]).sum();
    }
    else if (i == coarse.size() - 1) {
      coarse[i] = (F2C0_rev_*fine[std::slice(N-10,10,1)]).sum();
    }
    else {
      coarse[i] = (F2C_stencil_*
                   fine[std::slice(fine_ind-3,7,1)]).sum();
       fine_ind += 2;
    }
  }
  return coarse;
}

Array Interp42::Coarse2Fine(const Array& coarse) {

  Array fine(1, (coarse.size()-1)*2 + 1);
  int coarse_ind = 5;
  int N = coarse.size();
  for(int i = 0; i < fine.size(); i++) {
     if( i == 0) {
       fine[i] = (C2F0_*coarse[std::slice(0,3,1)]).sum();
     }
     else if (i == 1) {
       fine[i] = (C2F1_*coarse[std::slice(0,3,1)]).sum();
     }
     else if (i == 2) {
       fine[i] = (C2F2_*coarse[std::slice(0,3,1)]).sum();
     }
     else if (i == 3) {
       fine[i] = (C2F3_*coarse[std::slice(0,4,1)]).sum();
     }
     else if (i == 4) {
       fine[i] = (C2F4_*coarse[std::slice(0,3,1)]).sum();
     }
     else if (i == 5) {
       fine[i] = (C2F5_*coarse[std::slice(0,5,1)]).sum();
     }
     else if (i == 6) {
       fine[i] = (C2F6_*coarse[std::slice(0,4,1)]).sum();
     }
     else if (i == 7) {
       fine[i] = (C2F7_*coarse[std::slice(0,6,1)]).sum();
     }
     else if (i == 8) {
       fine[i] = (C2F8_*coarse[std::slice(0,5,1)]).sum();
     }
     else if (i == 9) {
       fine[i] = (C2F9_*coarse[std::slice(0,7,1)]).sum();
     }
     else if (i == fine.size() - 1) {
       fine[i] = (C2F0_rev_*coarse[std::slice(N-3,3,1)]).sum();
     }
     else if (i == fine.size() - 2) {
       fine[i] = (C2F1_rev_*coarse[std::slice(N-3,3,1)]).sum();
     }
     else if (i == fine.size() - 3) {
       fine[i] = (C2F2_rev_*coarse[std::slice(N-3,3,1)]).sum();
     }
     else if (i == fine.size() - 4) {
       fine[i] = (C2F3_rev_*coarse[std::slice(N-4,4,1)]).sum();
     }
     else if (i == fine.size() - 5) {
       fine[i] = (C2F4_rev_*coarse[std::slice(N-3,3,1)]).sum();
     }
     else if (i == fine.size() - 6) {
       fine[i] = (C2F5_rev_*coarse[std::slice(N-5,5,1)]).sum();
     }
     else if (i == fine.size() - 7) {
       fine[i] = (C2F6_rev_*coarse[std::slice(N-4,4,1)]).sum();
     }
     else if (i == fine.size() - 8) {
       fine[i] = (C2F7_rev_*coarse[std::slice(N-6,6,1)]).sum();
     }
     else if (i == fine.size() - 9) {
       fine[i] = (C2F8_rev_*coarse[std::slice(N-5,5,1)]).sum();
     }
     else if (i == fine.size() - 10) {
       fine[i] = (C2F9_rev_*coarse[std::slice(N-7,7,1)]).sum();
     }
     //interior stencil
     else if (i % 2 == 0) {
       fine[i] = coarse[coarse_ind];
     }
     else {
       fine[i] = (C2F_stencil_*
                   coarse[std::slice(coarse_ind-1,4,1)]).sum();
       coarse_ind ++;
     }
  }
  return fine;
}

void Interp42::SetInterpWeights() {

  // weights Coarse2Fine
  C2F0_ = {0.991574022857554, 0.016851954284894,
          -0.008425977142446};
  C2F0_rev_ = C2F0_;
  std::reverse(std::begin(C2F0_rev_), std::end(C2F0_rev_));

  C2F1_ = {0.406428503705774,0.687142992588452,
          -0.093571496294227};
  C2F1_rev_ = C2F1_;
  std::reverse(std::begin(C2F1_rev_), std::end(C2F1_rev_));

  C2F2_ = {-0.014867566418213,1.029735132836425,
           -0.014867566418212};
  C2F2_rev_ = C2F2_;
  std::reverse(std::begin(C2F2_rev_), std::end(C2F2_rev_));

  C2F3_ = {-0.064519970133486,0.566539940266974,
            0.560480029866513,-0.0625};
  C2F3_rev_ = C2F3_;
  std::reverse(std::begin(C2F3_rev_), std::end(C2F3_rev_));

  C2F4_ = {-0.009401996023082,0.018803992046164,
            0.990598003976917};
  C2F4_rev_ = C2F4_;
  std::reverse(std::begin(C2F4_rev_), std::end(C2F4_rev_));

  C2F5_ = {-0.020335987297944,-0.010109275404113,
            0.518726512702056,0.57421875,-0.0625};
  C2F5_rev_ = C2F5_;
  std::reverse(std::begin(C2F5_rev_), std::end(C2F5_rev_));

  C2F6_ = {-0.048040789592363,0.116914912518061,
           -0.089707456259030,1.020833333333333};
  C2F6_rev_ = C2F6_;
  std::reverse(std::begin(C2F6_rev_), std::end(C2F6_rev_));

  C2F7_ = {-0.013141490707368,0.038001731414735,
           -0.099078990707368,0.57421875,0.5625,-0.0625};
  C2F7_rev_ = C2F7_;
  std::reverse(std::begin(C2F7_rev_), std::end(C2F7_rev_));

  C2F8_ = {0.004986976682269,-0.009973953364539,
           0.004986976682268,0,1.0};
  C2F8_rev_ = C2F8_;
  std::reverse(std::begin(C2F8_rev_), std::end(C2F8_rev_));

  C2F9_ = {0.022698782465675,-0.046699648264684,
           0.025302949132342,-0.063802083333333,
            0.5625,0.5625,-0.0625};
  C2F9_rev_ = C2F9_;
  std::reverse(std::begin(C2F9_rev_), std::end(C2F9_rev_));

  C2F_stencil_ = {-0.0625,0.5625,0.5625,-0.0625};

  //Fine2Coarse
  F2C0_ = {0.495787011428777,0.705272991724725,
           -0.018803098705387,-0.092984662839436,
           -0.013273406150234,-0.028709629126509,
           -0.067822291189219,-0.018552692763343,
           0.007040437669086,0.032045339951541};
  F2C0_rev_ = F2C0_;
  std::reverse(std::begin(F2C0_rev_), std::end(F2C0_rev_));

  F2C1_ = {0.002427823922400,0.343571496294226,
           0.375242463660731,0.235258110788828,
           0.007649081510304,-0.004112247622012,
           0.047558608481923,0.015458331422943,
          -0.004057201368626,-0.018996467090719};
  F2C1_rev_ = F2C1_;
  std::reverse(std::begin(F2C1_rev_), std::end(F2C1_rev_));

  F2C2_ = {-0.001665600132809,-0.064194398620458,
           -0.007433783209106,0.319343272830920,
            0.552891909196419,0.289521774531380,
           -0.050069277912017,-0.055299901790159,
           0.002783428845917,0.014122576259912};
  F2C2_rev_ = F2C2_;
  std::reverse(std::begin(F2C2_rev_), std::end(F2C2_rev_));

  F2C_stencil_ = {-0.03125,0,0.28125,0.5,0.28125, 0, -0.03125};
}
