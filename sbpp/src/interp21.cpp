//============================================================
// Name        : interp21.cpp
// Date        : October 2022
// Author      : Fredrik Laurén
// Description : source file for the class interp21.
//============================================================

#include "interp21.hpp"
#include <stdexcept>

Interp21::Interp21(int Nc, int Nf) : Interp(Nc, Nf){}

Array Interp21::Coarse2Fine(const Array& coarse) {

  Array fine(1, (coarse.size()-1)*2 + 1);
  int coarse_ind = 0;
  for(int i = 0; i < fine.size(); ++i) {
     if(i == fine.size() - 1) {
       // last element out of bounds if not special treatment
       fine[fine.size()-1] = coarse[coarse.size()-1];
     }
     else if( i % 2 == 0) {
       fine[i] = coarse[coarse_ind];
     }
     else {
       ++coarse_ind;
       fine[i] = 0.5*(coarse[coarse_ind] + coarse[coarse_ind-1]);
     }
  }
  return fine;
}

Array Interp21::Fine2Coarse(const Array& fine) {
  Array coarse(1, (fine.size() + 1)/2);
  int fine_ind = 2;

  for(int i = 0; i < coarse.size(); ++i) {
    if( i == 0) {
      coarse[i] = 0.5*(fine[0] + fine[1]);
    }
    else if (i == coarse.size() - 1) {
      int N = fine.size();
      coarse[i] = 0.5*(fine[N-2] + fine[N-1]);
    }
    else {
      coarse[i] = 0.25*(fine[fine_ind-1] +
                        fine[fine_ind+1]) + 
                  0.5*fine[fine_ind];
      fine_ind += 2;
    }
  }
  return coarse;
}
