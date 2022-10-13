//============================================================
// Name        : interp.cpp
// Date        : May 2022
// Author      : Fredrik Laurén
// Description : source file for the class interp.
//============================================================


#include "interp.hpp"
#include <string.h>

Interp::Interp(int N1, int N2) {

  Nf_ = N1;
  Nc_ = N2;
  if(N1 < N2) {
    Nc_ = N1;
    Nf_ = N2;
  }
  if (Nc_ != Nf_ && ((Nf_-1) / (Nc_-1)) !=2) {
    int quotient = (Nf_ - 1)/(Nc_-1);
    std::string msg = "@Interp - (Nf-1)/(Nc-1) = " +
      std::to_string(quotient) + " != 2 \n";
    throw std::invalid_argument(msg);
  }
}

Array Interp::Interpolate(const Array &array) {
  if((Nf_ == Nc_) && array.size() == Nc_)
    return array; //same size, no need to interpolate
  else if (array.size() == Nf_)
    return Fine2Coarse(array);
  else if (array.size() == Nc_)
    return Coarse2Fine(array);
  else {
    std::string msg = "@Interpoation - wrong size: Nf_ = "
                      + std::to_string(Nf_) +
                      " Nc_ " + std::to_string(Nc_)  +
                      " array.size() = " +
                      std::to_string(array.size());
    //throw std::invalid_argument("@Interpolation - wrong size!");
    throw std::invalid_argument(msg);
  }
}
