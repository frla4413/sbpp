//============================================================
// Name        : interp.cpp
// Date        : May 2022
// Author      : Fredrik Laurén
// Description : source file for the class interp.
//============================================================


#include "interp.hpp"

Interp::Interp(int Nc, int Nf) {
  if ((static_cast<double>(Nf-1) / static_cast<double>(Nc-1)) !=2)
    throw std::invalid_argument("@Interp - (Nf-1)/(Nc-1) = 2 is required!");

  Nf_ = Nf;
  Nc_ = Nc;
}

Array Interp::Interpolate(const Array &array) {
  if((Nf_ == Nc_) && array.size() == Nc_)
    return array; //same size, no need to interpolate
  else if (array.size() == Nf_)
    return Fine2Coarse(array);
  else if (array.size() == Nc_)
    return Coarse2Fine(array);
  else
    throw std::invalid_argument("@Interpolation - wrong size!");
}
