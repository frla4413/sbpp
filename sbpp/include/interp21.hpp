//============================================================
// Name        : interp21.hpp
// Date        : October 2022
// Author      : Fredrik Laurén
// Description : hpp file for the class interp.
//============================================================
/*
 * Description
 * (2,1) order interpolation operator. Extends Interp.
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
#include <iostream>
#include "array.hpp"
#include "interp.hpp"

class Interp21:public Interp {
 public:
   Interp21(){};
   Interp21(int Nc, int Nf);
   ~Interp21(){};

 private:
   Array Coarse2Fine(const Array& coarse);
   Array Fine2Coarse(const Array& fine);
};
