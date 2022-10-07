//============================================================
// Name        : interp.hpp
// Date        : October 2022
// Author      : Fredrik Laurén
// Description : hpp file for the class interp.
//============================================================
/*
 * Description:
 * Base class for interpolation operators.
 * Use these operators to interpolate between a
 * fine and course grid.
 * The relation between fine and course must be 2:1.
 * The operators are taken from
 * "Mattsson and Karpenter,
 * Stable and Accurate Interpolation Operators
 * for High-Order Multi-Block Finite-Difference Methods"
 * The function Interpolate interpolates an array from grid
 * of size N1 to grid of size N2, using
 * either Course2Fine or Fine2Course.
 *
 * Member functions:
 *  o Interpolate(array, N1, N2)
 *  o Coarse2Fine(array) (protected, virtual)
 *  o Fine2Course(array) (protected, virtual)
 *
 * Member variables:
 *  o Nc_ - N course
 *  o Nf_ - N fine
 *
 * Comments:
 *  o
*/

#pragma once
#include <iostream>
#include "array.hpp"

class Interp {

  public:
    Interp(){};
    Interp(int Nc, int Nf);
    virtual ~Interp(){};
    Array Interpolate(const Array& array);

  protected:
     virtual Array Coarse2Fine(const Array& coarse) = 0;
     virtual Array Fine2Coarse(const Array& fine) = 0;
     int Nc_, Nf_;
};
