/*
 * interpolation.cpp
 *
 *  Created on: 16 Jan. 2021
 *      Author: Fredrik Laurén
 */

#include "interpolation.h"

Interpolation::Interpolation(int Nc, int Nf): Nc_(Nc), Nf_(Nf){};

Array2D Interpolation::Interpolate(const Array2D &array, int N1, int N2)
{
   if(N1 > N2)
      return Fine2Coarse(array);
   else if(N1 < N2)
      return Coarse2Fine(array);
   else //same size, no need to interpolate
      return array;
}
