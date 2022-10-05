/*
 * interpolation21.h
 *
 *  Created on: 16 Jan. 2021
 *      Author: Fredrik Laurén
 */

#pragma once
#include <iostream>
#include "array2D.h"
#include "interpolation.h"

class Interpolation21:public Interpolation
{
  public:
     Interpolation21(){};
     Interpolation21(int Nc, int Nf);
     ~Interpolation21(){};

  private:
     Array2D Coarse2Fine(const Array2D& coarse);
     Array2D Fine2Coarse(const Array2D& fine);
};
