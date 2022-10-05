/*
 * interpolation21.cpp
 *
 *  Created on: 16 Jan. 2021
 *      Author: Fredrik Laurén
 *      The operators are taken from Tomas Lunquist matlab-code sbpgeom, SBP_HREF.m
 */

#include "interpolation21.h"
#include <stdexcept>

Interpolation21::Interpolation21(int Nc, int Nf) : Interpolation(Nc, Nf){}

Array2D Interpolation21::Coarse2Fine(const Array2D& coarse)
{

   Array2D fine(1, (coarse.Size()-1)*2 + 1);
   int coarse_ind = 0;
   for(int i = 0; i < fine.Size(); i++)
   {
      if(i == fine.Size() - 1) // last element out of bounds if not special treatment
      {
         fine[fine.Size()-1] = coarse[coarse.Size()-1];
      }
      else if( i % 2 == 0)
      {
         fine[i] = coarse[coarse_ind];
      }
      else
      {
         coarse_ind++;
         fine[i] = 0.5*(coarse[coarse_ind] + coarse[coarse_ind-1]);
      }
   }
   return fine;
}

Array2D Interpolation21::Fine2Coarse(const Array2D& fine)
{
   Array2D coarse(1, (fine.Size() + 1)/2);
   int fine_ind = 2;

   for(int i = 0; i < coarse.Size(); i++)
   {
      if( i == 0)
      {
         coarse[i] = 0.5*(fine[0] + fine[1]);
      }
      else if (i == coarse.Size() - 1)
      {
         int N = fine.Size();
         coarse[i] = 0.5*(fine[N-2] + fine[N-1]);
      }
      else
      {
         coarse[i] = 0.25*(fine[fine_ind-1] + fine[fine_ind+1]) + 0.5*fine[fine_ind];
         fine_ind += 2;
      }
   }
   return coarse;
}
