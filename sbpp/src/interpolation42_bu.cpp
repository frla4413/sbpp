/*
 * interpolation42.cpp
 *
 *  Created on: 18 Jan. 2021
 *      Author: Fredrik Laurén
 */

#include "interpolation42.h"
#include <stdexcept>

Interpolation42::Interpolation42(int Nc, int Nf) : Interpolation(Nc, Nf){}

Array2D Interpolation42::Coarse2Fine(const Array2D& coarse)
{

   Array2D fine(1, (coarse.Size()-1)*2 + 1);
   int coarse_ind = 1;
   int N = coarse.Size();
   for(int i = 0; i < fine.Size(); i++)
   {
      if( i == 0)
      {fine[i] = (b_row1_*coarse[std::slice(0,3,1)]).sum();}

      else if (i == 1)
      {fine[i] = (b_row2_*coarse[std::slice(0,3,1)]).sum();}

      else if (i == 2)
      {fine[i] = (b_row3_*coarse[std::slice(0,4,1)]).sum();}

      else if (i == 3)
      {fine[i] = (b_row4_*coarse[std::slice(0,4,1)]).sum();}

      else if (i == 4)
      {fine[i] = (b_row5_*coarse[std::slice(0,5,1)]).sum();}

      else if (i == 5)
      {fine[i] = (b_row6_*coarse[std::slice(0,5,1)]).sum();}

      else if (i == 6)
      {fine[i] = (b_row7_*coarse[std::slice(0,6,1)]).sum();}

      else if (i == 7)
      {fine[i] = (b_row8_*coarse[std::slice(0,6,1)]).sum();}

      else if (i == 8)
      {fine[i] = (b_row9_*coarse[std::slice(0,7,1)]).sum();}

      else if (i == 9)
      {fine[i] = (b_row10_*coarse[std::slice(0,7,1)]).sum();}

      else if (i == 10)
      {fine[i] = (b_row11_*coarse[std::slice(0,8,1)]).sum();}

      else if (i == fine.Size() - 1)
      {fine[i] = (b_row1_rev_*coarse[std::slice(N-4,3,1)]).sum();}

      else if (i == fine.Size() - 2)
      {fine[i] = (b_row2_rev_*coarse[std::slice(N-4,3,1)]).sum();}

      else if (i == fine.Size() - 3)
      {fine[i] = (b_row3_rev_*coarse[std::slice(N-5,4,1)]).sum();}

      else if (i == fine.Size() - 4)
      {fine[i] = (b_row4_rev_*coarse[std::slice(N-5,4,1)]).sum();}

      else if (i == fine.Size() - 5)
      {fine[i] = (b_row5_rev_*coarse[std::slice(N-6,5,1)]).sum();}

      else if (i == fine.Size() - 6)
      {fine[i] = (b_row6_rev_*coarse[std::slice(N-6,5,1)]).sum();}

      else if (i == fine.Size() - 7)
      {fine[i] = (b_row7_rev_*coarse[std::slice(N-7,6,1)]).sum();}

      else if (i == fine.Size() - 8)
      {fine[i] = (b_row8_rev_*coarse[std::slice(N-7,6,1)]).sum();}

      else if (i == fine.Size() - 9)
      {fine[i] = (b_row9_rev_*coarse[std::slice(N-8,7,1)]).sum();}

      else if (i == fine.Size() - 10)
      {fine[i] = (b_row10_rev_*coarse[std::slice(N-8,7,1)]).sum();}

      else if (i == fine.Size() - 11)
      {fine[i] = (b_row11_*coarse[std::slice(N-9,8,1)]).sum();}

      //interior stencil
      else if (i % 2 == 1)
      {
         fine[i] = (b_stencil1_*coarse[std::slice(coarse_ind-5,4,1)]).sum();
      }
      else
      {
         fine[i] = (b_stencil2_*coarse[std::slice(coarse_ind-6,5,1)]).sum();
         coarse_ind ++;
      }
   }
   return fine;
}

Array2D Interpolation42::Fine2Coarse(const Array2D& fine)
{
   Array2D coarse(1, (fine.Size() + 1)/2);
   int fine_ind = 2;
   int N = fine.Size();

   for(int i = 0; i < coarse.Size(); i++)
   {

      std::cout << a_row1_.size() << " ";
      if( i == 0) {coarse[i] = (a_row1_*fine[std::slice(0,11,1)]).sum();}

      else if (i == 1)
      {coarse[i] = (a_row2_*fine[std::slice(0,11,1)]).sum();}

      else if (i == 2)
      {coarse[i] = (a_row3_*fine[std::slice(0,11,1)]).sum();}

      else if (i == coarse.Size() - 3)
      {coarse[i] = (a_row3_rev_*fine[std::slice(N-12,11,1)]).sum();}

      else if (i == coarse.Size() - 2)
      {coarse[i] = (a_row2_rev_*fine[std::slice(N-12,11,1)]).sum();}

      else if (i == coarse.Size() - 1)
      {coarse[i] = (a_row1_rev_*fine[std::slice(N-12,11,1)]).sum();}

      else
      {
         coarse[i] = (a_stencil_*fine[std::slice(fine_ind-5,9,1)]).sum();
         fine_ind += 2;
      }
   }
   return coarse;
}
