#pragma once
#include <iostream>
#include "gridfunction.h"
#include <vector>
#include "interpolation.h"
#include "array2D.h"

namespace DataTypes{

   //Block
   struct Block{ Array2D x; Array2D y; };
   Block SetBlock(const Array2D X, const Array2D Y);

   std::vector<DataTypes::Corner> GetCorners(const Block&);

   struct InsState{ Gridfunction u,v,p; };

   struct InsSolution
   {
      InsState state;
      double time;
   };

   enum class BdType 
   {
      Wall, SlipWall, Inflow, Pressure, Outflow
   };
}
