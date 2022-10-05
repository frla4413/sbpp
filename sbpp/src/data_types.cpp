#include "data_types.h"

namespace DataTypes
{
   // Block
   Block SetBlock(const Array2D x, const Array2D y)
   {
      Block block;
      block.x = x;
      block.y = y;
      return block;
   }
}
