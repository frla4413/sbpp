/*
 * gridfunction.h
 *
 *  Created on: September 2021
 *      Author: Fredrik Laurén
 *
 *  A Gridfunction is a functions that lives on the grid.
 *  It consists of a vector, where each element holds 
 *  the part of the function
 *  that lives on the specific grid block.
 *
 */

#pragma once
#include <functional>
#include <utility>
#include "array.hpp"

class Gridfunction
{
  public:

     Gridfunction();
     Gridfunction(int num_blocks);
     Gridfunction(const std::vector<Array2D>&  vector);
     Gridfunction(const std::valarray<double>& f,std::vector<std::pair<int,int>> shapes);
     ~Gridfunction(){};

     std::pair<int, int> GetShape(int block_idx) const;
     std::vector<std::pair<int, int>> GetShapes() const;

     int GetTotalSize() const;
     std::valarray<double> GridfunctionToValarray() const;

     double Sum();
     double Max();
     double Min();

     int GetNumBlocks();
     int GetNumBlocks() const;
     void Print() const;
     void PushBack(Array2D arr);

     std::vector<Array2D>::iterator begin();
     std::vector<Array2D>::iterator begin() const;
     std::vector<Array2D>::iterator end();
     std::vector<Array2D>::iterator end() const;
     
     // -------------------- operator functions --------------------------------------
     Gridfunction& operator+= (const Gridfunction& x);
     Gridfunction& operator-= (const Gridfunction& x);
     friend Gridfunction operator+  (const Gridfunction& x, const Gridfunction& y);
     friend Gridfunction operator+  (const Gridfunction& x, const double a);
     friend Gridfunction operator+  (const double a, const Gridfunction& x);
     friend Gridfunction operator-  (const Gridfunction& x, const Gridfunction& y);
     friend Gridfunction operator-  (const Gridfunction& x, const double a);
     friend Gridfunction operator-  (const double a, const Gridfunction& x);
     friend Gridfunction operator* (const Gridfunction& x, const Gridfunction& y);
     friend Gridfunction operator* (const double a, const Gridfunction& x);
     friend Gridfunction operator* (const Gridfunction& x, const double a);
     friend Gridfunction operator/ (const Gridfunction& x, const Gridfunction& y);
     friend Gridfunction operator/ (const double a, const Gridfunction& x);
     friend Gridfunction operator/ (const Gridfunction& x, const double a);

     Array2D& operator[] (const int block_idx);
     const Array2D& operator[] (const int block_idx) const;
     
  private:
     std::vector<Array2D> gridfunction_;
     std::vector<std::pair<int,int>> shapes_;
     int num_blocks_;
};
