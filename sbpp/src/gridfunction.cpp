/*
 * gridfunction2D.cpp
 *
 *  Created on: September 2021
 *      Author: Fredrik Laurén
 */

#include "gridfunction.h"
#include <iostream>
#include <assert.h>

Gridfunction::Gridfunction(){}

Gridfunction::Gridfunction(int num_blocks) : 
   gridfunction_(num_blocks), 
   shapes_(std::vector<std::pair<int,int>>(num_blocks)),
   num_blocks_(num_blocks)
{}

Gridfunction::Gridfunction(const std::vector<Array2D>&  vector)
{
   gridfunction_ = vector;
   for(auto& array2D : gridfunction_)
      shapes_.push_back(array2D.GetShape());
   num_blocks_ = gridfunction_.size();
}


Gridfunction::Gridfunction(const std::valarray<double>& f,
                           std::vector<std::pair<int,int>> shapes):
   shapes_(shapes), num_blocks_(shapes.size())
{
   int start = 0; 
   for(auto& shape : shapes_)
   {
      int size = shape.first*shape.second;
      gridfunction_.push_back(Array2D(shape.first,shape.second,
                              f[std::slice(start,size,1)]));
      start += size;
   }
}

int Gridfunction::GetNumBlocks() 
{ 
   return num_blocks_;
};

int Gridfunction::GetNumBlocks() const 
{ 
   return num_blocks_;
};

std::pair<int, int> Gridfunction::GetShape(int block_idx) const
{
   return shapes_[block_idx];
}

std::vector<std::pair<int, int>> Gridfunction::GetShapes() const 
{
   return shapes_;
};

int Gridfunction::GetTotalSize() const
{
   int tot_size = 0;
   for(auto& shape : shapes_)
      tot_size += shape.first * shape.second;

   return tot_size;
}


std::valarray<double> Gridfunction::GridfunctionToValarray() const
{
   std::valarray<double> array(GetTotalSize());

   int start = 0;
   for(auto& array2d : gridfunction_)
   {
      int size = array2d.Size();
      array[std::slice(start,size,1)] = array2d.GetArray();
      start += size;
   }
   return array;
}


void Gridfunction::Print() const
{
   for(auto& array2D : gridfunction_)
   {
      array2D.Print();
      std::cout<< std::endl;
   }
}

double Gridfunction::Sum()
{
   double sum = 0; 
   for(auto& array2D : gridfunction_)
      sum += array2D.Sum();
   return sum;
}

double Gridfunction::Max()
{
   double max_val = -10, max_candidate; 
   for(auto& array2D : gridfunction_)
   {
      max_candidate = array2D.Max();
      if(max_candidate > max_val)
      {
         max_val = max_candidate;
      }
   }
   return max_val;
}

double Gridfunction::Min()
{
   double min_val = 10, min_candidate; 
   for(auto& array2D : gridfunction_)
   {
      min_candidate = array2D.Max();
      if(min_candidate < min_val)
      {
         min_val = min_candidate;
      }
   }
   return min_val;
}

void Gridfunction::PushBack(Array2D arr)
{
   gridfunction_.push_back(arr);
   shapes_.push_back(arr.GetShape());
   num_blocks_ = gridfunction_.size();
}

std::vector<Array2D>::iterator Gridfunction::begin() 
{
   return gridfunction_.begin();
};

std::vector<Array2D>::iterator Gridfunction::begin() const
{
   return begin();
};

std::vector<Array2D>::iterator Gridfunction::end() 
{
   return gridfunction_.end();
};

std::vector<Array2D>::iterator Gridfunction::end() const
{
   return end();
};

// --------------------- operator functions --------------------------------------
Gridfunction& Gridfunction::operator+=(const Gridfunction&x)
{
   assert(this->num_blocks_ == x.num_blocks_);
   for(int i = 0; i < num_blocks_; i++)
      this->gridfunction_ [i]+= x.gridfunction_[i];
   this->shapes_ = x.shapes_;
   this->num_blocks_ = x.num_blocks_;
   return *this;
}

Gridfunction& Gridfunction::operator-=(const Gridfunction&x)
{
   assert(this->num_blocks_ == x.num_blocks_);
   for(int i = 0; i < num_blocks_; i++)
      this->gridfunction_ [i]-= x.gridfunction_[i];
   this->shapes_ = x.shapes_;
   this->num_blocks_ = x.num_blocks_;
   return *this;
}

Gridfunction operator+ (const Gridfunction& x, const Gridfunction& y)
{
   assert(x.num_blocks_ == y.num_blocks_);
   std::vector<Array2D> sum(x.num_blocks_);
   for(int i = 0; i < x.num_blocks_; i++)
      sum[i] = x.gridfunction_[i] + y.gridfunction_[i];
   return Gridfunction(sum);
}

Gridfunction operator+ (const Gridfunction& x, const double a)
{

   std::vector<Array2D> sum(x.num_blocks_);
   for(int i = 0; i < x.num_blocks_; i++)
      sum[i] = x.gridfunction_[i] + a;
   return Gridfunction(sum);
}

Gridfunction operator+ (const double a, const Gridfunction& x)
{
    return x + a;
}

Gridfunction operator- (const Gridfunction& x, const Gridfunction& y)
{
   assert(x.num_blocks_ == y.num_blocks_);
   std::vector<Array2D> diff(x.num_blocks_);
   for(int i = 0; i < x.num_blocks_; i++)
      diff[i] = x.gridfunction_[i] - y.gridfunction_[i];
   return Gridfunction(diff);
}

Gridfunction operator- (const Gridfunction& x, const double a)
{

   std::vector<Array2D> diff(x.num_blocks_);
   for(int i = 0; i < x.num_blocks_; i++)
      diff[i] = x.gridfunction_[i] - a;
   return Gridfunction(diff);
}

Gridfunction operator- (const double a, const Gridfunction& x)
{
    return x - a;
}

Gridfunction operator* (const Gridfunction& x, const Gridfunction& y)
{
   assert(x.num_blocks_ == y.num_blocks_);
   std::vector<Array2D> prod(x.num_blocks_);
   for(int i = 0; i < x.num_blocks_; i++)
      prod[i] = x.gridfunction_[i] * y.gridfunction_[i];
   return Gridfunction(prod);
}

Gridfunction operator* (const Gridfunction& x, const double a)
{

   std::vector<Array2D> prod(x.num_blocks_);
   for(int i = 0; i < x.num_blocks_; i++)
      prod[i] = x.gridfunction_[i] * a;
   return Gridfunction(prod);
}

Gridfunction operator* (const double a, const Gridfunction& x)
{
    return x * a;
}

Gridfunction operator/ (const Gridfunction& x, const Gridfunction& y)
{
   assert(x.num_blocks_ == y.num_blocks_);
   std::vector<Array2D> quoti(x.num_blocks_);
   for(int i = 0; i < x.num_blocks_; i++)
      quoti[i] = x.gridfunction_[i] / y.gridfunction_[i];
   return Gridfunction(quoti);
}

Gridfunction operator/ (const Gridfunction& x, const double a)
{

   std::vector<Array2D> quoti(x.num_blocks_);
   for(int i = 0; i < x.num_blocks_; i++)
      quoti[i] = x.gridfunction_[i] / a;
   return Gridfunction(quoti);
}

Gridfunction operator/ (const double a, const Gridfunction& x)
{
    return x / a;
}

Array2D& Gridfunction::operator[] (const int block_idx) 
{ 
   return gridfunction_[block_idx];
};
const Array2D& Gridfunction::operator[] (const int block_idx) const 
{ 
   return gridfunction_[block_idx];
};

