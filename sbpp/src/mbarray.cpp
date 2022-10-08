//============================================================
// Name        : mbarray.cpp
// Date        : May 2022
// Author      : Fredrik Laurén
// Description : source file for the class mbarray.
//============================================================

#include <iomanip>
#include "mbarray.hpp"

//alias to improve readability
using valarray = std::valarray<double>;
using Shapes = std::vector<Shape>;
using Arrays = std::vector<Array>;

MbArray::MbArray(const Shapes& shapes, const Arrays& arrays) {

  if(shapes.size() != arrays.size())
    throw std::invalid_argument("@MbArray invalid sizes");

  auto shapes_it = shapes.cbegin();
  auto arrays_it = arrays.cbegin();
  for(int i = 0; i < shapes.size(); ++i) {

    if((*shapes_it).Nx*(*shapes_it).Ny != (*arrays_it).size())
      throw std::invalid_argument("@MbArray invalid shape");
    ++shapes_it;
    ++arrays_it;
  }
  shapes_ = shapes;
  arrays_ = arrays;
}

MbArray::MbArray(const Arrays& arrays) {

  std::vector<Shape> shapes(arrays.size());
  for(int i = 0; i < arrays.size(); ++i) {
    shapes[i] = GetShape(arrays[i]);
  }
  shapes_ = shapes;
  arrays_ = arrays;
}

Arrays& MbArray::arrays() {
  return arrays_;
}

const Arrays& MbArray::arrays() const {
  return arrays_;
}

Array& MbArray::array(int idx) {
  return arrays_[idx];
}

const Shapes& MbArray::shapes() const {
  return shapes_;
}

const Shape& MbArray::shapes(int idx) const{
  return shapes_[idx];
}

int MbArray::Size() const {
  return shapes_.size();
}

std::vector<Array>::iterator MbArray::begin() {
  return arrays_.begin();
}

std::vector<Array>::iterator MbArray::end() {
  return arrays_.end();
}

std::vector<Array>::const_iterator MbArray::cbegin() const {
  return arrays_.cbegin();
}

std::vector<Array>::const_iterator MbArray::cend() const {
  return arrays_.cend();
}

Array MbArray::operator[] (const MbSlice slice) const {

  auto vals = arrays_[slice.block_idx][slice.slice];
  return {shapes_[slice.block_idx].Nx,
          shapes_[slice.block_idx].Ny, vals};
}

Array& MbArray::operator[] (int idx) {
  return arrays_[idx];
}

const Array& MbArray::operator[] (int idx) const {
  return arrays_[idx];
}

void CheckSize(const MbArray& x, const MbArray& y) {
  if (x.Size() != y.Size())
    throw std::invalid_argument("@MbArray: Size Error");
}

MbArray VectorOperation(const MbArray& x, const MbArray& y,
                        Operation op) {
  try {
    CheckSize(x, y);
  }
  catch (std::invalid_argument) {
    std::cout << "@MbArray VectorOperatation Size Error:";
    std::cout << "x.Size(): " << x.Size() << std::endl;
    std::cout << "y.Size(): " << y.Size() << std::endl;
    exit(1);
  }

  Arrays z(x.Size());

  switch (op) {
    case Operation::add :
      for(int i = 0; i < x.Size(); ++i)
        z[i] = x[i] + y[i];
      break;
    case Operation::sub :
      for(int i = 0; i < x.Size(); ++i)
        z[i] = x[i] - y[i];
      break;
    case Operation::mul :
      for(int i = 0; i < x.Size(); ++i)
        z[i] = x[i] * y[i];
      break;
    case Operation::div :
      for(int i = 0; i < x.Size(); ++i)
        z[i] = x[i] / y[i];
      break;
  }

  return {x.shapes(), z};
}

// ----------------- operator functions -------------------------
MbArray& MbArray::operator+=(double a) {
  for(auto & array : this->arrays_)
    array+=a;
  return *this;
}

MbArray& MbArray::operator-=(double a) {
  for(auto & array : this->arrays_)
    array-=a;
  return *this;
}

MbArray operator+ (const MbArray& x, const MbArray& y) {
  return VectorOperation(x,y,Operation::add);
}

MbArray operator- (const MbArray& x, const MbArray& y) {
  return VectorOperation(x,y,Operation::sub);
}

MbArray operator* (const MbArray& x, const MbArray& y) {
  return VectorOperation(x,y,Operation::mul);
}

MbArray operator/ (const MbArray& x, const MbArray& y) {
  return VectorOperation(x,y,Operation::div);
}

MbArray ScalarOperation(const MbArray& x, double a,
                        Operation op) {
  Arrays z(x.Size());
  switch (op) {
    case Operation::add :
      for(int i = 0; i < x.Size(); ++i)
        z[i] = x[i] + a;
      break;
    case Operation::sub :
      for(int i = 0; i < x.Size(); ++i)
        z[i] = x[i] - a;
      break;
    case Operation::mul :
      for(int i = 0; i < x.Size(); ++i)
        z[i] = x[i] * a;
      break;
    case Operation::div :
      for(int i = 0; i < x.Size(); ++i)
        z[i] = x[i] / a;
      break;
  }
  return {x.shapes(), z};
}

MbArray operator+ (const MbArray& x, double a) {
  Arrays z(x.Size());
  for(int i = 0; i < x.Size(); ++i)
    z[i] = x[i] + a;
  return {x.shapes(), z};
}

MbArray operator+ (double a, const MbArray& x) {
  Arrays z(x.Size());
  for(int i = 0; i < x.Size(); ++i)
    z[i] = a + x[i];
  return {x.shapes(), z};
}

MbArray operator- (const MbArray& x, double a) {
  Arrays z(x.Size());
  for(int i = 0; i < x.Size(); ++i)
    z[i] = x[i] - a;
  return {x.shapes(), z};
}

MbArray operator- (double a, const MbArray& x) {
  Arrays z(x.Size());
  for(int i = 0; i < x.Size(); ++i)
    z[i] = a - x[i];
  return {x.shapes(), z};
}

MbArray operator* (const MbArray& x, double a) {
  Arrays z(x.Size());
  for(int i = 0; i < x.Size(); ++i)
    z[i] = x[i] * a;
  return {x.shapes(), z};
}

MbArray operator* (double a, const MbArray& x) {
  Arrays z(x.Size());
  for(int i = 0; i < x.Size(); ++i)
    z[i] = a * x[i];
  return {x.shapes(), z};
}

MbArray operator/ (const MbArray& x, double a) {
  Arrays z(x.Size());
  for(int i = 0; i < x.Size(); ++i)
    z[i] = x[i] / a;
  return {x.shapes(), z};
}

MbArray operator/ (double a, const MbArray& x) {
  Arrays z(x.Size());
  for(int i = 0; i < x.Size(); ++i)
    z[i] = a / x[i];
  return {x.shapes(), z};
}

void Print(const MbArray& mb_array) {

  for(auto& array : mb_array.arrays()) {
    Print(array);
    std::cout << "\n";
  }
}

double Min(const MbArray mb_array) {

  auto min = mb_array[0][0];

  for(auto& array : mb_array.arrays()) {
    auto array_min = Min(array);
    if(array_min < min)
      min = array_min;
  }

  return min;
}

double Max(const MbArray mb_array) {

  auto max = mb_array[0][0];

  for(auto& array : mb_array.arrays()) {
    auto array_max = Max(array);
    if(array_max > max)
      max = array_max;
  }
  return max;
}

double Sum(const MbArray mb_array) {

  auto sum = 0.0;

  for(auto& array : mb_array.arrays()) {
    sum += Sum(array);
  }
  return sum;
}
