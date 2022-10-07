//============================================================
// Name        : array.cpp
// Date        : May 2022
// Author      : Fredrik Laurén
// Description : source file for the class array.
//============================================================

#include <iomanip>
#include "array.hpp"

//alias to improve readability
using valarray = std::valarray<double>;

bool Shape::operator !=(const Shape& other) const {
  return Nx != other.Nx && Ny != other.Ny;
}

void Print(Shape shape) {
  std::cout << "(" << shape.Nx << "," << shape.Ny << ")";
}

Array Ones(int Nx, int Ny) {
  valarray one_vec(1,Nx*Ny);
  return {Nx, Ny, one_vec};
};

Array Zeros(int Nx, int Ny) {
  valarray zero_vec(Nx*Ny);
  return {Nx, Ny, zero_vec};
};

Array::Array() :
  Nx_(0), Ny_(0), array_() {}

Array::Array(int Nx, int Ny) :
  Nx_(Nx), Ny_(Ny), size_(Nx*Ny), array_(valarray(Nx*Ny)) {}

Array::Array(int Nx, int Ny, const valarray& array) :
  Nx_(Nx), Ny_(Ny), size_(array.size()) {

  if(Nx*Ny != array.size())
    throw std::invalid_argument("@Array(Nx,Ny,x) Nx*Ny != x.size");
  array_ = array;
}

int Array::Nx() const {
  return Nx_;
};

int Array::Ny() const {
  return Ny_;
};
valarray& Array::array() {
  return array_;
};

const valarray& Array::array() const {
  return array_;
};

int Array::size() const {
  return size_;
}

Shape GetShape(const Array& array) {
  return {array.Nx(), array.Ny()};
};

void Print(const Array& array) {
  int Nx = array.Nx();
  int Ny = array.Ny();
  valarray row(Nx);
  for(int i = 0; i < Ny; i++) {
    row = array.array()[std::slice(i, Nx, Ny)];
    for (auto& val : row)
      std::cout << std::fixed << val << "  "
                << std::setw(9)
                << std::setprecision(3);
    std::cout << "\n";
  }
}

void CheckDimensions(const Array& x, const Array& y) {
  if (x.Nx_ != y.Nx_ || x.Ny_ != y.Ny_)
    throw std::invalid_argument("@Array: CheckDimensions");
}

double* Array::begin() {
  return std::begin(array_);
};

double* Array::end() {
  return std::end(array_);
};

const double* Array::cbegin() {
  return std::cbegin(array_);
};

const double* Array::cend() {
  return std::cend(array_);
};

void Array::Reverse() {
  std::reverse(std::begin(array_),std::end(array_));
}

double Max(const Array& array) {
  return array.array().max();
};

double Min(const Array& array) {
  return array.array().min();
};

double Sum(const Array& array) {
  return array.array().sum();
};

Array Abs(const Array& array) {
  return {array.Nx(), array.Ny(), std::abs(array.array())};
};

// ----------------- operator functions -------------------------
Array& Array::operator+=(const Array& x) {
  try {
    CheckDimensions(*this, x);
  }
  catch (std::invalid_argument) {
    std::cout << "Dimension error in Array: += \n";
  }

  this->array_ += x.array_;
  return *this;
}

Array& Array::operator-=(const Array& x) {
  try {
    CheckDimensions(*this, x);
  }
  catch (std::invalid_argument) {
    std::cout << "Dimension error in Array: -= \n";
  }

  this->array_ -= x.array_;
  return *this;
}

Array VectorOperation (const Array& x, const Array& y,
                      Operation op) {
  try {
    CheckDimensions(x, y);
  }
  catch (std::invalid_argument) {
    std::cout << "Dimension error in VectorOperation: \n";
    std::cout << "GetShape(x): ";
    Print(GetShape(x));

    std::cout << "\nGetShape(y): ";
    Print(GetShape(y));
    exit(1);
  }

  valarray z;
  switch (op) {
    case Operation::add :
      z = x.array_ + y.array_;
      break;
    case Operation::sub :
      z = x.array_ - y.array_;
      break;
    case Operation::mul :
      z = x.array_ * y.array_;
      break;
    case Operation::div :
      z = x.array_ / y.array_;
      break;
  }

  return {x.Nx_, x.Ny_, z};
}

Array operator+ (const Array& x, const Array& y) {
  return VectorOperation(x,y,Operation::add);
}

Array operator- (const Array& x, const Array& y) {
  return VectorOperation(x,y,Operation::sub);
}

Array operator* (const Array& x, const Array& y) {
  return VectorOperation(x,y,Operation::mul);
}

Array operator/ (const Array& x, const Array& y) {
  return VectorOperation(x,y,Operation::div);
}

Array& Array::operator+=(double a) {
  this->array_ += a;
  return *this;
}

Array& Array::operator-=(double a) {
  this->array_ -= a;
  return *this;
}

Array operator+ (const Array& x, double a) {
  return {x.Nx_, x.Ny_, x.array_ + a};
}

Array operator+ (double a, const Array& x) {
  return {x.Nx_, x.Ny_, a + x.array_};
}

Array operator- (const Array& x, double a) {
  return {x.Nx_, x.Ny_, x.array_ - a};
}

Array operator- (double a, const Array& x) {
  return {x.Nx_, x.Ny_, a - x.array_};
}

Array operator* (const Array& x, double a) {
  return {x.Nx_, x.Ny_, x.array_ * a};
}

Array operator* (double a, const Array& x) {
  return {x.Nx_, x.Ny_, a * x.array_ };
}

Array operator/ (const Array& x, double a) {
  return {x.Nx_, x.Ny_, x.array_ / a};
}

Array operator/ (double a, const Array& x) {
  return {x.Nx_, x.Ny_, a / x.array_ };
}

const double& Array::operator[] (const int pos) const {
  return array_[pos];
};

const valarray Array::operator[] (const std::slice idx) const {
  return array_[idx];
};

double& Array::operator[] (const int pos) {
  return array_[pos];
};

std::slice_array<double> Array::operator[](const std::slice idx) {
  return array_[idx];
};
