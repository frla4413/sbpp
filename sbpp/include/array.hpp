//============================================================
// Name        : array.hpp
// Date        : May 2022
// Author      : Fredrik Laurén
// Description : hpp file for the class array.
//============================================================

/*
 *  1D array with (a few) 2D features.
 *  Uses std::valarray.
 *  The elements are ordered column by column.
 *  Nx - number of points in x-direction --> number of columns
 *  Ny - number of points in y-direction --> number of rows
 *
 * Member functions:
 *  o size()
 *  o Nx(), Ny()
 *  o arrray()
 *  o begin(), cbegin(), end(), cend()
 *  o Reverse()
 *
 * Member variables:
 *  o Nx_, Ny_, size_ = Nx*Ny
 *  o array_ -- valarray with doubles
 *
 * Operator functions:
 *  o +=, -=, +, -, *, / -- element wise operations
 *  o [std::slice], [int] -- access value(s)
 *
 * Non-member functions:
 *  o Print(Array)
 *  o Shape(Array) --> (Nx, Ny)
 *  o Sum(array)
 *  o Max(array)
 *  o Min(array)
 *  o Zeros(Nx, Ny)
 *  o Ones(Nx, Ny)
 */

#pragma once
#include <functional>
#include <utility>
#include <valarray>
#include <iostream>

struct Shape {
  int Nx, Ny;
  bool operator !=(const Shape& other) const;
};

void Print(Shape shape);

enum class Operation {add, sub, div, mul};

class Array {

  //alias to improve readability
  using valarray = std::valarray<double>;

  public:

    Array();
    Array(int Nx, int Ny);
    Array(int Nx, int Ny, const valarray& array);

    int Nx() const;
    int Ny() const;
    valarray& array();
    const valarray& array() const;
    int size() const;

    double* begin();
    double* end();

    const double* cbegin();
    const double* cend();

    void Reverse();

// --------------- operator functions ---------------------------

    const double& operator[] (const int pos) const;
    const valarray operator[] (const std::slice index) const;

    double& operator[] (const int pos);
    std::slice_array<double> operator[] (const std::slice index);

    Array& operator+= (const Array& x);
    Array& operator-= (const Array& x);

    Array& operator+= (double a);
    Array& operator-= (double a);

    friend Array operator+ (const Array& x, const Array& y);
    friend Array operator- (const Array& x, const Array& y);
    friend Array operator* (const Array& x, const Array& y);
    friend Array operator/ (const Array& x, const Array& y);

    friend Array operator+ (const Array& x, double a);
    friend Array operator+ (double a, const Array& x);
    friend Array operator- (const Array& x, double a);
    friend Array operator- (double a, const Array& x);
    friend Array operator* (double a, const Array& x);
    friend Array operator* (const Array& x, double a);
    friend Array operator/ (double a, const Array& x);
    friend Array operator/ (const Array& x, double a);

  private:
    int Nx_ = 0; int Ny_ = 0; int size_ = 0;
    valarray array_;

    friend void CheckDimensions(const Array& x, const Array& y);
    friend Array VectorOperation(const Array& x, const Array& y,
                                 Operation op);
};

/*
 * Print the array in the terminal.
 */
void Print(const Array& array);

/*
 * Generate Array with 0/1 of size (Nx, Ny).
 */
Array Zeros(int Nx, int Ny);
Array Ones(int Nx, int Ny);

/*
 * Get the shape (Nx, Ny) of array.
 */
Shape GetShape(const Array& array);

/*
 * Compute Max, Min, Sum of array
 */
double Max(const Array& array);
double Min(const Array& array);
double Sum(const Array& array);

/*
 * Return element-wise abs of array.
 */
Array Abs(const Array& array);
