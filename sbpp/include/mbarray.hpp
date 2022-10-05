//============================================================
// Name        : mbarray.hpp
// Date        : May 2022
// Author      : Fredrik Laurén
// Description : hpp file for the class mbarray.
//============================================================

/*
 *  Multi-block array (MbArray) that extends Array.
 *  The elements are ordered column by column in each block.
 *  Nb - number of blocks
 *  Nx_i - points in x-direction of block i
 *  Ny_i - points in y-direction of block i
 *  Total number of elements = sum_i Nx_i*Ny_i, where i = 0, Nb-1
 *
 * Member functions:
 *
 * Member variables:
 *  o std::vector<int>Nx_, Ny_
 *  o std::vector<Shape> shapes_
 *  o arrays_ -- vector of Array
 *
 * Operator functions:
 *  o +=, -=, +, -, *, / -- element wise operations
 *  o [MbSlice] = [int block_idx, slice] -- access values
 *  o [int i] -- access Array i
 *
 * Non-member functions:
 *  o Print(MbArray)
 *  o Sum(MbArray)
 *  o Min(MbArray)
 *  o Max(MbArray)
 */

#pragma once
#include <valarray>
#include "array.hpp"

struct MbSlice {
  int block_idx;
  std::slice slice;
};

class MbArray {

  //alias to improve readability
  using Shapes = std::vector<Shape>;
  using Arrays = std::vector<Array>;

  public:
    MbArray(const Shapes& shapes, const Arrays& arrays);
    MbArray(const Arrays& arrays);
    MbArray(){};

    Arrays& arrays();
    const Arrays& arrays() const;

    Array& array(int idx);
    const Shapes& shapes() const;
    const Shape& shapes(int idx) const;
    int Size() const;

    std::vector<Array>::iterator begin();
    std::vector<Array>::iterator end();

    std::vector<Array>::const_iterator cbegin() const;
    std::vector<Array>::const_iterator cend() const;

// --------------- operator functions ---------------------------
    Array operator[] (const MbSlice slice) const;
    const Array& operator[] (int idx) const;
    Array& operator[] (int idx);
    friend MbArray operator+ (const MbArray& x, const MbArray& y);
    friend MbArray operator- (const MbArray& x, const MbArray& y);
    friend MbArray operator* (const MbArray& x, const MbArray& y);
    friend MbArray operator/ (const MbArray& x, const MbArray& y);

    friend MbArray operator+ (const MbArray& x, double a);
    friend MbArray operator+ (double a, const MbArray& x);
    friend MbArray operator- (const MbArray& x, double a);
    friend MbArray operator- (double a, const MbArray& x);
    friend MbArray operator* (double a, const MbArray& x);
    friend MbArray operator* (const MbArray& x, double a);
    friend MbArray operator/ (double a, const MbArray& x);
    friend MbArray operator/ (const MbArray& x, double a);

  private:
    Shapes shapes_;
    Arrays arrays_;
    void CheckSize(const MbArray& x, const MbArray& y);

    MbArray VectorOperation(const MbArray& x,
                            const MbArray& y,
                            Operation op);

    MbArray ScalarOperation(const MbArray& x, double a,
                            Operation op);

};

/*
 * Print the mbArray in the terminal.
 */
void Print(const MbArray& mb_array);

/*
 * Compute Max, Min, Sum of array
 */
double Max(const MbArray& mb_array);
double Min(const MbArray& mb_array);
double Sum(const MbArray& mb_array);
