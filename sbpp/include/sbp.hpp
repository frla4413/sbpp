//============================================================
// Name        : sbp.hpp
// Date        : June 2022
// Author      : Fredrik Laurén
// Description : hpp file for the base class Sbp
//============================================================

/*
 * Base class for Sbp operators.
 * The 2D derivative approximations uses the assumption that we
 * can use a Kronecker approach. The virtual function Diff
 * prodeuces a 1D approximation of the derivative, which is used in 
 * DXi and DEta.
 * 
 *
 * Member functions:
 *  o f_eta = DEta(const Array& f); (f_y on rectangle grid)
 *  o f_xi = DXi(const Array& f);   (f_x on rectangle grid)
 *
 * Protected variables:
 *  o Nx_, Ny_
 *  o d_xi_, d_eta_;
 *
 * Protected functions:
 *  o Diff - produces 1D approximation of derivative
 */

#pragma once
#include "array.hpp"

class Sbp {
  public:
    Sbp(){};
    Sbp(int Nx, int Ny, double d_xi, double d_eta);
    virtual ~Sbp() {};

    int Nx() const;
    int Ny() const;
    double d_xi() const;
    double d_eta() const;

    Array DEta(const Array& f) const;
    Array DXi(const Array& f) const;

//    Array DEtaT(const Array&);
//    Array DXiT(const Array&);

    virtual double Integrate(const Array&) const = 0;
    Array GetP() const;

    /*
     * Return 1D-quadrature weights.
     * Can be used when integrating along boundaries
     * and when generating the P-vector.
     */
    virtual Array GetQuadratureWeights(int size) const = 0;

  protected:
    int Nx_ = -1 , Ny_ = -1;
    double d_xi_ = -1, d_eta_ = -1;
    virtual void Diff(const double* f, double* df,
                      int,int,double) const = 0;
//    virtual void DT(const double* f,double* DTf,int,int,double) = 0;
};
