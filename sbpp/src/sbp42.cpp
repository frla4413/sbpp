#include "sbp42.hpp"
#include <stdexcept>
#include <string.h>

Sbp42::Sbp42(int Nx, int Ny, double dXi, double dEta) :
  Sbp(Nx,Ny,dXi,dEta) {

  // Minimum size is 8 for sbp42.
  if (Nx < 8 || Ny < 8) {
    std::string msg = "@Sbp42: Nx,Ny > 8 required!";
    throw std::invalid_argument(msg);
  }
}

void Sbp42::Diff(const double* f, double* df,
                 int stride, int length, double h) const {

  double q1 = 1/12.0;
  double q2 = 2.0/3.0;

  double q3 = 59.0/96.0;
  double q4 = 1/12.0;
  double q5 = 1/32.0;

  df[0] = (-0.5*f[0] + q3*f[stride] - q4*f[2*stride] -
           q5*f[3*stride])/h1_/h;

  df[stride] = (-q3*f[0] + q3*f[2*stride])/h2_/h;

  df[2*stride] = (q4*f[0] - q3*f[stride] + q3*f[3*stride] -
                  q4*f[4*stride])/(h3_*h);

  df[3*stride] = (q5*f[0] - q3*f[2*stride] + q2*f[4*stride] -
                  q1*f[5*stride])/h4_/h;

  df[(length-4)*stride] = (q3*f[(length-3)*stride] -
                           q5*f[(length-1)*stride] -
                           q2*f[(length-5)*stride] +
                           q1*f[(length-6)*stride])/h4_/h;

 df[(length-3)*stride] = (-q3*f[(length-4)*stride] +
                           q3*f[(length-2)*stride] -
                           q4*f[(length-1)*stride] +
                           q1*f[(length-5)*stride])/h3_/h;

  df[(length-2)*stride] = (- q3*f[(length-3)*stride] +
                             q3*f[(length-1)*stride])/h2_/h;

  df[(length-1)*stride] = (q5*f[(length-4)*stride] +
                           q4*f[(length-3)*stride] -
                           q3*f[(length-2)*stride] +
                           0.5*f[(length-1)*stride])/h1_/h;

  for(int i = 4; i < length-4; ++i)
    df[i*stride] = (q1*f[(i-2)*stride] - q2*f[(i-1)*stride] +
                    q2*f[(i+1)*stride] - q1*f[(i+2)*stride])/h;
}

/*
 * Use the SBP property 
 *    Q + QT = B  <--> QT = B  - Q <--> QT P^-1 = B P^-1 - Q P^-1 (1)
 * to get DT*f. 
 * Since  P is diagonal, we have that
 *    D = P^-1 Q <--> DT = QT P^-1 (2)
 *
 * Combining (1) and (2) yeilds
 *
 *    DT f = B P^-1 f - Q P^-1 f .
 *
 * We already have the routine for Q*v (Diff but remove P^-1)
 *
 * The routine is verified by comparing to results from sbpy
 */

//void Sbp42::DT(const double* f, double* dtf, int stride, int length, double h)
//{
//    double q1 = 1/12.0; 
//    double q2 = 2.0/3.0;
//
//    double q3 = 59.0/96.0;
//    double q4 = 1/12.0;
//    double q5 = 1/32.0;
//
//    Array2D pinv = 1.0/GetQuadratureWeights(length)/h;
//
//    // Compute -Q (pinv f)
//    dtf[0]        = -(-0.5*f[0]*pinv[0]  
//                      + q3*f[stride]*pinv[1] 
//                      - q4*f[2*stride]*pinv[2] 
//                      - q5*f[3*stride]*pinv[3]);
//
//    dtf[stride]   = -(-q3*f[0]*pinv[0] 
//                      + q3*f[2*stride]*pinv[2]);
//
//    dtf[2*stride] = -(q4*f[0]*pinv[0]
//                      - q3*f[stride]*pinv[1] 
//                      + q3*f[3*stride]*pinv[3] 
//                      - q4*f[4*stride]*pinv[4]);
//
//    dtf[3*stride] = -(q5*f[0]*pinv[0] 
//                      -q3*f[2*stride]*pinv[2]
//                      + q2*f[4*stride]*pinv[4] 
//                      - q1*f[5*stride]*pinv[5]);
//
//  dtf[(length-4)*stride] = -(q3*f[(length-3)*stride]*pinv[length-3] 
//                             - q5*f[(length-1)*stride]*pinv[length-1]
//                             - q2*f[(length-5)*stride]*pinv[length-5] 
//                             + q1*f[(length-6)*stride]*pinv[length-6]);
//
//  dtf[(length-3)*stride] = -(-q3*f[(length-4)*stride]*pinv[length-4] 
//                             + q3*f[(length-2)*stride]*pinv[length-2]
//                             - q4*f[(length-1)*stride]*pinv[length-1] 
//                             + q1*f[(length-5)*stride]*pinv[length-5]);
//
//  dtf[(length-2)*stride] = -(- q3*f[(length-3)*stride]*pinv[length-3] 
//                              + q3*f[(length-1)*stride]*pinv[length-1]);
//
//  dtf[(length-1)*stride] = -(q5*f[(length-4)*stride]*pinv[length-4]  
//                             + q4*f[(length-3)*stride]*pinv[length-3]
//                             - q3*f[(length-2)*stride]*pinv[length-2] 
//                             + 0.5*f[(length-1)*stride]*pinv[length-1]);
//
//  for(int i = 4; i < length-4; i++)
//   {
//       dtf[i*stride] = -(q1*f[(i-2)*stride]*pinv[i-2] - q2*f[(i-1)*stride]*pinv[i-1]
//                         + q2*f[(i+1)*stride]*pinv[i+1] - q1*f[(i+2)*stride]*pinv[i+2]);
//   }
//
//   // Add on boundary contribution: DT f = B (pinv f) - Q (pinv f)
//   dtf[0] -= pinv[0]*f[0];
//   dtf[(length-1)*stride] += pinv[length-1]*f[(length-1)*stride];
//}


Array Sbp42::GetQuadratureWeights(int N) const {
  std::valarray <double> weights(1,N);
  weights[0]  = h1_;
  weights[1]  = h2_;
  weights[2]  = h3_;
  weights[3]  = h4_;

  weights[N-1] = h1_;
  weights[N-2] = h2_;
  weights[N-3] = h3_;
  weights[N-4] = h4_;

  return {1, N, weights};
}

double Sbp42::Integrate(const Array& f) const {

  double sum_tot = 0;

  // integrate along columns
  for(int i = 0; i < Nx_; ++i) {
    int start = i*Ny_; //start at column i
    double sum  = 0.0;

    sum += h1_*f[start];
    sum += h1_*f[start+Ny_-1];
    sum += h2_*f[start+1];
    sum += h2_*f[start+Ny_-2];
    sum += h3_*f[start+2];
    sum += h3_*f[start+Ny_-3];
    sum += h4_*f[start+3];
    sum += h4_*f[start+Ny_-4];

    for (int j = 4; j < Ny_-4; ++j)
      sum += f[start+j];
    sum *= d_xi_*d_eta_;

    if (i == 0 || i == Nx_-1)
      sum_tot += sum*h1_;
    else if (i == 1 || i == Nx_-2)
      sum_tot += sum*h2_;
    else if (i == 2 || i == Nx_-3)
      sum_tot += sum*h3_;
    else if (i == 3 || i == Nx_-4)
      sum_tot += sum*h4_;
    else
      sum_tot += sum;
    }
    return sum_tot;
}
