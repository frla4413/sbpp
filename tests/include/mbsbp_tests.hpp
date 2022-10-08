#include "mbsbp.hpp"
#include "mesh.hpp"
#include "basics.hpp"

/* * Integrate f(x,y) = exp(-x^2 - y^2) * on the annulus with r0/r1 as its inner/outer radius.
 * The analytic value is F = pi*(exp(-r0^2) - exp(-r1^2)).
 */
double MbSbpTestIntegration(int N, int order) {

  double r0 = 0.1, r1 = 1;
  std::vector<Block> blocks{Annulus(N, r0, r1)};
  MbSbp sbp{blocks, order};

  auto F = [] (double x, double y){ return exp(-(x*x + y*y)); };
  auto F_vec{sbp.Evaluate(F)};

  double analytic_integral = M_PI*(exp(-r0*r0) - exp(-r1*r1));

  double error = std::abs(sbp.Integrate(F_vec) -
                          analytic_integral);
  return error;
}

/*
 * Print convergence tables for integration test on an annulus.
 * The final rate should be close to 4.
 */
bool RunMbSbpTestIntegration() {

  std::cout << "Test Integration. Expected rates 2 & 4."
            << std::endl;

  std::valarray<int> N_arr{21,41,81,161};
  std::valarray<double> err_arr(N_arr.size());

  Array rates;
  for(auto& order : {2,4}) {
    for(int i = 0; i < N_arr.size(); ++i)
       err_arr[i] = MbSbpTestIntegration(N_arr[i], order);
    rates = GetConvergenceRates(N_arr, err_arr);
//    PrintConvergence(N_arr, err_arr);
    err_arr.resize(N_arr.size());
  }
  double final_rate = std::abs(rates[rates.size()-1] - 4);
  return  final_rate < 0.2;;
}

double F(double x, double y){
  return sin(2*M_PI*x)*sin(2*M_PI*y);
}

double dFdx(double x, double y) {
  return 2*M_PI*cos(2*M_PI*x)*sin(2*M_PI*y);
}

double dFdy(double x, double y) {
  return 2*M_PI*sin(2*M_PI*x)*cos(2*M_PI*y);
}

double MbSbpTestDifferentiation(int N, int order,
                                std::string direction) {

  std::vector<Block> blocks = Annulus(N, 0.1, 1);

  blocks = {blocks[0]};
  MbSbp sbp{blocks, order};
  MbArray analytic_dF, computed_dF;

  auto F_vec{sbp.Evaluate(F)};

  if(direction == "X") {
    analytic_dF = sbp.Evaluate(dFdx);
    computed_dF = sbp.Dx(F_vec);
  }
  else {
    analytic_dF = sbp.Evaluate(dFdy);
    computed_dF = sbp.Dy(F_vec);
  }

  std::vector<Array> err2_vec;

  for(int i = 0; i < blocks.size(); ++i) {
    err2_vec.push_back(analytic_dF[i] - computed_dF[i]);
    err2_vec[i] = err2_vec[i]*err2_vec[i];
  }
  double L2_error = sqrt(sbp.Integrate(err2_vec));
  return L2_error;
}

bool RunMbSbpTestDifferentiaion() {

  std::cout << "Test Differentiation. Expected rates 1.5 & 2.5." 
            << std::endl;

  std::valarray<int> N_arr = {41,81,161};
  std::valarray<double> err_arr(N_arr.size());
  std::vector<double> theoretical_rates{1.5, 2.5};

  double final_rate = -1;
  int idx = 0;
  for(auto& order : {2,4}) {
    std::cout << "Test x-derivative:" << std::endl;
    for(int i = 0; i < N_arr.size(); i++)
      err_arr[i] = MbSbpTestDifferentiation(N_arr[i], order,"X");
    Array rates = GetConvergenceRates(N_arr, err_arr);
    final_rate = std::max(final_rate,
                 std::abs(rates[rates.size()-1] -
                          theoretical_rates[idx]));
 //   PrintConvergence(N_arr, err_arr);

    err_arr.resize(N_arr.size());
    std::cout << "Test y-derivative:" << std::endl;
    for(int i = 0; i < N_arr.size(); i++)
      err_arr[i] = MbSbpTestDifferentiation(N_arr[i], order,"Y");
    rates = GetConvergenceRates(N_arr, err_arr);
  //  PrintConvergence(N_arr, err_arr);
    err_arr.resize(N_arr.size());
    final_rate = std::max(final_rate,
                 std::abs(rates[rates.size()-1] -
                          theoretical_rates[idx]));
    ++idx;
  }
  return  final_rate < 0.2;;
}

bool MbSbpTestSbpProperty(int N, int order) {
  auto UFunc = [](double x, double y) {
    return exp(-x*y);
  };
  auto VFunc = [](double x, double y) { 
    return cos(2*M_PI*x)*sin(2*M_PI*y);
  };

  std::vector<Block> blocks = Annulus(N, 0.1, 1.0);

  MbSbp sbp = MbSbp(blocks, order);

  MbArray u = sbp.Evaluate(UFunc);
  MbArray v = sbp.Evaluate(VFunc);

  for(int i = 0; i < blocks.size(); ++i) {
    u += i;
    v += i;
  }

  MbArray dudx = sbp.DxAndInterface(u);
  MbArray dvdx = sbp.DxAndInterface(v);

  MbArray dudy = sbp.DyAndInterface(u);
  MbArray dvdy = sbp.DyAndInterface(v);

  double sbp_error = std::abs(sbp.Integrate(v*dudx) +
                              sbp.Integrate(dvdx*u) +
                              sbp.Integrate(v*dudy) +
                              sbp.Integrate(dvdy*u));

  auto boundaries = sbp.boundaries();

  Array u_bd, v_bd, bd_quad;
  ArrayPair normals;
  for(auto& boundary : boundaries) {
    u_bd = sbp.ToBlockBoundary(u, boundary.block_idx,
                               boundary.side);
    v_bd = sbp.ToBlockBoundary(v, boundary.block_idx,
                               boundary.side);
    bd_quad = sbp.GetBoundaryQuadrature(boundary.block_idx,
                                        boundary.side);
    normals = sbp.GetNormals(boundary.block_idx, boundary.side);
    sbp_error -= Sum(bd_quad*v_bd*
                    (normals.a1 + normals.a2)*u_bd);
  }
  return std::abs(sbp_error) < 1e-12;
}
