#include "basics.hpp"

Array GetConvergenceRates(const std::valarray<int>& n_array,
                       const std::valarray<double>& err_array) {
  std::valarray<double> rates(err_array.size() -1);
  for (int i = 0; i < err_array.size()-1; i++) {
    rates[i] = log(err_array[i]/err_array[i+1])/
               log(static_cast<double>(n_array[i+1]-1)/
                   static_cast<double>(n_array[i]-1));
  }
  return Array(1,rates.size(),rates);
}

void PrintConvergence(const std::valarray<int>& n_array,
                     const std::valarray<double>& err_array) {

  std::cout<< "N: \t Error: \t Rates\n";
  Array rates = GetConvergenceRates(n_array, err_array);
  std::cout<< n_array[0] << "\t" <<
  std::scientific << std::setprecision(3)
                  << err_array[0] << "\n";
  for(int i = 1; i < n_array.size(); i++)
    std::cout << n_array[i] << "\t" << std::scientific
              << std::setprecision(3) << err_array[i] << "\t"
              << std::fixed << std::setprecision(3)
              << rates[i-1] << "\n";
}
