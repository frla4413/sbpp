#pragma once
#include <iostream>
#include <fstream> // write to file
#include <cmath>
#include <fstream>
#include <string>
#include <valarray>
#include <iomanip> // for setprecision in convergence table
#include "mbgrid.hpp"


void PrintConvergence(const std::valarray<int>& n_vec,
                      const std::valarray<double>& err_vec);

Array GetConvergenceRates(const std::valarray<int>& n_vec,
                          const std::valarray<double>& err_vec);
