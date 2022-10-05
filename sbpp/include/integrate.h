#pragma once

#include <iostream> 
#include <vector>
#include "basics.h"
#include "array2D.h"
#include "gridfunction.h"
#include "data_types.h"

namespace Integrate
{
    // -------------------------------------------------------------------------
    // ---------------------------- Explicit time integration ------------------
    // -------------------------------------------------------------------------
    DataTypes::Solution ExplicitIntegration(std::vector< double >& tspan, 
                            const Gridfunction& y0, 
                            double dt,
                     std::function<Gridfunction(double t, const Gridfunction&)> odefun);

    DataTypes::Solution ExplicitIntegrationSaveSolution(std::vector< double >& tspan, 
                                                        const Gridfunction& y0, 
                                                        double dt,
                     std::function<Gridfunction(double t, const Gridfunction&)> odefun,
            std::function<void(const Gridfunction&, const std::string&)> export_to_tec,
            std::string name_base);

    Gridfunction RK4Step(double t, const Gridfunction& y, double dt, 
                         std::function<Gridfunction(double, const Gridfunction&)> odefun);

    // -------------------------------------------------------------------------
    // ---------------------------- Implicit time integration ------------------
    // -------------------------------------------------------------------------
    DataTypes::Solution ImplicitTimeIntegraction(std::vector< double >& tspan,
                                         const Gridfunction& y0, double dt,
            std::function<Gridfunction(double, const Gridfunction&)> odefun,
            std::function<Gridfunction(double, const Gridfunction&)> Jv);

    DataTypes::Solution ImplicitTimeIntegractionSaveSolution(std::vector< double >& tspan,
                                         const Gridfunction& y0, double dt,
            std::function<Gridfunction(double, const Gridfunction&)> odefun,
            std::function<Gridfunction(double, const Gridfunction&)> Jv,
            std::function<void(const Gridfunction&, const std::string&)> export_to_tec,
            std::string name_base);


    Gridfunction BDF1Step(const Gridfunction& y_prev, 
                          double dt, double t,
            std::function<Gridfunction(double, const Gridfunction&)> odefun,
            std::function<Gridfunction(double, const Gridfunction&)> Jv);

    Gridfunction BDF2Step(const Gridfunction& y_prev_prev, 
                          const Gridfunction& y_prev,
                          double dt, double t,
            std::function<Gridfunction(double, const Gridfunction&)> odefun,
            std::function<Gridfunction(double, const Gridfunction&)> Jv);

    Gridfunction BDF3Step(const Gridfunction& y_prev_prev_prev, 
                          const Gridfunction& y_prev_prev,
                          const Gridfunction& y_prev, 
                          double dt, double t,
            std::function<Gridfunction(double, const Gridfunction&)> odefun,
            std::function<Gridfunction(double, const Gridfunction&)> Jv);

    // documentation: 
    // https://software.intel.com/content/www/us/en/develop/documentation/onemkl-developer-reference-c/top.html
    //
    //The algorithm for MKL's GMRES is inspired by 
    // http://sep.stanford.edu/sep/claudio/Research/Prst_ExpRefl/ShtPSPI/intel/mkl/10.0.3.020/examples/solver/source/dcsrilu0_exampl1.c
    Gridfunction GMRESMKL(std::function<Gridfunction
                       (const Gridfunction&)> Ax,
                        const Gridfunction& b, const Gridfunction& x0); 
}
