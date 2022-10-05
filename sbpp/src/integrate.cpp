#include "integrate.h"
#include <math.h>
#include "mkl.h"
#include "tqdm.hpp"

namespace Integrate
{

    // -------------------------------------------------------------------------
    // ---------------------------- Explicit time integration ------------------
    // -------------------------------------------------------------------------

    /*
     * Explicit integration. 
     * Integrates y' = odefun(t,y) from tspan[0] to tspan[1].
     * Input: 
     *   o tspan[t0, t1] - limits of the integration
     *   o y0 - initial condition
     *   o dt - step size
     *   o odefun - std::function. Gridfunction y' = odefun(double t, Gridfunction y)
     *
     * Output:
     *   o DataTypes::Solution containing the approximation to y(t1) and t1.
     */

    DataTypes::Solution ExplicitIntegration(std::vector< double >& tspan,
                                            const Gridfunction& y0, double dt,
            std::function<Gridfunction(double, const Gridfunction&)> odefun)
    {
         int NoS  = (int) ceil((tspan[1] - tspan[0])/dt) + 1;

         Gridfunction y = y0;
         double t = tspan[0];

         for(int i : tq::trange(NoS))
         {
            y = RK4Step(t, y, dt, odefun);
            t += dt;
         }
         
         struct DataTypes::Solution solution;
         solution.sol = y;
         solution.time = t;
         return solution;
    }

    /* Explicit integration, save every XX solution step 
     * so that 15 snapshots are saved per time unit.
     *
     * Input and Output are similar to ExplicitIntegration except for
     *
     * Input: 
     *   o export to tec: std::function to write y to .tec format.
     *                    void export_to_tec(Gridfunction y, std::string name)
     *   o std::string name_base - name_base of the files to be saved. 
     *                             Each saved file will get an index extension.
     *                             Ex: name_base = "save/sol"
     *                             Then each file will be saved as 
     *                             "save/sol1, "save/sol2", ...
     */
   DataTypes::Solution ExplicitIntegrationSaveSolution(std::vector< double >& tspan,
                            const Gridfunction& y0, double dt,
            std::function<Gridfunction(double, const Gridfunction&)> odefun,
            std::function<void(const Gridfunction&, const std::string&)> export_to_tec,
            std::string name_base)
    {
         int NoS  = (int) ceil((tspan[1] - tspan[0])/dt) + 1;

         Gridfunction y = y0;
         double t = tspan[0];
         int steps_per_time = static_cast<int>( ceil(NoS/(tspan[1] - tspan[0])));
         int interval = ceil(steps_per_time/15);
         int tec_pos = 0;

         for(int i : tq::trange(NoS))
         {
            if(i % interval == 0)
            {
               std::string file_name = name_base + std::to_string(tec_pos);
               export_to_tec(y,file_name);
               tec_pos ++;
            }

            y = RK4Step(t, y, dt, odefun);
            t += dt;
         }
         
         struct DataTypes::Solution solution;
         solution.sol = y;
         solution.time = t;
         return solution;
    }

   /* Runge-Kutta 4th order step. 
    * To be used in explicit time integration.
    * Return y_{k+1} given y_k and t_k.
    *
    * Input: 
    *    o t - current time
    *    o y - current value of y_k
    *    o dt - time step
    *    o odefun - same as in ExplicitIntegration.
    *
    * Output:
    *    o y_{k+1}
    */
    Gridfunction RK4Step(double t, const Gridfunction& y, double dt, 
                         std::function<Gridfunction(double, const Gridfunction&)> odefun)
    {

       Gridfunction k1 = odefun(t,y);
       Gridfunction k2 = odefun(t + 0.5*dt, y + 0.5*dt*k1);
       Gridfunction k3 = odefun(t + 0.5*dt, y + 0.5*dt*k2);
       Gridfunction k4 = odefun(t + dt, y + dt*k3);

       return y + (dt/6.0)*(k1 + 2*k2 + 2*k3 + k4);
    }

    // -------------------------------------------------------------------------
    // ---------------------------- Implicit time integration ------------------
    // -------------------------------------------------------------------------
    //
    /*
     * Implicit integration. 
     * Integrates y' + odefun(t,y) = 0 from tspan[0] to tspan[1].
     * Input: 
     *   o tspan[t0, t1] - limits of the integration
     *   o y0 - initial condition
     *   o dt - step size
     *   o odefun - std::function. Gridfunction y' = odefun(double t, Gridfunction y)
     *   o Jv - std::function. The action of the Jacobian to odefun on a Gridfunction v.
     *
     * Output:
     *   o DataTypes::Solution containing the approximation to y(t1) and t1.
     */

    DataTypes::Solution ImplicitTimeIntegraction(std::vector< double >& tspan,
                                         const Gridfunction& y0, double dt,
            std::function<Gridfunction(double, const Gridfunction&)> odefun,
            std::function<Gridfunction(double, const Gridfunction&)> Jv)
    {
         int NoS  = static_cast<int>(ceil((tspan[1] - tspan[0])/dt));

         Gridfunction y = y0, y_prev = y0, y_prev_prev = y0;
         double t = tspan[0];

         for(int i : tq::trange(NoS))
         {

            t = tspan[0] + (i+1)*dt;
            if(i == 0)
            {
               y = BDF1Step(y_prev, dt, t, odefun, Jv);
            }
            else
            {
               y_prev_prev = y_prev;
               y_prev = y;
               y = BDF2Step(y_prev_prev,y_prev, dt, t, odefun, Jv);
            }

         }
         DataTypes::Solution solution;
         solution.sol = y;
         solution.time = t;
         return solution;
    }

    DataTypes::Solution ImplicitTimeIntegractionSaveSolution(std::vector< double >& tspan,
                                         const Gridfunction& y0, double dt,
            std::function<Gridfunction(double, const Gridfunction&)> odefun,
            std::function<Gridfunction(double, const Gridfunction&)> Jv,
            std::function<void(const Gridfunction&, const std::string&)> export_to_tec,
            std::string name_base)
    {
         int NoS  = static_cast<int>(ceil((tspan[1] - tspan[0])/dt));
         int steps_per_time = static_cast<int>( ceil(NoS/(tspan[1] - tspan[0])));
         int interval = ceil(steps_per_time/15);

         Gridfunction y = y0, y_prev = y0, y_prev_prev = y0, y_prev_prev_prev = y0;
         double t = tspan[0];

         int tec_pos = 0;
         for(int i : tq::trange(NoS))
         {
            if(i % interval == 0)
            {
               std::string file_name = name_base + std::to_string(tec_pos);
               export_to_tec(y,file_name);
               tec_pos ++;
            }

            t = tspan[0] + (i+1)*dt;
            if(i == 0)
            {
               y = BDF1Step(y_prev, dt, t, odefun, Jv);
            }
            else if(i == 1)
            {
               y_prev_prev = y_prev;
               y_prev = y;
               y = BDF2Step(y_prev_prev,y_prev, dt, t, odefun, Jv);
            }
            else
            {
               y_prev_prev_prev = y_prev_prev;
               y_prev_prev = y_prev;
               y_prev = y;
               y = BDF3Step(y_prev_prev_prev, y_prev_prev, y_prev, dt, t, odefun, Jv);
            }
         }
         DataTypes::Solution solution;
         solution.sol = y;
         solution.time = t;
         return solution;
    }

    Gridfunction BDF1Step(const Gridfunction& y_prev, double dt, double t,
            std::function<Gridfunction(double, const Gridfunction&)> odefun,
            std::function<Gridfunction(double, const Gridfunction&)> Jv)
    {
         Gridfunction y_new = y_prev; 

         auto F = [&dt, y_prev, odefun](double t, const Gridfunction& y_new)
         {
            return (y_new - y_prev)  + dt*odefun(t, y_new);
         };

         auto Jv_func = [&dt, Jv, t](const Gridfunction& v)
         {
            return v + dt*Jv(t, v);
         };

         Gridfunction b, delta, init;
         double l2err;

         b = F(t,y_new);
         for(int k = 0; k < 10; k++)
         {
            delta = GMRESMKL(Jv_func, b, init);

            y_new -= delta;
            b = F(t,y_new);
            l2err = Basics::L2Norm(b);
            if(l2err < 1e-10){ break; }
         }
         if (l2err > 1e-10) {std::cout << "BDF1: GMRES NOT CONVERGED " << 
                                          "l2err: " << l2err << std::endl;}

         return y_new;
    }

    Gridfunction BDF2Step(const Gridfunction& y_prev_prev, const Gridfunction& y_prev, 
            double dt, double t,
            std::function<Gridfunction(double, const Gridfunction&)> odefun,
            std::function<Gridfunction(double, const Gridfunction&)> Jv)
    {
         Gridfunction y_new = y_prev; 

         auto F = [&dt, y_prev_prev, y_prev, odefun](double t, 
                                 const Gridfunction& y_new)
         {
            return (y_new - (4.0/3.0)*y_prev + (1.0/3.0)*y_prev_prev)  
                   + dt*(2.0/3.0)*odefun(t, y_new);
         };

         auto Jv_func = [&dt, Jv, t](const Gridfunction& v)
         {
            return v + dt*(2.0/3.0)*Jv(t, v);
         };

         Gridfunction b, delta, init;
         double l2err;

         b = F(t,y_new);
         for(int k = 0; k < 10; k++)
         {
            delta = GMRESMKL(Jv_func, b, init);

            y_new -= delta;
            b = F(t,y_new);
            l2err = Basics::L2Norm(b);
            if(l2err < 1e-10){ break; }
         }
         if (l2err > 1e-10) {std::cout << "BDF2: GMRES NOT CONVERGED " << 
                                          "l2err: " << l2err << std::endl;}

         return y_new;
    }

    Gridfunction BDF3Step(const Gridfunction& y_prev_prev_prev, 
                          const Gridfunction& y_prev_prev, 
                          const Gridfunction& y_prev, 
                          double dt, double t,
            std::function<Gridfunction(double, const Gridfunction&)> odefun,
            std::function<Gridfunction(double, const Gridfunction&)> Jv)
    {
         Gridfunction y_new = y_prev;

         auto F = [&dt, y_prev_prev_prev,y_prev_prev, y_prev, odefun](double t, 
                                 const Gridfunction& y_new)
         {
            return y_new - (18.0/11.0)*y_prev + (9.0/11.0)*y_prev_prev  
                         - (2.0/11.0)*y_prev_prev_prev
                   + dt*(6.0/11.0)*odefun(t, y_new);
         };

         auto Jv_func = [&dt, Jv, t](const Gridfunction& v)
         {
            return v + dt*(6.0/11.0)*Jv(t, v);
         };

         Gridfunction b, delta, init;
         double l2err;

         b = F(t,y_new);
         for(int k = 0; k < 10; k++)
         {
            delta = GMRESMKL(Jv_func, b, init);

            y_new -= delta;
            b = F(t,y_new);
            l2err = Basics::L2Norm(b);
            if(l2err < 1e-10){ break; }
         }
         if (l2err > 1e-10) {std::cout << "BDF2: GMRES NOT CONVERGED" << 
                                          "l2err: " << l2err << std::endl;}
         return y_new;
    }


    /*
     * Wrapper of Intel MKL dfgmres. Solves the system Ax = b. 
     * Input: 
     *   Ax - A function returning the action of A on a Gridfunction
     *   b - right-hand side
     *   x0 - initial guess
     */
    Gridfunction GMRESMKL(std::function<Gridfunction(const Gridfunction&)> Ax,
                          const Gridfunction& b, const Gridfunction& x0) 
    {

      int N = b.GetTotalSize();
      int iter;
      int max_iter = 1000, subspace_size = 200;
      double rel_tol = 1e-6, abs_tol = 1e-6;

      std::valarray<double> sol(N), b_arr = b.GridfunctionToValarray();
      Gridfunction work;
      int RCI_request;
      std::valarray<int> ipar(128);
      std::vector<double> dpar(128);

      int tmp_size = (2*subspace_size + 1)*N +subspace_size*(subspace_size + 9)/2 + 1;
      std::valarray<double> tmp(tmp_size);

      dfgmres_init(&N,&sol[0],&b_arr[0],&RCI_request,&ipar[0],&dpar[0],&tmp[0]);

      ipar[0]  = N;
      ipar[1]  = 6; // error messages
//      ipar[2]  = 1;
//      ipar[3]  = 0; // current iteration number
      ipar[4]  = max_iter; // maximum iterations
//      ipar[5]  = 1; // error messages
//      ipar[6]  = 1; // error messages
//      ipar[7]  = 1; // 1--> check ipar[3] < ipar[4], defalut: 1
      ipar[8]  = 0;
      ipar[9]  = 0;
      ipar[10] = 0; // 0 --> not use preconditioner
      ipar[11] = 0;
      ipar[12] = 0;
      ipar[13] = 0;
      ipar[14] = subspace_size;

      dpar[0] = rel_tol; //relative tolerance
      dpar[1] = abs_tol; //absolute tolerance

      dfgmres_check(&N,&sol[0],&b_arr[0],&RCI_request,&ipar[0],&dpar[0],&tmp[0]);

GMRES: 
      {
         dfgmres(&N, &sol[0], &b_arr[0], &RCI_request, &ipar[0], &dpar[0], &tmp[0]);

         switch (RCI_request)
         {
            case 0:
               goto COMPLETE;
            case 1:
               {
                  work = Gridfunction(tmp[std::slice(ipar[21]-1,N,1)], 
                                      b.GetShapes());
                  work = Ax(work);
                  tmp[std::slice(ipar[22]-1,N,1)] = work.GridfunctionToValarray();
                  goto GMRES;
               }
            case 4:
               {
                  if(dpar[6] < 1e-12)
                     goto COMPLETE;
                  else
                     goto GMRES;
               }
         }
      }

COMPLETE:
      {
      ipar[12] = 0;
      dfgmres_get(&N, &sol[0], &b_arr[0], &RCI_request, 
                  &ipar[0], &dpar[0], &tmp[0], &iter);
      return Gridfunction(sol, b.GetShapes());
      }
    }
}
