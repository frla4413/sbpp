#include <memory> // for unique_ptr
#include "array.hpp"
#include "sbp.hpp"
#include "sbp21.hpp"
#include "sbp42.hpp"
#include "mesh.hpp"

/*
 * Test Sbp-operators on a square:
 *  o TestSbpIntegration - test integration routine
 *  o TestSbpP - test Generating P (including GetQuadratureWeights)
 *  o TestSbpDiff - test DXi and DEta
 */

int TestSbpIntegration(int Nx, int Ny, int order) {

  auto block = CartesianGrid(Nx, Ny);

  double dx = 1.0/static_cast<double>(Nx-1);
  double dy = 1.0/static_cast<double>(Ny-1);

  std::unique_ptr<Sbp> sbp;

  switch (order) {
    case 4:
      sbp = std::make_unique<Sbp42>(Nx, Ny, dx, dy);
      break;
    default:
      sbp = std::make_unique<Sbp21>(Nx, Ny, dx, dy);
      break;
  }

  auto ones{Ones(Nx, Ny)};
  auto area = sbp->Integrate(ones);
  auto error = std::abs(area - 1.0);

  area = sbp->Integrate(block.x);
  error = std::max(error, std::abs(area - 0.5));

  area = sbp->Integrate(block.y);
  error = std::max(error, std::abs(area - 0.5));

  if (order == 4) {
    area = sbp->Integrate(block.x*block.x);
    error = std::max(error, std::abs(area - 1.0/3.0));

    area = sbp->Integrate(block.y*block.y);
    error = std::max(error, std::abs(area - 1.0/3.0));

    area = sbp->Integrate(block.y*block.y*block.y);
    error = std::max(error, std::abs(area - 0.25));

    area = sbp->Integrate(block.x*block.x*block.x);
    error = std::max(error, std::abs(area - 0.25));
  }

  return (error < 1e-14);
};

int TestSbpDiff(int Nx, int Ny, int order) {
  auto block = CartesianGrid(Nx, Ny);

  double dx = 1.0/static_cast<double>(Nx-1);
  double dy = 1.0/static_cast<double>(Ny-1);

  std::unique_ptr<Sbp> sbp;

  switch (order) {
    case 4:
      sbp = std::make_unique<Sbp42>(Nx, Ny, dx, dy);
      break;
    default:
      sbp = std::make_unique<Sbp21>(Nx, Ny, dx, dy);
      break;
  }


  auto ones{Ones(Nx, Ny)};
  auto zeros{Zeros(Nx, Ny)};
  auto error = Max(Abs(sbp->DXi(ones) - zeros));
  error = std::max(error, Max(Abs(sbp->DEta(ones) - zeros)));
  error = std::max(error, Max(Abs(sbp->DXi(block.x) - ones)));
  error = std::max(error, Max(Abs(sbp->DEta(block.y) - ones)));

  if (order == 4) {

    error = std::max(error,
            Max(Abs(sbp->DXi(block.x*block.x) - 2*block.x)));

    error = std::max(error,
            Max(Abs(sbp->DEta(block.y*block.y) - 2*block.y)));
  }

  return (error < 1e-14);
};

int TestSbpP(int Nx, int Ny, int order) {

  auto block = CartesianGrid(Nx, Ny);

  double dx = 1.0/static_cast<double>(Nx-1);
  double dy = 1.0/static_cast<double>(Ny-1);

  std::unique_ptr<Sbp> sbp;

  switch (order) {
    case 4:
      sbp = std::make_unique<Sbp42>(Nx, Ny, dx, dy);
      break;
    default:
      sbp = std::make_unique<Sbp21>(Nx, Ny, dx, dy);
      break;
  }

  Array P = sbp->GetP();

  double area = Sum(P);
  double error = std::abs(area - 1);

  area = Sum(P*block.x);
  error = std::max(error, std::abs(area - 0.5));

  area = Sum(P*block.y);
  error = std::max(error, std::abs(area - 0.5));

  if(order == 4) {
    area = Sum(P*block.x*block.x);
    error = std::max(error, std::abs(area - 1.0/3.0));

    area = Sum(P*block.y*block.y);
    error = std::max(error, std::abs(area - 1.0/3.0));

    area = Sum(P*block.x*block.x*block.x);
    error = std::max(error, std::abs(area - 0.25));

    area = Sum(P*block.y*block.y*block.y);
    error = std::max(error, std::abs(area - 0.25));
  }

  return (error < 1e-13);
}
