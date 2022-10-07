#include <gtest/gtest.h>
#include "sbp_tests.hpp"
#include "mb_grid_tests.hpp"
#include "mbsbp_tests.hpp"
#include "interp_tests.hpp"

/*
 * A set of unit tests that always should hold.
 * Check the .hpp files for more information of the tests.
 */
TEST(SbpTest, BasicAssertions) {

  int Nx = 31, Ny = 11;
  for (auto order : {2,4} ) {
    EXPECT_EQ(TestSbpIntegration(Nx,Ny, order), 1);
    EXPECT_EQ(TestSbpP(Nx,Ny, order), 1);
    EXPECT_EQ(TestSbpDiff(Nx,Ny, order), 1);
  }
}

TEST(MbGridTest, BasicAssertions) {
  EXPECT_EQ(TestMbGridAnnulusInterfaces(), 1);
  EXPECT_EQ(TestMbGridSquareInterfaces(),1);
}

/*
 * Test integration on annulus (consisting of four sub-grids)
 * Expect:
 *    o Sbp21: 2
 *    o Sbp42: 4
 * Test derivatives on annulus
 * Expect:
 *    o Sbp21: 1.5 in both directions
 *    o Sbp42: 2.5 in both directions
 */
TEST(MbSbpTest, BasicAssertions) {
  EXPECT_EQ(RunMbSbpTestIntegration(),1);
  EXPECT_EQ(RunMbSbpTestDifferentiaion(),1);
}

/*
 *
 * Test interpolation operators
 * interp21 is exact for 1st degree polynomial.
 * interp42 is exact for 2nd degree polynomial.
 */

TEST(InterpTest, BasicAssertions) {
  int Nc = 11, Nf = 21;
  std::vector<int> orders{2,4};
  for (auto& order : orders)
    EXPECT_EQ(TestInterpolation(11, 21, order),1);
}
