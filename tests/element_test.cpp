#include <array>
#include <cmath>
#include <functional>
#include <list>
#include <stdexcept>

#include <gtest/gtest.h>
// skips -Woverloaded-virtual warnings for capd library
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Woverloaded-virtual"
#include <capd/capdlib.h>
#pragma GCC diagnostic pop

#include "affine/affine.h"

namespace {

using Point1D = capd::vectalg::Vector<double, 1>;
using Point2D = capd::vectalg::Vector<double, 2>;
using Point3D = capd::vectalg::Vector<double, 3>;

const auto element_test_1d = std::mem_fn(
  &affine::Affine_space<double, 1>::element<affine::Equal_to_precision>);
const auto element_test_2d = std::mem_fn(
  &affine::Affine_space<double, 2>::element<affine::Equal_to_precision>);
const auto element_test_3d = std::mem_fn(
  &affine::Affine_space<double, 3>::element<affine::Equal_to_precision>);

} // namespace

TEST(ElementTest, point)
{
  const Point1D p0_1d{ 0.0 };
  const Point1D p1_1d{ 1.0 };
  const Point3D p0_3d{ 0.15937, -7.87418, -0.0342971 };
  const Point3D p1_3d{ -0.364325, 4.09101, 9.47693 };
  const affine::Affine_space space_1d = p0_1d;
  const affine::Affine_space space_3d = p0_3d;

  EXPECT_PRED3(element_test_1d, space_1d, p0_1d, affine::Equal_to_precision());
  EXPECT_PRED3(std::not_fn(element_test_1d),
               space_1d,
               p1_1d,
               affine::Equal_to_precision());
  EXPECT_PRED3(element_test_3d, space_3d, p0_3d, affine::Equal_to_precision());
  EXPECT_PRED3(std::not_fn(element_test_3d),
               space_3d,
               p1_3d,
               affine::Equal_to_precision());
}

TEST(ElementTest, line)
{
  const Point1D p0_1d{ 3.7452 };
  const Point2D p0_2d{ -8.49156, 1.51593 };
  const Point3D p0_3d{ -8.62318, 0.216789, -6.73864 };
  const Point1D p1_1d{ 0.15937 };
  const Point2D p1_2d{ -7.87418, -0.0342971 };
  const Point3D p1_3d{ -0.364325, 4.09101, 9.47693 };
  const Point2D p2_2d{ -8.33719, -0.341978 };
  const Point3D p2_3d{ -5.91821, 7.95807, -9.30484 };
  const double t_1 = 6.02605;
  const double t_2 = -4.34583;
  const double t_3 = 6.31083;
  const double s_2 = -0.613258;
  const double s_3 = 9.56608;
  const affine::Affine_space space_1d(
    p0_1d, std::vector{ p1_1d }, affine::Equal_to_precision());
  const affine::Affine_space space_2d(
    p0_2d, std::vector{ p1_2d }, affine::Equal_to_precision());
  const affine::Affine_space space_3d(
    p0_3d, std::vector{ p1_3d }, affine::Equal_to_precision());

  EXPECT_PRED3(element_test_1d,
               space_1d,
               p0_1d + t_1 * p1_1d,
               affine::Equal_to_precision());
  EXPECT_PRED3(element_test_2d,
               space_2d,
               p0_2d + t_2 * p1_2d,
               affine::Equal_to_precision());
  EXPECT_PRED3(std::not_fn(element_test_2d),
               space_2d,
               p0_2d + s_2 * p2_2d,
               affine::Equal_to_precision());
  EXPECT_PRED3(element_test_3d,
               space_3d,
               p0_3d + t_3 * p1_3d,
               affine::Equal_to_precision(1e-14));
  EXPECT_PRED3(std::not_fn(element_test_3d),
               space_3d,
               p0_3d + s_3 * p2_3d,
               affine::Equal_to_precision());
}

TEST(ElementTest, plane)
{
  const Point3D p0{ -6.74407, 9.95479, -8.01248 };
  const Point3D p1{ 9.03308, 2.52597, 9.46997 };
  const Point3D p2{ -0.829302, -0.334162, 4.82537 };
  const Point3D p3{ 2.75133, 2.39972, -7.47072 };
  const double t_0 = -6.00478;
  const double s_0 = 9.07631;
  const double t_1 = -4.82656;
  const double s_1 = -6.92954;
  const double r = 2.03547;
  const affine::Affine_space space(
    p0, std::vector{ p1, p2 }, affine::Equal_to_precision());
  EXPECT_PRED3(element_test_3d,
               space,
               p0 + t_0 * p1 + s_0 * p2,
               affine::Equal_to_precision(1e-14));
  EXPECT_PRED3(element_test_3d,
               space,
               p0 + t_1 * p1 + s_1 * p2,
               affine::Equal_to_precision(1e-14));
  EXPECT_PRED3(std::not_fn(element_test_3d),
               space,
               p0 + r * p3,
               affine::Equal_to_precision());
}