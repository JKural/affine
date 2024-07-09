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

#include "affine.h"
#include "comparisons.h"
#include "intersection.h"

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

TEST(SpanningSpaceTest, onePoint)
{
  const Point1D p0_1d{ 0.0 };
  const Point2D p0_2d{ 0.0, 0.0 };
  const Point3D p0_3d{ 0.0, 0.0, 0.0 };
  const auto space_1d = affine::Affine_space<double, 1>::spanning_space(
    { p0_1d }, affine::Equal_to_precision());
  const auto space_2d = affine::Affine_space<double, 2>::spanning_space(
    { p0_2d }, affine::Equal_to_precision());
  const auto space_3d = affine::Affine_space<double, 3>::spanning_space(
    { p0_3d }, affine::Equal_to_precision());
  EXPECT_PRED3(element_test_1d, space_1d, p0_1d, affine::Equal_to_precision());
  EXPECT_PRED3(element_test_2d, space_2d, p0_2d, affine::Equal_to_precision());
  EXPECT_PRED3(element_test_3d, space_3d, p0_3d, affine::Equal_to_precision());
  EXPECT_EQ(space_1d.dimension(), 0);
  EXPECT_EQ(space_2d.dimension(), 0);
  EXPECT_EQ(space_3d.dimension(), 0);
}

TEST(SpanningSpaceTest, twoPoints)
{
  const Point1D p0_1d{ 0.0 };
  const Point2D p0_2d{ 1.0, 0.0 };
  const Point3D p0_3d{ 3.0, -1.0, -2.0 };
  const Point1D p1_1d{ 1.0 };
  const Point2D p1_2d{ 0.0, 1.0 };
  const Point3D p1_3d{ 4.0, 8.0, 3.0 };
  const auto space_1d = affine::Affine_space<double, 1>::spanning_space(
    { p0_1d, p1_1d }, affine::Equal_to_precision());
  const auto space_2d = affine::Affine_space<double, 2>::spanning_space(
    { p0_2d, p1_2d }, affine::Equal_to_precision());
  const auto space_3d = affine::Affine_space<double, 3>::spanning_space(
    { p0_3d, p1_3d }, affine::Equal_to_precision());
  EXPECT_PRED3(element_test_1d, space_1d, p0_1d, affine::Equal_to_precision());
  EXPECT_PRED3(element_test_2d, space_2d, p0_2d, affine::Equal_to_precision());
  EXPECT_PRED3(element_test_3d, space_3d, p0_3d, affine::Equal_to_precision());
  EXPECT_PRED3(element_test_1d, space_1d, p1_1d, affine::Equal_to_precision());
  EXPECT_PRED3(element_test_2d, space_2d, p1_2d, affine::Equal_to_precision());
  EXPECT_PRED3(
    element_test_3d, space_3d, p1_3d, affine::Equal_to_precision(1e-14));
  EXPECT_EQ(space_1d.dimension(), 1);
  EXPECT_EQ(space_2d.dimension(), 1);
  EXPECT_EQ(space_3d.dimension(), 1);
}

TEST(SpanningSpaceTest, threePoints)
{
  const Point3D p0_3d{ 3.0, -1.0, -2.0 };
  const Point3D p1_3d{ 4.0, 8.0, 3.0 };
  const Point3D p2_3d{ -12.0, 6.0, 6.0 };
  const auto space_3d = affine::Affine_space<double, 3>::spanning_space(
    { p0_3d, p1_3d, p2_3d }, affine::Equal_to_precision());
  EXPECT_PRED3(element_test_3d, space_3d, p0_3d, affine::Equal_to_precision());
  EXPECT_PRED3(element_test_3d, space_3d, p1_3d, affine::Equal_to_precision());
  EXPECT_PRED3(element_test_3d, space_3d, p2_3d, affine::Equal_to_precision());
  EXPECT_EQ(space_3d.dimension(), 2);
}
