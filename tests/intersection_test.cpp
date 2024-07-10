#include <gtest/gtest.h>

#include <array>
#include <cmath>
#include <functional>
#include <list>
#include <stdexcept>

// skips -Woverloaded-virtual warnings for capd library
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Woverloaded-virtual"
#include <capd/capdlib.h>
#pragma GCC diagnostic pop

#include "affine/affine.h"
#include "./test_utils.h"


TEST(IntersectionTest, pointWithPoint)
{
  const Point1D p0_1d{ 0.0 };
  const Point1D p1_1d{ 1.0 };
  const Point2D p0_2d{ 0.0, 0.0 };
  const Point2D p1_2d{ 1.0, 0.0 };
  const affine::Affine_space space0_1d(p0_1d);
  const affine::Affine_space space1_1d(p1_1d);
  const affine::Affine_space space0_2d(p0_2d);
  const affine::Affine_space space1_2d(p1_2d);

  const auto intersection0_1d =
    intersection(space0_1d, space0_1d, affine::Equal_to_precision());
  const auto intersection1_1d =
    intersection(space0_1d, space1_1d, affine::Equal_to_precision());
  const auto intersection0_2d =
    intersection(space0_2d, space0_2d, affine::Equal_to_precision());
  const auto intersection1_2d =
    intersection(space0_2d, space1_2d, affine::Equal_to_precision());
  EXPECT_PRED1(has_value_test, intersection0_1d);
  EXPECT_EQ(intersection0_1d->point(), space0_1d.point());
  EXPECT_PRED1(std::not_fn(has_value_test), intersection1_1d);
  EXPECT_PRED1(has_value_test, intersection0_2d);
  EXPECT_EQ(intersection0_2d->point(), space0_2d.point());
  EXPECT_PRED1(std::not_fn(has_value_test), intersection1_2d);
}

TEST(IntersectionTest, lineWithLine)
{
  const Point2D p0_2d{ 0.0, 1.0 };
  const Point2D p1_2d{ 1.0, 0.1 };
  const Point2D p2_2d{ 0.0, 0.0 };
  const Point2D p3_2d{ 1.0, 1.0 };
  const Point3D p0_3d{ 1.0, 0.0, 0.0 };
  const Point3D p1_3d{ 0.0, 1.0, 1.0 };
  const Point3D p2_3d{ 0.0, 0.0, 0.0 };
  const Point3D p3_3d{ 1.0, 1.0, 1.0 };
  const auto space0_2d = affine::Affine_space<double, 2>::spanning_space(
    { p0_2d, p1_2d }, affine::Equal_to_precision());
  const auto space1_2d = affine::Affine_space<double, 2>::spanning_space(
    { p2_2d, p3_2d }, affine::Equal_to_precision());
  const auto space0_3d = affine::Affine_space<double, 3>::spanning_space(
    { p0_3d, p1_3d }, affine::Equal_to_precision());
  const auto space1_3d = affine::Affine_space<double, 3>::spanning_space(
    { p2_3d, p3_3d }, affine::Equal_to_precision());

  const auto intersection_2d =
    intersection(space0_2d, space1_2d, affine::Equal_to_precision());
  const auto intersection_3d =
    intersection(space0_3d, space1_3d, affine::Equal_to_precision());
  EXPECT_PRED1(has_value_test, intersection_2d);
  EXPECT_PRED3(element_test,
               space0_2d,
               intersection_2d->point(),
               affine::Equal_to_precision());
  EXPECT_PRED3(element_test,
               space1_2d,
               intersection_2d->point(),
               affine::Equal_to_precision());
  EXPECT_EQ(intersection_2d->dimension(), 0);

  EXPECT_PRED1(has_value_test, intersection_3d);
  EXPECT_PRED3(element_test,
               space0_3d,
               intersection_3d->point(),
               affine::Equal_to_precision());
  EXPECT_PRED3(element_test,
               space1_3d,
               intersection_3d->point(),
               affine::Equal_to_precision());
  EXPECT_EQ(intersection_3d->dimension(), 0);
}

TEST(IntersectionTest, lineWithPlane)
{
  const Point3D p0{ 0.0, 0.0, 0.0 };
  const Point3D p1{ 1.0, 1.0, 1.0 };
  const Point3D p2{ 1.0, 0.0, 0.0 };
  const Point3D p3{ 0.0, 1.0, 0.0 };
  const Point3D p4{ 0.0, 0.0, 1.0 };
  const auto space0 = affine::Affine_space<double, 3>::spanning_space(
    { p0, p1 }, affine::Equal_to_precision());
  const auto space1 = affine::Affine_space<double, 3>::spanning_space(
    { p2, p3, p4 }, affine::Equal_to_precision());

  const auto intersection =
    affine::intersection(space0, space1, affine::Equal_to_precision());
  EXPECT_PRED1(has_value_test, intersection);
  EXPECT_PRED3(element_test,
               space0,
               intersection->point(),
               affine::Equal_to_precision());
  EXPECT_PRED3(element_test,
               space1,
               intersection->point(),
               affine::Equal_to_precision());
  EXPECT_EQ(intersection->dimension(), 0);
}

TEST(IntersectionTest, planeWithPlane)
{

  const Point3D p0{ 7.06976, 3.01336, 0.573289 };
  const Point3D p1{ 0.555447, 9.69536, -1.08835 };
  const Point3D p2{ -3.03373, 3.06169, -5.20502 };
  const Point3D p3{ -2.43233, -8.75036, 4.67991 };
  const Point3D p4{ 5.13406, -5.94389, 1.66657 };
  const Point3D p5{ -0.629713, 0.869835, -4.97337 };
  const auto space0 = affine::Affine_space<double, 3>::spanning_space(
    { p0, p1, p2 }, affine::Equal_to_precision());
  const auto space1 = affine::Affine_space<double, 3>::spanning_space(
    { p3, p4, p5 }, affine::Equal_to_precision());

  const auto intersection =
    affine::intersection(space0, space1, affine::Equal_to_precision());
  EXPECT_PRED1(has_value_test, intersection);
  EXPECT_PRED3(element_test,
               space0,
               intersection->point(),
               affine::Equal_to_precision(1e-14));
  EXPECT_PRED3(element_test,
               space0,
               intersection->point() + intersection->base(0),
               affine::Equal_to_precision(1e-14));
  EXPECT_PRED3(element_test,
               space1,
               intersection->point(),
               affine::Equal_to_precision(1e-14));
  EXPECT_PRED3(element_test,
               space1,
               intersection->point() + intersection->base(0),
               affine::Equal_to_precision(1e-14));
  EXPECT_EQ(intersection->dimension(), 1);
}

TEST(IntersectionTest, intervals)
{
  using IPoint2D = capd::vectalg::Vector<capd::DInterval, 2>;
  const auto ielement_test_2d = std::mem_fn(&affine::Affine_space<capd::DInterval, 2>::element<affine::IApprox_equal>);
  const auto ihas_value_test_2d = std::mem_fn(&std::optional<affine::Affine_space<capd::DInterval, 2>>::has_value);
  const IPoint2D p0{0.0, 1.0};
  const IPoint2D p1{1.0, 0.1};
  const IPoint2D p2{0.0, 0.0};
  const IPoint2D p3{1.0, 1.0};
  const auto space0 = affine::Affine_space<capd::DInterval, 2>::spanning_space({p0, p1}, affine::IApprox_equal());
  const auto space1 = affine::Affine_space<capd::DInterval, 2>::spanning_space({p2, p3}, affine::IApprox_equal());

  const auto intersection = affine::intersection(space0, space1, affine::IApprox_equal());
  EXPECT_PRED1(ihas_value_test_2d, intersection);
  EXPECT_PRED3(ielement_test_2d, space0, intersection->point(), affine::IApprox_equal());
  EXPECT_PRED3(ielement_test_2d, space1, intersection->point(), affine::IApprox_equal());
  EXPECT_EQ(intersection->dimension(), 0);
}