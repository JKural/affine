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
using Point4D = capd::vectalg::Vector<double, 4>;

bool
implicit_convertion_test(const affine::Affine_space<double, 4>& /* unused */)
{
  return true;
}

} // namespace

TEST(ConstructorTest, defaultConstructor)
{
  const affine::Affine_space<double, 1> space_1d;
  const affine::Affine_space<double, 3> space_3d;

  EXPECT_EQ(space_1d.dimension(), 0);
  EXPECT_EQ(space_1d.ambient_dimension(), 1);
  EXPECT_EQ(space_1d.point(), Point1D());
  EXPECT_PRED1(std::mem_fn(&std::vector<Point1D>::empty), space_1d.base());

  EXPECT_EQ(space_3d.dimension(), 0);
  EXPECT_EQ(space_3d.ambient_dimension(), 3);
  EXPECT_EQ(space_3d.point(), Point3D());
  EXPECT_PRED1(std::mem_fn(&std::vector<Point3D>::empty), space_3d.base());
}

TEST(ConstructorTest, fromPoint)
{
  const Point2D p_2d{ 1.0, 2.0 };
  const affine::Affine_space space_1d(Point1D{ 1.0 });
  const affine::Affine_space space_2d(p_2d);
  const affine::Affine_space space_3d = Point3D{ 1.0, 2.0, 3.0 };

  EXPECT_EQ(space_1d.dimension(), 0);
  EXPECT_EQ(space_1d.ambient_dimension(), 1);
  EXPECT_EQ(space_1d.point(), Point1D{ 1.0 });
  EXPECT_PRED1(std::mem_fn(&std::vector<Point1D>::empty), space_1d.base());

  EXPECT_EQ(space_2d.dimension(), 0);
  EXPECT_EQ(space_2d.ambient_dimension(), 2);
  EXPECT_EQ(space_2d.point(), p_2d);
  EXPECT_PRED1(std::mem_fn(&std::vector<Point2D>::empty), space_2d.base());

  EXPECT_EQ(space_3d.dimension(), 0);
  EXPECT_EQ(space_3d.ambient_dimension(), 3);
  EXPECT_EQ(space_3d.point(), (Point3D{ 1.0, 2.0, 3.0 }));
  EXPECT_PRED1(std::mem_fn(&std::vector<Point3D>::empty), space_3d.base());

  EXPECT_TRUE(implicit_convertion_test(Point4D{ 1.0, 2.0, 3.0, 4.0 }));
}

TEST(ConstructorTest, fromRange)
{
  const Point1D p_1d{ 0.0 };
  const std::array<Point1D, 3> a_1d{ Point1D{ 1.0 },
                                     Point1D{ 2.0 },
                                     Point1D{ 3.0 } };
  const affine::Affine_space space_1d(p_1d, a_1d, affine::Equal_to_precision());
  const Point4D p_4d{ 0.491742, -4.73389, -6.07428, 0.246362 };
  const std::list<Point4D> l_4d{
    { -4.85797, 6.30975, -0.959427, -5.05184 },
    { -3.53352, 6.14748, 6.88908, -8.43334 },
    { -7.37983, -0.708072, -0.999901, -9.75365 },
  };
  const affine::Affine_space space_4d(p_4d, l_4d, affine::Equal_to_precision());

  EXPECT_EQ(space_1d.dimension(), 1);
  EXPECT_EQ(space_1d.ambient_dimension(), 1);
  EXPECT_EQ(space_1d.point(), p_1d);
  EXPECT_DOUBLE_EQ(capd::vectalg::euclNorm(space_1d.base(0)), 1);
  EXPECT_THROW(space_1d.base(1), std::out_of_range);

  EXPECT_EQ(space_4d.dimension(), 3);
  EXPECT_EQ(space_4d.ambient_dimension(), 4);
  EXPECT_EQ(space_4d.point(), p_4d);
  EXPECT_DOUBLE_EQ(capd::vectalg::euclNorm(space_4d.base(0)), 1);
  EXPECT_DOUBLE_EQ(capd::vectalg::euclNorm(space_4d.base(1)), 1);
  EXPECT_DOUBLE_EQ(capd::vectalg::euclNorm(space_4d.base(2)), 1);
  EXPECT_THROW(space_4d.base(3), std::out_of_range);
  EXPECT_NEAR(space_4d.base(0) * space_4d.base(1), 0, 1e-14);
  EXPECT_NEAR(space_4d.base(0) * space_4d.base(2), 0, 1e-14);
  EXPECT_NEAR(space_4d.base(1) * space_4d.base(2), 0, 1e-14);
}

TEST(ConstructorTest, fromVector)
{
  const Point4D p_4d{ 0.491742, -4.73389, -6.07428, 0.246362 };
  std::vector<Point4D> v_4d{
    { -4.85797, 6.30975, -0.959427, -5.05184 },
    { -3.53352, 6.14748, 6.88908, -8.43334 },
    { -7.37983, -0.708072, -0.999901, -9.75365 },
  };
  const affine::Affine_space space_4d(
    p_4d, std::move(v_4d), affine::Equal_to_precision());

  EXPECT_PRED1(std::mem_fn(&std::vector<Point4D>::empty), v_4d);
  EXPECT_EQ(space_4d.dimension(), 3);
  EXPECT_EQ(space_4d.ambient_dimension(), 4);
  EXPECT_EQ(space_4d.point(), p_4d);
  EXPECT_DOUBLE_EQ(capd::vectalg::euclNorm(space_4d.base(0)), 1);
  EXPECT_DOUBLE_EQ(capd::vectalg::euclNorm(space_4d.base(1)), 1);
  EXPECT_DOUBLE_EQ(capd::vectalg::euclNorm(space_4d.base(2)), 1);
  EXPECT_THROW(space_4d.base(3), std::out_of_range);
  EXPECT_NEAR(space_4d.base(0) * space_4d.base(1), 0, 1e-14);
  EXPECT_NEAR(space_4d.base(0) * space_4d.base(2), 0, 1e-14);
  EXPECT_NEAR(space_4d.base(1) * space_4d.base(2), 0, 1e-14);
}