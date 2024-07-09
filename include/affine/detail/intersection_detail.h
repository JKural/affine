#pragma once

#include <optional>

#include "affine/affine_space.h"
#include "affine/detail/affine_detail.h"

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Woverloaded-virtual"
#include <capd/capdlib.h>
#pragma GCC diagnostic pop

namespace affine {

namespace detail {

template<class ScalarT, unsigned AmbientDim, class F>
std::optional<Affine_space<ScalarT, AmbientDim>>
space_with_space_intersection(
  const Affine_space<ScalarT, AmbientDim>& /* unused */,
  const Affine_space<ScalarT, AmbientDim>& /* unused */,
  F /* unused */)
{
  throw std::logic_error("Not implemented");
}

template<class ScalarT, unsigned AmbientDim, class F>
std::optional<Affine_space<ScalarT, AmbientDim>>
space_with_full_space_intersection(const Affine_space<ScalarT, AmbientDim>& lhs,
                                   const Affine_space<ScalarT, AmbientDim>& rhs,
                                   F /* unused */)
{
  assert(rhs.dimension() == rhs.ambient_dimension());
  return lhs;
}

template<class ScalarT, class F>
std::optional<Affine_space<ScalarT, 1>>
point_with_point_intersection(const Affine_space<ScalarT, 1>& lhs,
                              const Affine_space<ScalarT, 1>& rhs,
                              F compare)
{
  assert(lhs.dimension() == 0 && rhs.dimension() == 0);
  if (compare(lhs.point()[0], rhs.point()[0]))
    return lhs;
  else
    return std::nullopt;
}

template<class ScalarT, unsigned AmbientDim, class F>
std::optional<Affine_space<ScalarT, AmbientDim>>
point_with_point_intersection(const Affine_space<ScalarT, AmbientDim>& lhs,
                              const Affine_space<ScalarT, AmbientDim>& rhs,
                              F compare)
{
  assert(lhs.dimension() == 0 && rhs.dimension() == 0);
  if (compare(capd::vectalg::euclNorm(lhs.point() - rhs.point()), 0))
    return lhs;
  else
    return std::nullopt;
}

template<class ScalarT, unsigned AmbientDim, class F>
std::optional<Affine_space<ScalarT, AmbientDim>>
point_with_space_intersection(const Affine_space<ScalarT, AmbientDim>& lhs,
                              const Affine_space<ScalarT, AmbientDim>& rhs,
                              F compare)
{
  assert(lhs.dimension() == 0);
  if (rhs.element(lhs.point(), std::move(compare)))
    return lhs;
  else
    return std::nullopt;
}

template<class ScalarT, class F>
std::optional<Affine_space<ScalarT, 2>>
line_with_line_intersection(const Affine_space<ScalarT, 2>& lhs,
                            const Affine_space<ScalarT, 2>& rhs,
                            F compare)
{
  assert(lhs.dimension() == 1 && rhs.dimension() == 1);
  const auto& a = lhs.point();
  const auto& v = lhs.base(0);
  const auto& b = rhs.point();
  const auto& w = rhs.base(0);
  const auto det = -v[0] * w[1] + v[1] * w[0];
  const auto diff = b - a;
  if (compare(det, 0)) {
    if (compare(capd::vectalg::euclNorm(diff), 0))
      return lhs;
    else
      return std::nullopt;
  } else {
    const auto detx = -diff[0] * w[1] + diff[1] * w[0];
    return Affine_space<ScalarT, 2>(a + v * detx / det);
  }
}

template<class ScalarT, unsigned AmbientDim, class F>
std::optional<Affine_space<ScalarT, AmbientDim>>
line_with_line_intersection(const Affine_space<ScalarT, AmbientDim>& lhs,
                            const Affine_space<ScalarT, AmbientDim>& rhs,
                            F compare)
{
  assert(lhs.dimension() == 1 && rhs.dimension() == 1);
  const auto& w = rhs.base(0);
  auto v = lhs.base(0);
  const auto prod = v * w;
  v -= prod * w;
  const auto norm = capd::vectalg::euclNorm(v);
  if (compare(norm, 0)) {
    if (rhs.element(lhs.point(), std::move(compare)))
      return lhs;
    else
      return std::nullopt;
  } else {
    v /= norm;
    const auto diff = lhs.point() - rhs.point();
    const auto t = -(v * diff) / norm;
    const auto s = w * diff + prod * t;
    if (compare(capd::vectalg::euclNorm(w * s - lhs.base(0) * t - diff), 0)) {
      return Affine_space<ScalarT, AmbientDim>(lhs.point() + lhs.base(0) * t);
    } else {
      return std::nullopt;
    }
  }
}

template<class ScalarT, class F>
std::optional<Affine_space<ScalarT, 2>>
line_with_codim_1_space_intersection(const Affine_space<ScalarT, 2>& lhs,
                                     const Affine_space<ScalarT, 2>& rhs,
                                     F compare)
{
  return line_with_line_intersection(lhs, rhs, std::move(compare));
}

template<class ScalarT, unsigned AmbientDim, class F>
std::optional<Affine_space<ScalarT, AmbientDim>>
line_with_codim_1_space_intersection(
  const Affine_space<ScalarT, AmbientDim>& lhs,
  const Affine_space<ScalarT, AmbientDim>& rhs,
  F compare)
{
  assert(lhs.dimension() == 1 &&
         rhs.dimension() == rhs.ambient_dimension() - 1);
  auto v = lhs.base(0);
  for (auto&& w : rhs.base()) {
    v -= (v * w) * w;
  }
  if (compare(capd::vectalg::euclNorm(v), 0)) {
    if (rhs.element(lhs.point(), std::move(compare)))
      return lhs;
    else
      return std::nullopt;
  } else {
    v.normalize();
    const auto diff = lhs.point() - rhs.point();
    const auto t = -(diff * v) / (lhs.base(0) * v);
    return Affine_space<ScalarT, AmbientDim>(lhs.point() + lhs.base(0) * t);
  }
}

template<class ScalarT, class F>
std::optional<Affine_space<ScalarT, 3>>
plane_with_plane_intersection(const Affine_space<ScalarT, 3>& lhs,
                              const Affine_space<ScalarT, 3>& rhs,
                              F compare)
{
  assert(lhs.dimension() == 2 && rhs.dimension() == 2);
  const auto n_1 = detail::cross_product(lhs.base(0), lhs.base(1));
  const auto n_2 = detail::cross_product(rhs.base(0), rhs.base(1));
  const auto h_1 = lhs.point() * n_1;
  const auto h_2 = rhs.point() * n_2;
  const auto prod = n_1 * n_2;
  if (compare(prod * prod, 1)) {
    if (rhs.element(lhs.point(), std::move(compare)))
      return lhs;
    else
      return std::nullopt;
  } else {
    const auto c_1 = (h_1 - h_2 * prod) / (1 - prod * prod);
    const auto c_2 = (h_2 - h_1 * prod) / (1 - prod * prod);
    return Affine_space<ScalarT, 3>(
      c_1 * n_1 + c_2 * n_2,
      std::vector{ detail::cross_product(n_1, n_2) },
      std::move(compare));
  }
}

template<class ScalarT, unsigned AmbientDim, class F>
std::optional<Affine_space<ScalarT, AmbientDim>>
plane_with_plane_intersection(const Affine_space<ScalarT, AmbientDim>& lhs,
                              const Affine_space<ScalarT, AmbientDim>& rhs,
                              F compare)
{
  return space_with_space_intersection(lhs, rhs, std::move(compare));
}

} // namespace detail

} // namespace affine