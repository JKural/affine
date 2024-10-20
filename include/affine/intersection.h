#pragma once

#include <algorithm>
#include <functional>
#include <optional>

#include "affine/affine_space.h"
#include "affine/detail/intersection_detail.h"

namespace affine {

template<class ScalarT, unsigned AmbientDim, class F>
std::optional<Affine_space<ScalarT, AmbientDim>>
intersection(const Affine_space<ScalarT, AmbientDim>& lhs,
             const Affine_space<ScalarT, AmbientDim>& rhs,
             F compare)
{
  namespace rs = std::ranges;
  using Aff_space = Affine_space<ScalarT, AmbientDim>;
  const auto& [first, second] =
    rs::minmax(lhs, rhs, {}, std::mem_fn(&Aff_space::dimension));
  if (second.dimension() == AmbientDim)
    return detail::space_with_full_space_intersection(
      first, second, std::move(compare));
  switch (first.dimension()) {
    case 0: {
      switch (second.dimension()) {
        case 0:
          return detail::point_with_point_intersection(
            first, second, std::move(compare));
        default:
          return detail::point_with_space_intersection(
            first, second, std::move(compare));
      }
    }
    case 1: {
      if (second.dimension() == AmbientDim - 1)
        return detail::line_with_codim_1_space_intersection(
          first, second, std::move(compare));
      switch (second.dimension()) {
        case 1:
          return detail::line_with_line_intersection(
            first, second, std::move(compare));
        default:
          return detail::space_with_space_intersection(
            first, second, std::move(compare));
      }
    }
    case 2: {
      switch (second.dimension()) {
        case 2:
          return detail::plane_with_plane_intersection(
            first, second, std::move(compare));
        default:
          return detail::space_with_space_intersection(
            first, second, std::move(compare));
      }
    }
    default:
      return detail::space_with_space_intersection(
        first, second, std::move(compare));
  }
}

} // namespace affine
