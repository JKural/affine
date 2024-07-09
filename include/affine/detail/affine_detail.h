#pragma once

#include <algorithm>
#include <ranges>
#include <stdexcept>

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Woverloaded-virtual"
#include <capd/capdlib.h>
#pragma GCC diagnostic pop

namespace affine {

namespace detail {

template<class T, int N>
using Point = capd::vectalg::Vector<T, N>;

template<std::ranges::forward_range R, class F>
std::ptrdiff_t
gram_schmidt_orthonormalization(R&& vectors, F compare)
{
  namespace rs = std::ranges;
  for (auto outer_iter = rs::begin(vectors); outer_iter != rs::end(vectors);
       ++outer_iter) {
    auto& v_i = *outer_iter;
    const auto max_element = std::ranges::max_element(
      std::ranges::subrange(outer_iter, rs::end(vectors)), {}, [](auto&& x) {
        return x * x;
      });
    rs::iter_swap(max_element, outer_iter);
    if (compare(capd::vectalg::euclNorm(v_i), 0))
      return rs::distance(rs::begin(vectors), outer_iter);
    v_i.normalize();
    for (auto inner_iter = rs::next(outer_iter); inner_iter != rs::end(vectors);
         ++inner_iter) {
      auto& v_j = *inner_iter;
      v_j -= (v_j * v_i) * v_i;
    }
  }
  return rs::ssize(vectors);
}

template<class T, std::ranges::input_range R, class F>
bool
in_range_of_orthogonal_vectors(T point, R&& vectors, F compare)
{
  for (auto&& v : vectors)
    point -= (v * point) * v;
  const auto norm = capd::vectalg::euclNorm(point);
  return compare(norm, 0);
}

template<class ScalarT>
Point<ScalarT, 3>
cross_product(const Point<ScalarT, 3>& lhs, const Point<ScalarT, 3>& rhs)
{
  return Point<ScalarT, 3>{ lhs[1] * rhs[2] - lhs[2] * rhs[1],
                            lhs[2] * rhs[0] - lhs[0] * rhs[2],
                            lhs[0] * rhs[1] - lhs[1] * rhs[0] };
}

} // namespace detail

} // namespace affine
