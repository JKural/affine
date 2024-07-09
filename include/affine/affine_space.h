#pragma once

#include <concepts>
#include <initializer_list>
#include <ranges>
#include <span>
#include <stdexcept>
#include <vector>

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Woverloaded-virtual"
#include <capd/capdlib.h>
#pragma GCC diagnostic pop

#include "affine/detail/affine_detail.h"

namespace affine {

template<class ScalarT, unsigned AmbientDim>
class Affine_space
{
public:
  using Scalar = ScalarT;
  using Point = detail::Point<ScalarT, AmbientDim>;

  Affine_space() = default;

  Affine_space(const Point& particular_point)
    : particular_point_{ particular_point }
  {
  }

  template<std::ranges::input_range R, class F>
    requires std::convertible_to<std::ranges::range_value_t<R>, Point> &&
             std::ranges::sized_range<R>
  explicit Affine_space(const Point& particular_point,
                        R&& generators,
                        F compare)
    : particular_point_(particular_point)
  {
    namespace rs = std::ranges;
    base_.resize(rs::size(generators));
    rs::copy(generators, rs::begin(base_));
    const auto dim =
      detail::gram_schmidt_orthonormalization(base_, std::move(compare));
    base_.resize(dim);
  }

  template<class F>
  explicit Affine_space(const Point& particular_point,
                        std::vector<Point> generators,
                        F compare)
    : particular_point_(particular_point)
    , base_(std::move(generators))
  {
    const auto dim = detail::gram_schmidt_orthonormalization(base_, std::move(compare));
    base_.resize(dim);
  }

  constexpr std::size_t ambient_dimension() const { return AmbientDim; }

  std::size_t dimension() const { return base_.size(); }

  const Point& point() const { return particular_point_; }

  const std::vector<Point>& base() const { return base_; }

  const Point& base(std::size_t i) const { return base_.at(i); }

  template<class F>
  bool element(const Point& point, F compare) const
  {
    return detail::in_range_of_orthogonal_vectors(
      particular_point_ - point, base_, std::move(compare));
  }

  template<std::ranges::input_range R, class F>
    requires std::convertible_to<std::ranges::range_value_t<R>, Point> &&
             std::ranges::sized_range<R> &&
             (!std::same_as<R, std::initializer_list<Point>>)
  static Affine_space spanning_space(R&& points, F compare)
  {
    return spanning_space_impl(std::move(points), std::move(compare));
  }

  template<class F>
  static Affine_space spanning_space(std::initializer_list<Point> il, F compare)
  {
    return spanning_space_impl(std::move(il), std::move(compare));
  }

private:
  template<std::ranges::input_range R, class F>
    requires std::convertible_to<std::ranges::range_value_t<R>, Point> &&
             std::ranges::sized_range<R>
  static Affine_space spanning_space_impl(R&& points, F compare)
  {
    namespace rs = std::ranges;
    if (rs::empty(points))
      throw std::invalid_argument("'points' cannot be empty");
    const auto& head = *rs::begin(points);
    const auto rest = points | rs::views::drop(1);
    auto subtract_head = [&head](auto&& x) { return x - head; };
    return Affine_space(
      head, rest | rs::views::transform(subtract_head), std::move(compare));
  }

  Point particular_point_{};
  std::vector<Point> base_{};
};

template<class ScalarT, unsigned AmbientDim>
std::ostream&
operator<<(std::ostream& out,
           const Affine_space<ScalarT, AmbientDim>& affine_space)
{
  namespace rs = std::ranges;
  out << "[" << affine_space.point() << ", ";
  if (affine_space.base().empty())
    out << "{}]";
  else {
    out << "{" << affine_space.base(0);
    for (auto&& v : affine_space.base() | rs::views::drop(1))
      out << ", " << v;
    out << "}]";
  }
  return out;
}

} // namespace affine