#pragma once

#include <cmath>
#include <stdexcept>

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Woverloaded-virtual"
#include <capd/capdlib.h>
#pragma GCC diagnostic pop

namespace affine {

class Equal_to_precision
{
public:
  constexpr Equal_to_precision() = default;

  constexpr Equal_to_precision(double eps)
    : eps_{ eps }
  {
    if (eps_ <= 0.0)
      throw std::domain_error("precision must be positive real number");
  }

  template<class T, class U>
  constexpr bool operator()(const T& lhs, const U& rhs) const
  {
    if (lhs == rhs)
      return true;

    auto lhs_abs = capd::abs(lhs);
    auto rhs_abs = capd::abs(rhs);
    auto diff_abs = capd::abs(lhs - rhs);
    auto norm = std::min(lhs_abs + rhs_abs, std::numeric_limits<double>::max());
    return diff_abs < std::max(eps_, eps_ * norm);
  }

  constexpr double eps() const { return eps_; }

private:
  double eps_ = 1e-15;
};

class IApprox_equal
{
public:
  IApprox_equal() = default;

  template<class T, class U>
  bool operator()(const T& lhs, const U& rhs) const
  {
    if (lhs == rhs)
      return true;
    
    return capd::intervals::isSingular(lhs - rhs);
  }
};

} // namespace affine