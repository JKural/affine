#pragma once

// skips -Woverloaded-virtual warnings for capd library
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Woverloaded-virtual"
#include <capd/capdlib.h>
#pragma GCC diagnostic pop

#include "affine/affine.h"

using Point1D = capd::vectalg::Vector<double, 1>;
using Point2D = capd::vectalg::Vector<double, 2>;
using Point3D = capd::vectalg::Vector<double, 3>;
using Point4D = capd::vectalg::Vector<double, 4>;

using IPoint1D = capd::vectalg::Vector<capd::DInterval, 1>;
using IPoint2D = capd::vectalg::Vector<capd::DInterval, 2>;
using IPoint3D = capd::vectalg::Vector<capd::DInterval, 3>;
using IPoint4D = capd::vectalg::Vector<capd::DInterval, 4>;

inline constexpr auto element_test =
  [](auto&& space, auto&& elem, auto&& pred) {
    return std::forward<decltype(space)>(space).element(
      std::forward<decltype(elem)>(elem), std::forward<decltype(pred)>(pred));
  };

inline constexpr auto has_value_test = [](auto&& opt) {
  return std::forward<decltype(opt)>(opt).has_value();
};