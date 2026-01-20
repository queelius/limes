#pragma once

/**
 * @file limes.hpp
 * @brief Main entry point for the limes library.
 *
 * limes is a C++20 header-only library for composable calculus expressions
 * with symbolic differentiation and numerical integration.
 *
 * @section quick_start Quick Start
 *
 * @code{.cpp}
 * #include <limes/limes.hpp>
 * using namespace limes::expr;
 *
 * auto x = arg<0>;
 * auto f = sin(x * x);                       // f(x) = sin(x²)
 * auto df = derivative(f).wrt<0>();          // df/dx = 2x·cos(x²)
 *
 * auto I = integral(f).over<0>(0.0, 1.0);    // ∫₀¹ sin(x²) dx
 * auto result = I.eval();                    // ≈ 0.3103
 * @endcode
 *
 * @section namespaces Namespaces
 *
 * - `limes::expr` - Expression layer (user-facing API)
 * - `limes::methods` - Integration method objects
 * - `limes::algorithms` - Low-level numerical backend
 *
 * @author Alex Towell (queelius@gmail.com)
 * @see https://metafunctor.com
 * @copyright MIT License
 */

#include "fwd.hpp"
#include "algorithms/algorithms.hpp"
#include "expr/expr.hpp"

/**
 * @namespace limes
 * @brief Root namespace for the limes library.
 */

/// Namespace alias for convenience
namespace li = limes;
