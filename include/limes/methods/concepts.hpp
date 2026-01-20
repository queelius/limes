#pragma once

#include <concepts>
#include <functional>
#include <span>
#include <array>
#include "../algorithms/core/result.hpp"

namespace limes::methods {

// =============================================================================
// Method Concepts
// =============================================================================

/// Concept for 1D integration methods
/// A method must be callable with (function, lower, upper) and return integration_result
template<typename M, typename T>
concept IntegrationMethod = requires(M const& m, std::function<T(T)> f, T a, T b) {
    { m(f, a, b) } -> std::convertible_to<algorithms::integration_result<T>>;
};

/// Concept for methods that can be configured with tolerance
template<typename M, typename T>
concept ToleranceAware = requires(M const& m, T tol) {
    { m.with_tolerance(tol) } -> std::same_as<M>;
};

/// Concept for adaptive methods (that refine until convergence)
template<typename M>
concept Adaptive = requires(M const& m) {
    { m.tolerance } -> std::convertible_to<double>;
    { m.max_subdivisions } -> std::convertible_to<std::size_t>;
};

/// Concept for Monte Carlo methods (that use random sampling)
template<typename M>
concept MonteCarlo = requires(M const& m) {
    { m.samples } -> std::convertible_to<std::size_t>;
};

/// Concept for quadrature methods (fixed nodes and weights)
template<typename M>
concept Quadrature = requires(M const& m) {
    { M::order } -> std::convertible_to<std::size_t>;
};

// =============================================================================
// Type traits for method detection
// =============================================================================

template<typename M>
struct is_integration_method : std::false_type {};

template<typename M>
inline constexpr bool is_integration_method_v = is_integration_method<M>::value;

template<typename M>
struct is_adaptive_method : std::false_type {};

template<typename M>
inline constexpr bool is_adaptive_method_v = is_adaptive_method<M>::value;

template<typename M>
struct is_monte_carlo_method : std::false_type {};

template<typename M>
inline constexpr bool is_monte_carlo_method_v = is_monte_carlo_method<M>::value;

template<typename M>
struct is_quadrature_method : std::false_type {};

template<typename M>
inline constexpr bool is_quadrature_method_v = is_quadrature_method<M>::value;

} // namespace limes::methods
