#pragma once

/**
 * @file methods.hpp
 * @brief Integration method objects for composable numerical integration.
 *
 * This module provides first-class method objects that encapsulate integration
 * strategies. Methods can be passed to `.eval()` on integrals to control how
 * the integration is performed.
 *
 * @section methods_usage Usage
 *
 * @code{.cpp}
 * using namespace limes::expr;
 * using namespace limes::methods;
 *
 * auto x = arg<0>;
 * auto I = integral(sin(x)).over<0>(0.0, 3.14159);
 *
 * // Gauss-Legendre quadrature (7 nodes)
 * auto r1 = I.eval(gauss<7>());
 *
 * // Adaptive integration with custom tolerance
 * auto r2 = I.eval(adaptive_method().with_tolerance(1e-12));
 *
 * // Monte Carlo with reproducible seed
 * auto r3 = I.eval(monte_carlo_method(100000).with_seed(42));
 *
 * // Simpson's rule
 * auto r4 = I.eval(simpson_method<100>());
 * @endcode
 *
 * @section methods_available Available Methods
 *
 * | Method | Factory | Description |
 * |--------|---------|-------------|
 * | gauss_legendre | `gauss<N>()` | N-point Gauss-Legendre quadrature |
 * | adaptive | `adaptive_method(tol)` | Adaptive subdivision until convergence |
 * | monte_carlo | `monte_carlo_method(n)` | Random sampling with n points |
 * | simpson | `simpson_method<N>()` | Simpson's rule with N subdivisions |
 * | trapezoidal | `trapezoidal_method<N>()` | Trapezoidal rule with N subdivisions |
 * | adaptive_composed | `make_adaptive(base, tol)` | Wrap any method with adaptive refinement |
 *
 * @defgroup methods Integration Methods
 * @brief Pluggable integration method objects.
 */

#include <cstddef>
#include <optional>
#include <random>
#include <functional>
#include "concepts.hpp"
#include "../algorithms/integrators/integrators.hpp"
#include "../algorithms/quadrature/quadrature.hpp"

/**
 * @namespace limes::methods
 * @brief Integration method objects.
 * @ingroup methods
 *
 * This namespace contains method objects that can be passed to `Integral::eval()`
 * to control the integration algorithm. Methods are composable and can be
 * configured with builder-style methods like `with_tolerance()` or `with_seed()`.
 */
namespace limes::methods {

// =============================================================================
// Gauss-Legendre Quadrature Method
// =============================================================================

/**
 * @brief Gauss-Legendre quadrature method with N nodes.
 * @ingroup methods
 *
 * Gauss-Legendre quadrature is optimal for smooth functions, achieving
 * degree 2N-1 polynomial exactness with only N function evaluations.
 *
 * @tparam N Number of quadrature nodes (higher = more accurate)
 * @tparam T Value type (default: double)
 *
 * @par Example
 * @code{.cpp}
 * auto I = integral(exp(x)).over<0>(0.0, 1.0);
 * auto result = I.eval(gauss<7>());  // 7-point quadrature
 * @endcode
 *
 * @see gauss() Factory function
 */
template<std::size_t N, typename T = double>
struct gauss_legendre {
    static constexpr std::size_t order = N;
    using value_type = T;

    constexpr gauss_legendre() noexcept = default;

    /// Integrate function f over [a, b]
    [[nodiscard]] algorithms::integration_result<T>
    operator()(std::function<T(T)> const& f, T a, T b) const {
        algorithms::quadrature_integrator<T, algorithms::quadrature::gauss_legendre<T, N>> integrator;
        return integrator(f, a, b);
    }
};

// Specialization: mark as quadrature method
template<std::size_t N, typename T>
struct is_quadrature_method<gauss_legendre<N, T>> : std::true_type {};

template<std::size_t N, typename T>
struct is_integration_method<gauss_legendre<N, T>> : std::true_type {};

/// Factory function for gauss-legendre quadrature
template<std::size_t N, typename T = double>
[[nodiscard]] constexpr auto gauss() {
    return gauss_legendre<N, T>{};
}

// =============================================================================
// Adaptive Integration Method
// =============================================================================

/**
 * @brief Adaptive integration method with configurable tolerance.
 * @ingroup methods
 *
 * Uses recursive interval subdivision until the error estimate falls below
 * the specified tolerance. This is the recommended default for most integrands,
 * especially those with localized features or varying smoothness.
 *
 * @tparam T Value type (default: double)
 *
 * @par Example
 * @code{.cpp}
 * auto I = integral(sin(x*x)).over<0>(0.0, 10.0);
 *
 * // Default tolerance (1e-10)
 * auto r1 = I.eval(adaptive_method());
 *
 * // Custom tolerance
 * auto r2 = I.eval(adaptive_method(1e-14));
 *
 * // Builder pattern
 * auto r3 = I.eval(adaptive_method().with_tolerance(1e-12).with_max_subdivisions(2000));
 * @endcode
 *
 * @see adaptive_method() Factory function
 */
template<typename T = double>
struct adaptive {
    using value_type = T;

    T tolerance;
    std::size_t max_subdivisions;

    constexpr adaptive(T tol = T(1e-10), std::size_t max_sub = 1000) noexcept
        : tolerance{tol}, max_subdivisions{max_sub} {}

    /// Integrate function f over [a, b]
    [[nodiscard]] algorithms::integration_result<T>
    operator()(std::function<T(T)> const& f, T a, T b) const {
        algorithms::adaptive_integrator<T> integrator;
        return integrator(f, a, b, tolerance);
    }

    /// Create a new method with different tolerance
    [[nodiscard]] constexpr adaptive with_tolerance(T tol) const noexcept {
        return adaptive{tol, max_subdivisions};
    }

    /// Create a new method with different max subdivisions
    [[nodiscard]] constexpr adaptive with_max_subdivisions(std::size_t max_sub) const noexcept {
        return adaptive{tolerance, max_sub};
    }
};

// Specialization: mark as adaptive method
template<typename T>
struct is_adaptive_method<adaptive<T>> : std::true_type {};

template<typename T>
struct is_integration_method<adaptive<T>> : std::true_type {};

/// Factory function for adaptive integration
template<typename T = double>
[[nodiscard]] constexpr auto adaptive_method(T tol = T(1e-10)) {
    return adaptive<T>{tol};
}

// =============================================================================
// Monte Carlo Integration Method
// =============================================================================

/**
 * @brief Monte Carlo integration using random sampling.
 * @ingroup methods
 *
 * Estimates the integral by averaging function values at random points.
 * The error decreases as O(1/√n) regardless of dimension, making this
 * method especially useful for high-dimensional integrals and irregular regions.
 *
 * @tparam T Value type (default: double)
 *
 * @par Example
 * @code{.cpp}
 * // 1D Monte Carlo
 * auto I = integral(f).over<0>(0.0, 1.0);
 * auto r = I.eval(monte_carlo_method(100000));
 *
 * // With seed for reproducibility
 * auto r2 = I.eval(monte_carlo_method(100000).with_seed(42));
 *
 * // N-dimensional box integration
 * auto J = integral(g).over_box({{0,1}, {0,1}, {0,1}});  // 3D unit cube
 * auto r3 = J.eval(monte_carlo_method(1000000));
 * @endcode
 *
 * @see monte_carlo_method() Factory function
 * @see BoxIntegral For N-dimensional Monte Carlo integration
 */
template<typename T = double>
struct monte_carlo {
    using value_type = T;

    std::size_t samples;
    std::optional<std::size_t> seed;

    constexpr monte_carlo(std::size_t n, std::optional<std::size_t> s = std::nullopt) noexcept
        : samples{n}, seed{s} {}

    /// Integrate function f over [a, b] using random sampling
    [[nodiscard]] algorithms::integration_result<T>
    operator()(std::function<T(T)> const& f, T a, T b) const {
        std::mt19937_64 rng;
        if (seed) {
            rng.seed(*seed);
        } else {
            std::random_device rd;
            rng.seed(rd());
        }

        std::uniform_real_distribution<T> dist(a, b);

        T sum = T(0);
        T sum_sq = T(0);

        for (std::size_t i = 0; i < samples; ++i) {
            T x = dist(rng);
            T y = f(x);
            sum += y;
            sum_sq += y * y;
        }

        T interval_length = b - a;
        T mean = sum / T(samples);
        T variance = (sum_sq / T(samples) - mean * mean) / T(samples);
        T value = mean * interval_length;
        T error = std::sqrt(variance) * interval_length;

        algorithms::integration_result<T> result{value, error, samples, samples};
        result.variance_ = variance;
        return result;
    }

    /// Create method with specified seed for reproducibility
    [[nodiscard]] constexpr monte_carlo with_seed(std::size_t s) const noexcept {
        return monte_carlo{samples, s};
    }

    /// Create method with different sample count
    [[nodiscard]] constexpr monte_carlo with_samples(std::size_t n) const noexcept {
        return monte_carlo{n, seed};
    }
};

// Specialization: mark as Monte Carlo method
template<typename T>
struct is_monte_carlo_method<monte_carlo<T>> : std::true_type {};

template<typename T>
struct is_integration_method<monte_carlo<T>> : std::true_type {};

/// Factory function for Monte Carlo integration
template<typename T = double>
[[nodiscard]] constexpr auto monte_carlo_method(std::size_t n) {
    return monte_carlo<T>{n};
}

// =============================================================================
// Simpson's Rule Method
// =============================================================================

/**
 * @brief Simpson's rule with N subdivisions.
 * @ingroup methods
 *
 * Classic Simpson's 1/3 rule using parabolic interpolation.
 * Achieves O(h⁴) convergence for smooth functions.
 *
 * @tparam N Number of subdivisions (must be even)
 * @tparam T Value type (default: double)
 *
 * @par Example
 * @code{.cpp}
 * auto I = integral(exp(x)).over<0>(0.0, 1.0);
 * auto result = I.eval(simpson_method<100>());
 * @endcode
 *
 * @see simpson_method() Factory function
 */
template<std::size_t N, typename T = double>
struct simpson {
    static_assert(N % 2 == 0, "Simpson's rule requires even number of subdivisions");

    static constexpr std::size_t subdivisions = N;
    using value_type = T;

    constexpr simpson() noexcept = default;

    /// Integrate function f over [a, b]
    [[nodiscard]] algorithms::integration_result<T>
    operator()(std::function<T(T)> const& f, T a, T b) const {
        T h = (b - a) / T(N);
        T sum = f(a) + f(b);

        // Sum odd indices (coefficient 4)
        for (std::size_t i = 1; i < N; i += 2) {
            sum += T(4) * f(a + T(i) * h);
        }

        // Sum even indices (coefficient 2)
        for (std::size_t i = 2; i < N; i += 2) {
            sum += T(2) * f(a + T(i) * h);
        }

        T value = sum * h / T(3);

        // Error estimate (rough): O(h^4)
        T error = std::abs(value) * std::pow(h, 4);

        return algorithms::integration_result<T>{value, error, N + 1, N + 1};
    }
};

// Specialization: mark as integration method
template<std::size_t N, typename T>
struct is_integration_method<simpson<N, T>> : std::true_type {};

/// Factory function for Simpson's rule
template<std::size_t N, typename T = double>
[[nodiscard]] constexpr auto simpson_method() {
    return simpson<N, T>{};
}

// =============================================================================
// Trapezoidal Rule Method
// =============================================================================

/**
 * @brief Trapezoidal rule with N subdivisions.
 * @ingroup methods
 *
 * Classic trapezoidal rule using linear interpolation.
 * Simple but effective for periodic functions over complete periods.
 * Achieves O(h²) convergence for smooth functions.
 *
 * @tparam N Number of subdivisions
 * @tparam T Value type (default: double)
 *
 * @par Example
 * @code{.cpp}
 * auto I = integral(sin(x)).over<0>(0.0, 2*pi);
 * auto result = I.eval(trapezoidal_method<100>());
 * @endcode
 *
 * @see trapezoidal_method() Factory function
 */
template<std::size_t N, typename T = double>
struct trapezoidal {
    static constexpr std::size_t subdivisions = N;
    using value_type = T;

    constexpr trapezoidal() noexcept = default;

    /// Integrate function f over [a, b]
    [[nodiscard]] algorithms::integration_result<T>
    operator()(std::function<T(T)> const& f, T a, T b) const {
        T h = (b - a) / T(N);
        T sum = (f(a) + f(b)) / T(2);

        for (std::size_t i = 1; i < N; ++i) {
            sum += f(a + T(i) * h);
        }

        T value = sum * h;

        // Error estimate (rough): O(h^2)
        T error = std::abs(value) * h * h;

        return algorithms::integration_result<T>{value, error, N + 1, N + 1};
    }
};

// Specialization: mark as integration method
template<std::size_t N, typename T>
struct is_integration_method<trapezoidal<N, T>> : std::true_type {};

/// Factory function for trapezoidal rule
template<std::size_t N, typename T = double>
[[nodiscard]] constexpr auto trapezoidal_method() {
    return trapezoidal<N, T>{};
}

// =============================================================================
// Composed Adaptive Method (wraps another method for adaptive refinement)
// =============================================================================

/**
 * @brief Compose an adaptive scheme around any base method.
 * @ingroup methods
 *
 * Wraps another integration method with adaptive interval subdivision.
 * The base method is used for initial estimates, and the interval is
 * recursively subdivided until the error estimate is below tolerance.
 *
 * @tparam BaseMethod The wrapped integration method
 * @tparam T Value type (default: double)
 *
 * @par Example
 * @code{.cpp}
 * // Adaptive Simpson's rule
 * auto method = make_adaptive(simpson_method<10>(), 1e-10);
 * auto result = I.eval(method);
 *
 * // Adaptive Gauss-Legendre
 * auto method2 = make_adaptive(gauss<5>(), 1e-12);
 * @endcode
 *
 * @see make_adaptive() Factory function
 */
template<typename BaseMethod, typename T = double>
struct adaptive_composed {
    using value_type = T;

    BaseMethod base;
    T tolerance;
    std::size_t max_depth;

    constexpr adaptive_composed(BaseMethod m, T tol = T(1e-10), std::size_t depth = 20) noexcept
        : base{m}, tolerance{tol}, max_depth{depth} {}

    /// Integrate function f over [a, b] with adaptive refinement
    [[nodiscard]] algorithms::integration_result<T>
    operator()(std::function<T(T)> const& f, T a, T b) const {
        return integrate_adaptive(f, a, b, 0);
    }

    [[nodiscard]] constexpr adaptive_composed with_tolerance(T tol) const noexcept {
        return adaptive_composed{base, tol, max_depth};
    }

private:
    [[nodiscard]] algorithms::integration_result<T>
    integrate_adaptive(std::function<T(T)> const& f, T a, T b, std::size_t depth) const {
        auto full = base(f, a, b);

        if (depth >= max_depth) {
            full.converged_ = false;
            return full;
        }

        T mid = (a + b) / T(2);
        auto left = base(f, a, mid);
        auto right = base(f, mid, b);
        auto refined = left + right;

        T error = std::abs(full.value() - refined.value());

        if (error < tolerance) {
            refined.error_ = error;
            return refined;
        }

        // Recursive refinement
        auto left_refined = integrate_adaptive(f, a, mid, depth + 1);
        auto right_refined = integrate_adaptive(f, mid, b, depth + 1);

        return left_refined + right_refined;
    }
};

// Specialization: mark as adaptive method
template<typename M, typename T>
struct is_adaptive_method<adaptive_composed<M, T>> : std::true_type {};

template<typename M, typename T>
struct is_integration_method<adaptive_composed<M, T>> : std::true_type {};

/// Factory function for adaptive composition
template<typename M, typename T = double>
[[nodiscard]] constexpr auto make_adaptive(M base_method, T tol = T(1e-10)) {
    return adaptive_composed<M, T>{base_method, tol};
}

// =============================================================================
// Default method
// =============================================================================

/// The default integration method (adaptive Gauss-Legendre)
template<typename T = double>
using default_method = adaptive<T>;

} // namespace limes::methods
