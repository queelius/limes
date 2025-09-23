#pragma once

#include <thread>
#include <string_view>

// Core components
#include "concepts/integrator_concepts.hpp"
#include "core/integration_result.hpp"

// Accumulators for precision control
#include "accumulators/accumulators.hpp"

// Quadrature rules
#include "quadrature/quadrature_rules.hpp"

// Integration algorithms
#include "integrators/univariate_integrator.hpp"

// Parallel execution
#include "parallel/parallel_integration.hpp"

// Coordinate transforms
#include "transforms/coordinate_transforms.hpp"

// Main namespace alias
namespace ai = algebraic_integrators;

namespace algebraic_integrators {

// High-level integration interface with sensible defaults
template<concepts::Field T = double>
class integrate {
public:
    using value_type = T;
    using result_type = integration_result<T>;

    // Simple adaptive integration with automatic method selection
    template<concepts::UnivariateFunction<T> F>
    static result_type adaptive(F&& f, T a, T b, T tol = default_tolerance()) {
        if (std::isinf(a) || std::isinf(b)) {
            // Use tanh-sinh for infinite intervals
            auto integrator = make_tanh_sinh<T>();
            return integrator(std::forward<F>(f), a, b, tol);
        } else {
            // Use adaptive Gauss-Kronrod for finite intervals
            auto integrator = make_adaptive_integrator<T>();
            return integrator(std::forward<F>(f), a, b, tol);
        }
    }

    // Robust integration with fallback methods
    template<concepts::UnivariateFunction<T> F>
    static result_type robust(F&& f, T a, T b, T tol = default_tolerance()) {
        // Try Gauss-Kronrod first
        auto gk = make_adaptive_integrator<T>();
        auto result = gk(f, a, b, tol);

        if (!result.converged() || result.error() > tol) {
            // Fallback to Romberg
            auto romberg = make_romberg<T>();
            auto romberg_result = romberg(f, a, b, tol);

            // Choose better result
            if (romberg_result.error() < result.error()) {
                result = romberg_result;
            }
        }

        return result;
    }

    // High-precision integration with compensated summation
    template<concepts::UnivariateFunction<T> F>
    static result_type precise(F&& f, T a, T b, T tol = precise_tolerance()) {
        using acc = accumulators::klein_accumulator<T>;
        using rule = quadrature::gauss_kronrod_15<T>;

        quadrature_integrator<T, rule, acc> integrator{rule{}, acc{}};
        return integrator(std::forward<F>(f), a, b, tol);
    }

    // Parallel integration for expensive functions
    template<concepts::UnivariateFunction<T> F>
    static result_type parallel(F&& f, T a, T b, T tol = default_tolerance()) {
        auto integrator = parallel::make_parallel_integrator<T>();
        return integrator(std::forward<F>(f), a, b, tol);
    }

    // Monte Carlo integration
    static auto monte_carlo(std::size_t samples = 1000000) {
        return parallel::make_parallel_monte_carlo<T>(samples);
    }

    // Transform-based integration for special cases
    template<concepts::UnivariateFunction<T> F, typename Transform>
    static result_type with_transform(
        F&& f, T a, T b, Transform&& transform, T tol = default_tolerance()
    ) {
        auto transformed_f = [&](T t) -> T {
            T x = transform.forward(t);
            return f(x) * transform.jacobian(t);
        };

        T ta = transform.inverse(a);
        T tb = transform.inverse(b);

        return adaptive(transformed_f, ta, tb, tol);
    }

    // Convenience method for oscillatory integrands
    template<concepts::UnivariateFunction<T> F>
    static result_type oscillatory(F&& f, T a, T b, T frequency, T tol = default_tolerance()) {
        // Use Filon's method or adaptive subdivision
        // This is a simplified implementation
        std::size_t subdivisions = std::max(
            std::size_t(10),
            std::size_t(frequency * (b - a) / std::numbers::pi_v<T>)
        );

        T h = (b - a) / T(subdivisions);
        result_type total{};

        for (std::size_t i = 0; i < subdivisions; ++i) {
            T ai = a + i * h;
            T bi = (i == subdivisions - 1) ? b : ai + h;
            total += adaptive(f, ai, bi, tol / subdivisions);
        }

        return total;
    }

private:
    static constexpr T default_tolerance() {
        return T(1e-8);
    }

    static constexpr T precise_tolerance() {
        return T(100) * std::numeric_limits<T>::epsilon();
    }
};

// Convenience free functions
template<typename F, typename T>
auto integrate_adaptive(F&& f, T a, T b, T tol = T(1e-8)) {
    return integrate<T>::adaptive(std::forward<F>(f), a, b, tol);
}

template<typename F, typename T>
auto integrate_robust(F&& f, T a, T b, T tol = T(1e-8)) {
    return integrate<T>::robust(std::forward<F>(f), a, b, tol);
}

template<typename F, typename T>
auto integrate_precise(F&& f, T a, T b, T tol = T(1e-12)) {
    return integrate<T>::precise(std::forward<F>(f), a, b, tol);
}

template<typename F, typename T>
auto integrate_parallel(F&& f, T a, T b, T tol = T(1e-8)) {
    return integrate<T>::parallel(std::forward<F>(f), a, b, tol);
}

// Builder pattern for custom integrators
template<concepts::Field T = double>
class integrator_builder {
public:
    using value_type = T;

    integrator_builder& with_quadrature(std::string_view name) {
        quadrature_name_ = name;
        return *this;
    }

    integrator_builder& with_accumulator(std::string_view name) {
        accumulator_name_ = name;
        return *this;
    }

    integrator_builder& with_parallel(bool enable = true) {
        use_parallel_ = enable;
        return *this;
    }

    integrator_builder& with_threads(std::size_t n) {
        num_threads_ = n;
        return *this;
    }

    integrator_builder& with_tolerance(T tol) {
        tolerance_ = tol;
        return *this;
    }

    template<typename F>
    auto integrate(F&& f, T a, T b) const {
        if (use_parallel_) {
            parallel::execution_policy policy;
            policy.num_threads = num_threads_;
            auto integrator = parallel::make_parallel_integrator<T>(policy);
            return integrator(std::forward<F>(f), a, b, tolerance_);
        } else {
            return algebraic_integrators::integrate<T>::adaptive(std::forward<F>(f), a, b, tolerance_);
        }
    }

private:
    std::string_view quadrature_name_ = "gauss-kronrod";
    std::string_view accumulator_name_ = "neumaier";
    bool use_parallel_ = false;
    std::size_t num_threads_ = std::thread::hardware_concurrency();
    T tolerance_ = T(1e-8);
};

// Create a builder
template<concepts::Field T = double>
auto make_integrator() {
    return integrator_builder<T>{};
}

} // namespace algebraic_integrators