#pragma once

#include <array>
#include <vector>
#include <cmath>
#include "../concepts/multivariate_concepts.hpp"
#include "../accumulators/accumulators.hpp"
#include "../quadrature/quadrature_rules.hpp"
#include "regions.hpp"

namespace calckit::multivariate {

// Result type for multivariate integration
template <typename T, int Dim>
struct integration_result_nd {
    T value{};           // Integral value
    T error{};          // Error estimate
    int dimension{Dim}; // Dimension
    size_t evaluations{}; // Function evaluations
    bool converged{false}; // Convergence flag

    // Default constructor
    integration_result_nd() = default;

    // Value constructor
    integration_result_nd(T val, T err, size_t evals, bool conv = true)
        : value(val), error(err), evaluations(evals), converged(conv) {}

    // Comparison
    friend bool operator==(const integration_result_nd&, const integration_result_nd&) = default;
};

// Tensor product integrator: extends 1D rule to multiple dimensions
template <typename T, int Dim, typename Rule1D, typename Accumulator>
class tensor_product_integrator {
public:
    using value_type = T;
    using rule_type = Rule1D;
    using accumulator_type = Accumulator;
    static constexpr int dimension_value = Dim;

private:
    Rule1D rule_1d;
    Accumulator acc;

public:
    // Constructor
    tensor_product_integrator(Rule1D rule = Rule1D{}, Accumulator accumulator = Accumulator{})
        : rule_1d(rule), acc(accumulator) {}

    // Integrate over hyperrectangle
    template <concepts::MultivariateFunction<T, Dim> F>
    integration_result_nd<T, Dim> operator()(
        F f,
        const hyperrectangle<T, Dim>& region
    ) {
        // Get 1D nodes and weights for each dimension
        auto nodes_1d = rule_1d.nodes();
        auto weights_1d = rule_1d.weights();
        size_t n = nodes_1d.size();

        // Total number of quadrature points = n^Dim
        size_t total_points = 1;
        for (int d = 0; d < Dim; ++d) {
            total_points *= n;
        }

        // Reset accumulator
        acc = Accumulator{};

        // Generate all tensor product nodes and weights
        integrate_recursive(f, region, nodes_1d, weights_1d, 0, {}, T(1));

        // Scale by Jacobian (region volume / reference volume)
        T jacobian = region.volume();
        T value = acc.result() * jacobian;

        // Error estimate (placeholder - needs refinement)
        T error = std::abs(value) * std::numeric_limits<T>::epsilon() * std::sqrt(total_points);

        return {value, error, total_points, true};
    }

    // Adaptive integration with subdivision
    template <concepts::MultivariateFunction<T, Dim> F>
    integration_result_nd<T, Dim> adaptive(
        F f,
        const hyperrectangle<T, Dim>& region,
        T tolerance = T(1e-8),
        size_t max_subdivisions = 1000
    ) {
        // Initial estimate on full region
        auto result = (*this)(f, region);

        if (result.error < tolerance) {
            return result;
        }

        // Subdivide and integrate
        return adaptive_recursive(f, region, tolerance, max_subdivisions, 0);
    }

private:
    // Recursive tensor product evaluation
    template <typename F>
    void integrate_recursive(
        F f,
        const hyperrectangle<T, Dim>& region,
        const std::vector<T>& nodes_1d,
        const std::vector<T>& weights_1d,
        int current_dim,
        std::array<T, Dim> point,
        T weight
    ) {
        if (current_dim == Dim) {
            // All dimensions filled - evaluate function
            acc.add(f(point) * weight);
            return;
        }

        // Iterate over quadrature points in current dimension
        size_t n = nodes_1d.size();
        T lower = region.lower[current_dim];
        T upper = region.upper[current_dim];

        for (size_t i = 0; i < n; ++i) {
            // Transform node from [-1,1] to [lower, upper]
            T xi = nodes_1d[i];
            T x = lower + (xi + T(1)) * (upper - lower) / T(2);

            // Update point and weight
            point[current_dim] = x;
            T new_weight = weight * weights_1d[i];

            // Recurse to next dimension
            integrate_recursive(f, region, nodes_1d, weights_1d, current_dim + 1, point, new_weight);
        }
    }

    // Adaptive recursive subdivision
    template <typename F>
    integration_result_nd<T, Dim> adaptive_recursive(
        F f,
        const hyperrectangle<T, Dim>& region,
        T tolerance,
        size_t max_subdivisions,
        size_t current_subdivisions
    ) {
        // Base case: reached subdivision limit
        if (current_subdivisions >= max_subdivisions) {
            return (*this)(f, region);
        }

        // Integrate on full region
        auto result_full = (*this)(f, region);

        // Subdivide along longest dimension
        auto [left, right] = region.subdivide();

        // Integrate on subregions
        auto result_left = (*this)(f, left);
        auto result_right = (*this)(f, right);

        // Combined result
        T value_combined = result_left.value + result_right.value;
        T error_combined = result_left.error + result_right.error;
        size_t evals_combined = result_left.evaluations + result_right.evaluations;

        // Check error difference
        T error_diff = std::abs(value_combined - result_full.value);

        if (error_diff < tolerance) {
            // Converged!
            return {value_combined, error_diff, evals_combined, true};
        }

        // Need more subdivision
        auto result_left_refined = adaptive_recursive(
            f, left, tolerance / T(2), max_subdivisions, current_subdivisions + 1
        );
        auto result_right_refined = adaptive_recursive(
            f, right, tolerance / T(2), max_subdivisions, current_subdivisions + 1
        );

        return {
            result_left_refined.value + result_right_refined.value,
            result_left_refined.error + result_right_refined.error,
            result_left_refined.evaluations + result_right_refined.evaluations,
            result_left_refined.converged && result_right_refined.converged
        };
    }
};

// Convenience factory functions
template <typename T, int Dim, typename Rule1D>
auto make_tensor_product_integrator(Rule1D rule = Rule1D{}) {
    using Acc = accumulators::kahan_accumulator<T>;
    return tensor_product_integrator<T, Dim, Rule1D, Acc>(rule, Acc{});
}

// 2D specialization with explicit interface
template <typename T, typename Rule1D, typename Accumulator = accumulators::kahan_accumulator<T>>
class integrator_2d {
private:
    tensor_product_integrator<T, 2, Rule1D, Accumulator> impl;

public:
    integrator_2d(Rule1D rule = Rule1D{}, Accumulator acc = Accumulator{})
        : impl(rule, acc) {}

    // Integrate over rectangle [x0, x1] × [y0, y1]
    template <typename F>
    auto operator()(F f, T x0, T x1, T y0, T y1) {
        hyperrectangle<T, 2> region({x0, y0}, {x1, y1});
        return impl(f, region);
    }

    // Integrate over rectangle with tolerance
    template <typename F>
    auto adaptive(F f, T x0, T x1, T y0, T y1, T tol = T(1e-8)) {
        hyperrectangle<T, 2> region({x0, y0}, {x1, y1});
        return impl.adaptive(f, region, tol);
    }

    // Integrate over rectangle region
    template <typename F>
    auto operator()(F f, const hyperrectangle<T, 2>& region) {
        return impl(f, region);
    }
};

// 3D specialization
template <typename T, typename Rule1D, typename Accumulator = accumulators::kahan_accumulator<T>>
class integrator_3d {
private:
    tensor_product_integrator<T, 3, Rule1D, Accumulator> impl;

public:
    integrator_3d(Rule1D rule = Rule1D{}, Accumulator acc = Accumulator{})
        : impl(rule, acc) {}

    // Integrate over box [x0, x1] × [y0, y1] × [z0, z1]
    template <typename F>
    auto operator()(F f, T x0, T x1, T y0, T y1, T z0, T z1) {
        hyperrectangle<T, 3> region({x0, y0, z0}, {x1, y1, z1});
        return impl(f, region);
    }

    // Integrate with tolerance
    template <typename F>
    auto adaptive(F f, T x0, T x1, T y0, T y1, T z0, T z1, T tol = T(1e-8)) {
        hyperrectangle<T, 3> region({x0, y0, z0}, {x1, y1, z1});
        return impl.adaptive(f, region, tol);
    }

    // Integrate over box region
    template <typename F>
    auto operator()(F f, const hyperrectangle<T, 3>& region) {
        return impl(f, region);
    }
};

// Convenience free functions
template <typename F, typename T>
auto integrate_2d(F f, T x0, T x1, T y0, T y1, T tol = T(1e-8)) {
    using Rule = quadrature::gauss_legendre_rule<T, 7>;  // 7-point rule
    integrator_2d<T, Rule> integrator;
    return integrator.adaptive(f, x0, x1, y0, y1, tol);
}

template <typename F, typename T>
auto integrate_3d(F f, T x0, T x1, T y0, T y1, T z0, T z1, T tol = T(1e-8)) {
    using Rule = quadrature::gauss_legendre_rule<T, 5>;  // 5-point rule (5^3 = 125 points)
    integrator_3d<T, Rule> integrator;
    return integrator.adaptive(f, x0, x1, y0, y1, z0, z1, tol);
}

} // namespace calckit::multivariate
