#pragma once

#include <concepts>
#include <array>
#include <ranges>
#include "integrator_concepts.hpp"

namespace calckit::concepts {

// Multivariate function: f: R^Dim → R
template <typename F, typename T, int Dim>
concept MultivariateFunction = requires(F f, std::array<T, Dim> x) {
    { f(x) } -> std::convertible_to<T>;
};

// Vector field: F: R^Dim → R^Dim
template <typename F, typename T, int Dim>
concept VectorField = requires(F f, std::array<T, Dim> x) {
    { f(x) } -> std::convertible_to<std::array<T, Dim>>;
};

// Integration region in R^Dim
template <typename R, typename T, int Dim>
concept IntegrationRegion = requires(R region) {
    // Dimension query
    { region.dimension() } -> std::same_as<int>;
    requires region.dimension() == Dim;

    // Volume/measure
    { region.volume() } -> std::convertible_to<T>;

    // Bounding box
    { region.lower_bounds() } -> std::convertible_to<std::array<T, Dim>>;
    { region.upper_bounds() } -> std::convertible_to<std::array<T, Dim>>;

    // Point containment
    { region.contains(std::declval<std::array<T, Dim>>()) } -> std::convertible_to<bool>;
};

// Region that supports uniform sampling (for Monte Carlo)
template <typename R, typename T, int Dim>
concept SampleableRegion = IntegrationRegion<R, T, Dim> && requires(R region) {
    // Generate random point uniformly in region
    { region.sample_uniform() } -> std::convertible_to<std::array<T, Dim>>;
};

// Multivariate quadrature rule: provides nodes and weights for integration
template <typename Q, typename T, int Dim>
concept MultivariateQuadratureRule = requires(Q rule) {
    typename Q::value_type;
    requires std::same_as<typename Q::value_type, T>;

    // Nodes (quadrature points)
    { rule.nodes() } -> std::ranges::range;
    requires std::same_as<
        std::ranges::range_value_t<decltype(rule.nodes())>,
        std::array<T, Dim>
    >;

    // Weights corresponding to nodes
    { rule.weights() } -> std::ranges::range;
    requires std::same_as<
        std::ranges::range_value_t<decltype(rule.weights())>,
        T
    >;

    // Rule order/precision
    { rule.order() } -> std::convertible_to<int>;

    // Number of quadrature points
    { rule.size() } -> std::convertible_to<size_t>;
};

// Tensor product of 1D rules
template <typename T1D, typename T, int Dim>
concept TensorProductCompatible = requires(T1D rule) {
    typename T1D::value_type;
    requires std::same_as<typename T1D::value_type, T>;

    // Must be a 1D quadrature rule
    { rule.nodes() } -> std::ranges::range;
    { rule.weights() } -> std::ranges::range;
};

// Multivariate integration result
template <typename R, typename T, int Dim>
concept MultivariateIntegrationResult = requires(R result) {
    // Value of integral
    { result.value } -> std::convertible_to<T>;

    // Error estimate
    { result.error } -> std::convertible_to<T>;

    // Dimension
    { result.dimension } -> std::convertible_to<int>;
    requires result.dimension == Dim;

    // Function evaluations
    { result.evaluations } -> std::convertible_to<size_t>;

    // Convergence flag
    { result.converged } -> std::convertible_to<bool>;
};

// Coordinate transformation: bijection R^Dim → R^Dim
template <typename T, typename Src, typename Dst, int Dim>
concept CoordinateTransform = requires(T transform, std::array<Src, Dim> x) {
    // Forward map: source → destination
    { transform.forward(x) } -> std::convertible_to<std::array<Dst, Dim>>;

    // Inverse map: destination → source
    { transform.inverse(std::declval<std::array<Dst, Dim>>()) }
        -> std::convertible_to<std::array<Src, Dim>>;

    // Jacobian determinant (for change of variables)
    { transform.jacobian(x) } -> std::convertible_to<Src>;
};

// Diffeomorphism: smooth invertible map with Jacobian
template <typename T, typename V, int Dim>
concept Diffeomorphism = CoordinateTransform<T, V, V, Dim> && requires(T transform, std::array<V, Dim> x) {
    // Must have non-zero Jacobian
    requires !std::is_void_v<decltype(transform.jacobian(x))>;
};

// Accumulator for multivariate integration (extends univariate concept)
template <typename A, typename T>
concept MultivariateAccumulator = Accumulator<A> && requires(A acc) {
    typename A::value_type;
    requires std::same_as<typename A::value_type, T>;

    // Can accumulate weighted values
    { acc.add(std::declval<T>(), std::declval<T>()) } -> std::same_as<A&>;

    // Result extraction
    { acc.result() } -> std::convertible_to<T>;

    // Reset state
    { acc.reset() } -> std::same_as<void>;
};

// Parallel execution policy for multivariate integration
template <typename P>
concept ParallelPolicy = requires(P policy) {
    // Number of threads
    { policy.num_threads } -> std::convertible_to<size_t>;

    // Thread pool (if applicable)
    // { policy.thread_pool() } -> /* thread pool concept */;

    // Execution type (sequential, parallel, etc.)
    { policy.is_parallel() } -> std::convertible_to<bool>;
};

// Monte Carlo sampler
template <typename S, typename T, int Dim>
concept MonteCarloSampler = requires(S sampler, size_t n) {
    // Generate n samples in [0,1]^Dim
    { sampler.generate(n) } -> std::ranges::range;
    requires std::same_as<
        std::ranges::range_value_t<decltype(sampler.generate(n))>,
        std::array<T, Dim>
    >;

    // Type of sampler (random, quasi-random, stratified)
    { sampler.type() } -> std::convertible_to<std::string_view>;

    // Reset sampler state
    { sampler.reset() } -> std::same_as<void>;
};

// Quasi-random sequence (low-discrepancy)
template <typename Q, typename T, int Dim>
concept QuasiRandomSequence = MonteCarloSampler<Q, T, Dim> && requires(Q seq) {
    // Quasi-random sequences should be deterministic
    // Skip ahead in sequence
    { seq.skip(std::declval<size_t>()) } -> std::same_as<void>;
};

// Adaptive subdivision strategy
template <typename S, typename T, int Dim>
concept SubdivisionStrategy = requires(S strategy, IntegrationRegion<T, Dim> auto region) {
    // Subdivide region into multiple subregions
    { strategy.subdivide(region) } -> std::ranges::range;
    requires IntegrationRegion<
        std::ranges::range_value_t<decltype(strategy.subdivide(region))>,
        T,
        Dim
    >;

    // Priority/score for subdivision (higher = more important)
    { strategy.priority(
        std::declval<T>(),  // value
        std::declval<T>(),  // error
        region
    ) } -> std::convertible_to<T>;
};

// Error estimator for integration
template <typename E, typename T>
concept ErrorEstimator = requires(E estimator, T value1, T value2) {
    // Estimate error from two approximations
    { estimator.estimate(value1, value2) } -> std::convertible_to<T>;

    // Relative or absolute error
    { estimator.is_relative() } -> std::convertible_to<bool>;
};

} // namespace calckit::concepts
