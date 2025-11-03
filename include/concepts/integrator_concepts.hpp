#pragma once

#include <concepts>
#include <type_traits>
#include <functional>

namespace calckit::concepts {

// Fundamental numeric type concepts
template<typename T>
concept Arithmetic = std::is_arithmetic_v<T>;

template<typename T>
concept FloatingPoint = std::floating_point<T>;

template<typename T>
concept Field = requires(T a, T b) {
    { a + b } -> std::convertible_to<T>;
    { a - b } -> std::convertible_to<T>;
    { a * b } -> std::convertible_to<T>;
    { a / b } -> std::convertible_to<T>;
    { -a } -> std::convertible_to<T>;
    { T(0) };
    { T(1) };
};

template<typename T>
concept RealField = Field<T> && requires(T a) {
    { std::abs(a) } -> std::convertible_to<T>;
    { std::sqrt(a) } -> std::convertible_to<T>;
    { std::exp(a) } -> std::convertible_to<T>;
    { std::log(a) } -> std::convertible_to<T>;
    { std::sin(a) } -> std::convertible_to<T>;
    { std::cos(a) } -> std::convertible_to<T>;
};

// Function concepts
template<typename F, typename... Args>
concept Invocable = std::invocable<F, Args...>;

template<typename F, typename R, typename... Args>
concept InvocableR = std::invocable<F, Args...> &&
    std::convertible_to<std::invoke_result_t<F, Args...>, R>;

template<typename F, typename T>
concept UnivariateFunction = Invocable<F, T> &&
    requires { typename std::invoke_result_t<F, T>; };

template<typename F, typename T, typename R>
concept UnivariateFunctionR = UnivariateFunction<F, T> &&
    std::convertible_to<std::invoke_result_t<F, T>, R>;

template<typename F, typename T>
concept MultivariateFunction = requires(F f, T* begin, T* end) {
    { f(begin, end) } -> std::convertible_to<T>;
};

// Accumulator concepts
template<typename A, typename T>
concept Accumulator = requires(A acc, T value) {
    { acc += value } -> std::same_as<A&>;
    { acc() } -> std::convertible_to<T>;
    { A{} };
    { A{T{}} };
};

template<typename A, typename T>
concept CompensatedAccumulator = Accumulator<A, T> && requires(A acc) {
    { acc.error() } -> std::convertible_to<T>;
    { acc.correction() } -> std::convertible_to<T>;
};

// Integration result concept
template<typename R, typename T>
concept IntegrationResult = requires(R result) {
    { result.value() } -> std::convertible_to<T>;
    { result.error() } -> std::convertible_to<T>;
    { result.iterations() } -> std::convertible_to<std::size_t>;
    { static_cast<T>(result) } -> std::convertible_to<T>;
};

// Quadrature rule concept
template<typename Q, typename T>
concept QuadratureRule = requires(Q rule) {
    typename Q::value_type;
    typename Q::size_type;
    { rule.size() } -> std::convertible_to<std::size_t>;
    { rule.weight(std::size_t{}) } -> std::convertible_to<T>;
    { rule.abscissa(std::size_t{}) } -> std::convertible_to<T>;
};

// Adaptive strategy concept
template<typename S, typename T>
concept AdaptiveStrategy = requires(S strategy, T error, T tolerance, std::size_t depth) {
    { strategy.should_refine(error, tolerance, depth) } -> std::convertible_to<bool>;
    { strategy.subdivide(T{}, T{}) } -> std::convertible_to<std::pair<T, T>>;
    { strategy.combine(T{}, T{}) } -> std::convertible_to<T>;
};

// Integrator concept
template<typename I, typename T>
concept Integrator = requires {
    typename I::value_type;
    typename I::result_type;
    requires IntegrationResult<typename I::result_type, T>;
};

template<typename I, typename F, typename T>
concept UnivariateIntegrator = Integrator<I, T> &&
    requires(I integrator, F func, T a, T b, T tol) {
    { integrator(func, a, b) } -> std::same_as<typename I::result_type>;
    { integrator(func, a, b, tol) } -> std::same_as<typename I::result_type>;
};

// Sampler concept for Monte Carlo
template<typename S, typename T>
concept Sampler = requires(S sampler, std::size_t n) {
    typename S::value_type;
    { sampler.sample() } -> std::convertible_to<T>;
    { sampler.samples(n) } -> std::ranges::range;
};

// Transform concept for change of variables
template<typename T, typename U>
concept CoordinateTransform = requires(T transform, U x) {
    { transform(x) } -> std::convertible_to<U>;
    { transform.jacobian(x) } -> std::convertible_to<U>;
    { transform.inverse(x) } -> std::convertible_to<U>;
};

// Parallel execution policy concept
template<typename P>
concept ExecutionPolicy = requires(P policy) {
    { policy.threads() } -> std::convertible_to<std::size_t>;
    { policy.chunk_size() } -> std::convertible_to<std::size_t>;
    requires std::is_copy_constructible_v<P>;
};

} // namespace calckit::concepts