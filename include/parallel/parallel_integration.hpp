#pragma once

#include <algorithm>
#include <execution>
#include <thread>
#include <vector>
#include <numeric>
#include <future>
#include <functional>
#include <atomic>
#include <random>

#ifdef __AVX2__
#include <immintrin.h>
#endif

#include "../concepts/integrator_concepts.hpp"
#include "../core/integration_result.hpp"
#include "../integrators/univariate_integrator.hpp"

namespace calckit::parallel {

// Execution policy for parallel integration
struct execution_policy {
    std::size_t num_threads = std::thread::hardware_concurrency();
    std::size_t min_chunk_size = 1000;
    bool use_simd = true;

    constexpr std::size_t threads() const noexcept { return num_threads; }
    constexpr std::size_t chunk_size() const noexcept { return min_chunk_size; }
};

// SIMD-optimized evaluation for vectorizable functions
template<concepts::Field T>
class simd_evaluator {
public:
    static constexpr std::size_t simd_width =
#ifdef __AVX512F__
        8;  // AVX-512: 8 doubles
#elif defined(__AVX2__)
        4;  // AVX2: 4 doubles
#elif defined(__SSE2__)
        2;  // SSE2: 2 doubles
#else
        1;  // No SIMD
#endif

    template<typename F>
    static void evaluate_batch(F&& f, const T* x, T* y, std::size_t n) {
#ifdef __AVX2__
        if constexpr (std::is_same_v<T, double>) {
            evaluate_avx2(std::forward<F>(f), x, y, n);
        } else {
            evaluate_scalar(std::forward<F>(f), x, y, n);
        }
#else
        evaluate_scalar(std::forward<F>(f), x, y, n);
#endif
    }

private:
#ifdef __AVX2__
    template<typename F>
    static void evaluate_avx2(F&& f, const double* x, double* y, std::size_t n) {
        std::size_t simd_end = n - (n % 4);

        // Process 4 values at a time with AVX2
        for (std::size_t i = 0; i < simd_end; i += 4) {
            __m256d xvec = _mm256_loadu_pd(&x[i]);

            // Extract individual values and evaluate
            alignas(32) double vals[4];
            _mm256_store_pd(vals, xvec);

            alignas(32) double results[4];
            for (int j = 0; j < 4; ++j) {
                results[j] = f(vals[j]);
            }

            __m256d yvec = _mm256_load_pd(results);
            _mm256_storeu_pd(&y[i], yvec);
        }

        // Handle remaining elements
        for (std::size_t i = simd_end; i < n; ++i) {
            y[i] = f(x[i]);
        }
    }
#endif

    template<typename F>
    static void evaluate_scalar(F&& f, const T* x, T* y, std::size_t n) {
        for (std::size_t i = 0; i < n; ++i) {
            y[i] = f(x[i]);
        }
    }
};

// Parallel adaptive integration
template<concepts::Field T, typename Integrator>
class parallel_integrator {
public:
    using value_type = T;
    using result_type = integration_result<T>;
    using integrator_type = Integrator;

    parallel_integrator(Integrator base_integrator, execution_policy policy = {})
        : integrator_{std::move(base_integrator)}, policy_{policy} {}

    template<concepts::UnivariateFunction<T> F>
    result_type operator()(F&& f, T a, T b, T tol = default_tolerance()) const {
        if (b - a < T(policy_.min_chunk_size)) {
            // Small interval - use single-threaded integration
            return integrator_(std::forward<F>(f), a, b, tol);
        }

        // Divide interval into subintervals for parallel processing
        std::size_t n_chunks = std::min(policy_.num_threads,
            static_cast<std::size_t>((b - a) / policy_.min_chunk_size));

        if (n_chunks <= 1) {
            return integrator_(std::forward<F>(f), a, b, tol);
        }

        return parallel_adaptive(std::forward<F>(f), a, b, tol, n_chunks);
    }

    // Parallel integration over multiple disjoint intervals
    template<concepts::UnivariateFunction<T> F>
    result_type integrate_intervals(
        F&& f,
        const std::vector<std::pair<T, T>>& intervals,
        T tol = default_tolerance()
    ) const {
        std::vector<std::future<result_type>> futures;
        futures.reserve(intervals.size());

        for (const auto& [ai, bi] : intervals) {
            futures.emplace_back(
                std::async(std::launch::async,
                    [this, &f, ai, bi, tol, n = intervals.size()]() {
                        return integrator_(f, ai, bi, tol / n);
                    }
                )
            );
        }

        result_type total{};
        for (auto& future : futures) {
            total += future.get();
        }

        return total;
    }

private:
    Integrator integrator_;
    execution_policy policy_;

    static constexpr T default_tolerance() {
        return T(1000) * std::numeric_limits<T>::epsilon();
    }

    template<typename F>
    result_type parallel_adaptive(F&& f, T a, T b, T tol, std::size_t n_chunks) const {
        T h = (b - a) / T(n_chunks);
        std::vector<std::pair<T, T>> intervals;
        intervals.reserve(n_chunks);

        for (std::size_t i = 0; i < n_chunks; ++i) {
            T ai = a + i * h;
            T bi = (i == n_chunks - 1) ? b : ai + h;
            intervals.emplace_back(ai, bi);
        }

        // Parallel integration with work stealing
        std::vector<result_type> results(n_chunks);
        std::atomic<std::size_t> next_interval{0};

        auto worker = [&]() {
            while (true) {
                std::size_t i = next_interval.fetch_add(1);
                if (i >= n_chunks) break;

                const auto& [ai, bi] = intervals[i];
                results[i] = integrator_(f, ai, bi, tol / n_chunks);
            }
        };

        std::vector<std::thread> threads;
        threads.reserve(policy_.num_threads);

        for (std::size_t i = 0; i < policy_.num_threads; ++i) {
            threads.emplace_back(worker);
        }

        for (auto& thread : threads) {
            thread.join();
        }

        // Combine results
        return std::reduce(results.begin(), results.end(), result_type{});
    }
};

// Parallel Monte Carlo integration
template<concepts::Field T>
class parallel_monte_carlo {
public:
    using value_type = T;
    using result_type = integration_result<T>;

    parallel_monte_carlo(std::size_t samples = 1000000, execution_policy policy = {})
        : n_samples_{samples}, policy_{policy} {}

    template<concepts::UnivariateFunction<T> F>
    result_type operator()(F&& f, T a, T b) const {
        const std::size_t samples_per_thread = n_samples_ / policy_.num_threads;
        std::vector<std::future<std::pair<T, T>>> futures;

        for (std::size_t t = 0; t < policy_.num_threads; ++t) {
            std::size_t thread_samples = (t == policy_.num_threads - 1)
                ? n_samples_ - t * samples_per_thread
                : samples_per_thread;

            futures.emplace_back(
                std::async(std::launch::async,
                    [this, &f, a, b, thread_samples, seed = t]() {
                        return monte_carlo_worker(f, a, b, thread_samples, seed);
                    }
                )
            );
        }

        T sum = T(0);
        T sum_sq = T(0);

        for (auto& future : futures) {
            auto [partial_sum, partial_sum_sq] = future.get();
            sum += partial_sum;
            sum_sq += partial_sum_sq;
        }

        T mean = sum / T(n_samples_);
        T variance = sum_sq / T(n_samples_) - mean * mean;
        T width = b - a;
        T integral = width * mean;
        T error = width * std::sqrt(variance / T(n_samples_));

        return {integral, error, 1, n_samples_};
    }

    // Multivariate Monte Carlo
    template<typename F>
    result_type integrate_multivariate(
        F&& f,
        const std::vector<std::pair<T, T>>& bounds
    ) const {
        const std::size_t dim = bounds.size();
        const std::size_t samples_per_thread = n_samples_ / policy_.num_threads;
        std::vector<std::future<std::pair<T, T>>> futures;

        for (std::size_t t = 0; t < policy_.num_threads; ++t) {
            std::size_t thread_samples = (t == policy_.num_threads - 1)
                ? n_samples_ - t * samples_per_thread
                : samples_per_thread;

            futures.emplace_back(
                std::async(std::launch::async,
                    [this, &f, &bounds, thread_samples, seed = t]() {
                        return monte_carlo_multivariate_worker(
                            f, bounds, thread_samples, seed
                        );
                    }
                )
            );
        }

        T sum = T(0);
        T sum_sq = T(0);

        for (auto& future : futures) {
            auto [partial_sum, partial_sum_sq] = future.get();
            sum += partial_sum;
            sum_sq += partial_sum_sq;
        }

        T volume = std::accumulate(bounds.begin(), bounds.end(), T(1),
            [](T v, const auto& bound) { return v * (bound.second - bound.first); }
        );

        T mean = sum / T(n_samples_);
        T variance = sum_sq / T(n_samples_) - mean * mean;
        T integral = volume * mean;
        T error = volume * std::sqrt(variance / T(n_samples_));

        return {integral, error, 1, n_samples_};
    }

private:
    std::size_t n_samples_;
    execution_policy policy_;

    template<typename F>
    std::pair<T, T> monte_carlo_worker(
        F&& f, T a, T b, std::size_t samples, std::size_t seed
    ) const {
        std::mt19937_64 rng(seed);
        std::uniform_real_distribution<T> dist(a, b);

        T sum = T(0);
        T sum_sq = T(0);

        for (std::size_t i = 0; i < samples; ++i) {
            T x = dist(rng);
            T y = f(x);
            sum += y;
            sum_sq += y * y;
        }

        return {sum, sum_sq};
    }

    template<typename F>
    std::pair<T, T> monte_carlo_multivariate_worker(
        F&& f,
        const std::vector<std::pair<T, T>>& bounds,
        std::size_t samples,
        std::size_t seed
    ) const {
        std::mt19937_64 rng(seed);
        std::vector<std::uniform_real_distribution<T>> dists;

        for (const auto& [a, b] : bounds) {
            dists.emplace_back(a, b);
        }

        std::vector<T> point(bounds.size());
        T sum = T(0);
        T sum_sq = T(0);

        for (std::size_t i = 0; i < samples; ++i) {
            for (std::size_t d = 0; d < bounds.size(); ++d) {
                point[d] = dists[d](rng);
            }

            T y = f(point.data(), point.data() + point.size());
            sum += y;
            sum_sq += y * y;
        }

        return {sum, sum_sq};
    }
};

// Factory functions
template<concepts::Field T>
auto make_parallel_integrator(execution_policy policy = {}) {
    auto base = make_adaptive_integrator<T>();
    return parallel_integrator<T, decltype(base)>{std::move(base), policy};
}

template<concepts::Field T>
auto make_parallel_monte_carlo(std::size_t samples = 1000000, execution_policy policy = {}) {
    return parallel_monte_carlo<T>{samples, policy};
}

} // namespace calckit::parallel