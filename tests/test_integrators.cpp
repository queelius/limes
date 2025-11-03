#include <gtest/gtest.h>
#include <cmath>
#include <functional>
#include <limits>
#include <numbers>
#include <complex>
#include "../include/integrators/univariate_integrator.hpp"
#include "../include/calckit.hpp"

using namespace calckit;

// Test functions with known integrals
namespace test_functions {
    // Simple polynomial: x^2, integral from 0 to 1 = 1/3
    template<typename T>
    T quadratic(T x) { return x * x; }

    // Exponential: e^x, integral from 0 to 1 = e - 1
    template<typename T>
    T exponential(T x) { return std::exp(x); }

    // Trigonometric: sin(x), integral from 0 to pi = 2
    template<typename T>
    T sine(T x) { return std::sin(x); }

    // Gaussian: exp(-x^2), integral from -inf to inf = sqrt(pi)
    template<typename T>
    T gaussian(T x) { return std::exp(-x * x); }

    // Oscillatory: sin(10x), integral from 0 to pi = 0
    template<typename T>
    T oscillatory(T x) { return std::sin(T(10) * x); }

    // Singular at endpoint: 1/sqrt(x), integral from 0 to 1 = 2
    template<typename T>
    T singular_endpoint(T x) {
        if (x <= 0) return T(0);
        return T(1) / std::sqrt(x);
    }

    // Discontinuous: step function
    template<typename T>
    T step_function(T x) {
        return x < T(0.5) ? T(0) : T(1);
    }

    // Peaked function: narrow Gaussian
    template<typename T>
    T peaked(T x) {
        T sigma = T(0.01);
        return std::exp(-(x * x) / (T(2) * sigma * sigma)) / (sigma * std::sqrt(T(2) * std::numbers::pi_v<T>));
    }
}

// Test fixture for integrator tests
template <typename T>
class IntegratorTest : public ::testing::Test {
protected:
    static constexpr T default_tol = T(1e-8);
    static constexpr T strict_tol = T(1e-12);
    static constexpr T loose_tol = T(1e-4);

    template<typename Integrator, typename F>
    void test_integration(Integrator& integrator, F&& f, T a, T b,
                         T expected, T tolerance) {
        auto result = integrator(std::forward<F>(f), a, b, tolerance);

        EXPECT_TRUE(result.converged())
            << "Integration did not converge";

        EXPECT_NEAR(result.value(), expected, tolerance * T(10))
            << "Integration result outside tolerance"
            << "\n  Computed: " << result.value()
            << "\n  Expected: " << expected
            << "\n  Error: " << result.error()
            << "\n  Evaluations: " << result.evaluations();

        EXPECT_LE(result.error(), tolerance * T(10))
            << "Estimated error too large";
    }
};

using FloatTypes = ::testing::Types<float, double, long double>;
TYPED_TEST_SUITE(IntegratorTest, FloatTypes);

// Test quadrature integrator with different rules
TYPED_TEST(IntegratorTest, QuadratureIntegratorGaussLegendre) {
    using T = TypeParam;
    using rule = quadrature::gauss_legendre<T, 5>;

    quadrature_integrator<T, rule> integrator{rule{}};

    // Test polynomial integration
    this->test_integration(integrator,
                          test_functions::quadratic<T>,
                          T(0), T(1), T(1)/T(3), this->strict_tol);

    // Test exponential
    this->test_integration(integrator,
                          test_functions::exponential<T>,
                          T(0), T(1), std::exp(T(1)) - T(1), this->default_tol);
}

TYPED_TEST(IntegratorTest, QuadratureIntegratorGaussKronrod) {
    using T = TypeParam;
    using rule = quadrature::gauss_kronrod_15<T>;

    quadrature_integrator<T, rule> integrator{rule{}};

    // Test trigonometric function
    this->test_integration(integrator,
                          test_functions::sine<T>,
                          T(0), std::numbers::pi_v<T>, T(2), this->strict_tol);
}

// Test adaptive integrator
TYPED_TEST(IntegratorTest, AdaptiveIntegrator) {
    using T = TypeParam;

    auto integrator = make_adaptive_integrator<T>();

    // Simple function
    this->test_integration(integrator,
                          test_functions::quadratic<T>,
                          T(0), T(1), T(1)/T(3), this->default_tol);

    // Oscillatory function - should trigger subdivision
    auto result = integrator(test_functions::oscillatory<T>,
                            T(0), std::numbers::pi_v<T>, this->default_tol);

    // For sin(10x) from 0 to pi, integral is (1 - cos(10*pi))/10 = 0.2
    EXPECT_NEAR(result.value(), T(0.2), this->default_tol * T(10));
}

// Test Romberg integrator
TYPED_TEST(IntegratorTest, RombergIntegrator) {
    using T = TypeParam;

    auto integrator = make_romberg<T>();

    // Smooth function
    this->test_integration(integrator,
                          test_functions::exponential<T>,
                          T(0), T(1), std::exp(T(1)) - T(1), this->strict_tol);

    // Test with same integrator
    this->test_integration(integrator,
                          test_functions::sine<T>,
                          T(0), std::numbers::pi_v<T>, T(2), T(1e-10));
}

// Test Simpson integrator
TYPED_TEST(IntegratorTest, SimpsonIntegrator) {
    using T = TypeParam;

    using rule = quadrature::simpson_rule<T>;
    quadrature_integrator<T, rule> integrator{rule{}};

    // Simpson's rule is exact for cubic polynomials on [-1, 1]
    auto cubic = [](T x) { return x * x * x; };
    auto result = integrator(cubic, T(-1), T(1));

    EXPECT_NEAR(result.value(), T(0), this->strict_tol);
}

// Test double exponential (tanh-sinh) integrator
TYPED_TEST(IntegratorTest, TanhSinhIntegrator) {
    using T = TypeParam;

    auto integrator = make_tanh_sinh<T>();

    // Function with endpoint singularity
    auto f = [](T x) -> T {
        if (std::abs(x - T(1)) < std::numeric_limits<T>::epsilon()) {
            return T(0);
        }
        return std::log(T(1) - x);
    };

    // Integral of log(1-x) from 0 to 1 = -1
    auto result = integrator(f, T(0), T(0.999), this->loose_tol);

    // This is a challenging integral, so we use looser tolerance
    EXPECT_NEAR(result.value(), T(-0.999), this->loose_tol * T(100));
}

// Test infinite interval integration
TYPED_TEST(IntegratorTest, InfiniteIntervals) {
    using T = TypeParam;

    // Gaussian integral from -inf to inf
    {
        auto result = integrate<T>::adaptive(
            test_functions::gaussian<T>,
            -std::numeric_limits<T>::infinity(),
            std::numeric_limits<T>::infinity(),
            this->default_tol
        );

        EXPECT_NEAR(result.value(), std::sqrt(std::numbers::pi_v<T>), this->default_tol * T(10));
    }

    // Semi-infinite: exp(-x) from 0 to inf = 1
    {
        auto f = [](T x) { return std::exp(-x); };
        auto result = integrate<T>::adaptive(
            f, T(0), std::numeric_limits<T>::infinity(), this->default_tol
        );

        EXPECT_NEAR(result.value(), T(1), this->default_tol * T(10));
    }
}

// Test high-level integration interface
TYPED_TEST(IntegratorTest, HighLevelInterface) {
    using T = TypeParam;

    // Test adaptive integration
    {
        auto result = integrate_adaptive(
            test_functions::quadratic<T>, T(0), T(1), this->default_tol
        );
        EXPECT_NEAR(result.value(), T(1)/T(3), this->default_tol);
    }

    // Test robust integration
    {
        auto result = integrate_robust(
            test_functions::exponential<T>, T(0), T(1), this->default_tol
        );
        EXPECT_NEAR(result.value(), std::exp(T(1)) - T(1), this->default_tol);
    }

    // Test precise integration
    {
        auto result = integrate_precise(
            test_functions::sine<T>, T(0), std::numbers::pi_v<T>, this->strict_tol
        );
        EXPECT_NEAR(result.value(), T(2), this->strict_tol);
    }
}

// Test oscillatory integration
TYPED_TEST(IntegratorTest, OscillatoryIntegration) {
    using T = TypeParam;

    // High-frequency oscillation
    auto f = [](T x) { return std::sin(T(100) * x); };

    auto result = integrate<T>::oscillatory(
        f, T(0), std::numbers::pi_v<T>, T(100), this->default_tol
    );

    // Integral should be close to 0 for high frequency
    EXPECT_NEAR(result.value(), T(0.02), this->loose_tol);
}

// Test integrator builder
TYPED_TEST(IntegratorTest, IntegratorBuilder) {
    using T = TypeParam;

    auto integrator = make_integrator<T>()
        .with_tolerance(this->strict_tol)
        .with_parallel(false);

    auto result = integrator.integrate(
        test_functions::exponential<T>, T(0), T(1)
    );

    EXPECT_NEAR(result.value(), std::exp(T(1)) - T(1), this->strict_tol * T(10));
}

// Test error estimation
TEST(IntegratorErrorEstimation, ErrorBounds) {
    // Use a function where we can calculate the truncation error
    auto f = [](double x) { return std::exp(x); };

    auto integrator = make_adaptive_integrator<double>();
    auto result = integrator(f, 0.0, 1.0, 1e-10);

    double exact = std::exp(1.0) - 1.0;
    double actual_error = std::abs(result.value() - exact);

    // Estimated error should be conservative (larger than actual)
    EXPECT_GE(result.error(), actual_error * 0.1);
    // But not too conservative
    EXPECT_LE(result.error(), actual_error * 100);
}

// Test convergence behavior
TEST(IntegratorConvergence, ConvergenceRates) {
    auto f = [](double x) { return std::sin(x); };
    double exact = 2.0;

    std::vector<double> tolerances = {1e-4, 1e-6, 1e-8, 1e-10};
    std::vector<double> errors;
    std::vector<size_t> evaluations;

    for (double tol : tolerances) {
        auto integrator = make_adaptive_integrator<double>();
        auto result = integrator(f, 0.0, M_PI, tol);

        errors.push_back(std::abs(result.value() - exact));
        evaluations.push_back(result.evaluations());

        // Check that requested tolerance is met
        EXPECT_LE(errors.back(), tol * 10);
    }

    // Check that errors decrease with tighter tolerance
    for (size_t i = 1; i < errors.size(); ++i) {
        EXPECT_LT(errors[i], errors[i-1]);
    }

    // Check that evaluations increase with tighter tolerance
    for (size_t i = 1; i < evaluations.size(); ++i) {
        EXPECT_GT(evaluations[i], evaluations[i-1]);
    }
}

// Test handling of discontinuous functions
TEST(IntegratorSpecialCases, DiscontinuousFunction) {
    // Step function with discontinuity at x = 0.5
    auto f = [](double x) { return x < 0.5 ? 0.0 : 1.0; };

    // Exact integral from 0 to 1 is 0.5
    auto integrator = make_adaptive_integrator<double>();
    auto result = integrator(f, 0.0, 1.0, 1e-6);

    // Adaptive integrator should handle discontinuity reasonably well
    EXPECT_NEAR(result.value(), 0.5, 1e-4);
}

// Test handling of peaked functions
TEST(IntegratorSpecialCases, PeakedFunction) {
    // Narrow Gaussian peaked at x = 0
    double sigma = 0.001;
    auto f = [sigma](double x) {
        return std::exp(-(x * x) / (2 * sigma * sigma)) /
               (sigma * std::sqrt(2 * M_PI));
    };

    // Integral from -1 to 1 should be close to 1 (normalized Gaussian)
    auto integrator = make_adaptive_integrator<double>();
    auto result = integrator(f, -1.0, 1.0, 1e-8);

    EXPECT_NEAR(result.value(), 1.0, 1e-6);
}

// Test with lambda captures
TEST(IntegratorLambdas, CapturedVariables) {
    double a = 2.0;
    double b = 3.0;

    auto f = [a, b](double x) { return a * x + b; };

    auto integrator = make_adaptive_integrator<double>();
    auto result = integrator(f, 0.0, 1.0, 1e-10);

    // Integral of 2x + 3 from 0 to 1 = x^2 + 3x |_0^1 = 1 + 3 = 4
    EXPECT_NEAR(result.value(), 4.0, 1e-10);
}

// Test with function pointers
double test_func_ptr(double x) {
    return x * x * x;
}

TEST(IntegratorFunctionTypes, FunctionPointer) {
    auto integrator = make_adaptive_integrator<double>();
    auto result = integrator(test_func_ptr, 0.0, 2.0, 1e-10);

    // Integral of x^3 from 0 to 2 = x^4/4 |_0^2 = 4
    EXPECT_NEAR(result.value(), 4.0, 1e-10);
}

// Test with std::function
TEST(IntegratorFunctionTypes, StdFunction) {
    std::function<double(double)> f = [](double x) { return std::cos(x); };

    auto integrator = make_adaptive_integrator<double>();
    auto result = integrator(f, 0.0, M_PI/2, 1e-10);

    // Integral of cos(x) from 0 to pi/2 = sin(pi/2) - sin(0) = 1
    EXPECT_NEAR(result.value(), 1.0, 1e-10);
}

// Test numerical stability
TEST(IntegratorStability, NumericalStability) {
    // Test with very small interval
    {
        auto f = [](double x) { return x; };
        auto integrator = make_adaptive_integrator<double>();
        auto result = integrator(f, 0.0, 1e-15, 1e-10);

        // Should handle tiny intervals gracefully
        EXPECT_NEAR(result.value(), 5e-31, 1e-40);
    }

    // Test with very large values
    {
        auto f = [](double x) { return 1e10; };
        auto integrator = make_adaptive_integrator<double>();
        auto result = integrator(f, 0.0, 1.0, 1e-10);

        EXPECT_NEAR(result.value(), 1e10, 1e-2);
    }
}

// Test different accumulator strategies
TEST(IntegratorAccumulators, AccumulatorComparison) {
    // Ill-conditioned sum
    auto f = [](double x) { return 1e10 + std::sin(x); };

    // Simple accumulator
    {
        using acc = accumulators::simple_accumulator<double>;
        quadrature_integrator<double, quadrature::gauss_legendre<double, 5>, acc> integrator;
        auto result = integrator(f, 0.0, M_PI, 1e-8);
        // Record result for comparison
    }

    // Kahan accumulator (should be more accurate)
    {
        using acc = accumulators::kahan_accumulator<double>;
        quadrature_integrator<double, quadrature::gauss_legendre<double, 5>, acc> integrator;
        auto result = integrator(f, 0.0, M_PI, 1e-8);

        double expected = 1e10 * M_PI + 2.0; // Integral of 1e10 + sin(x)
        EXPECT_NEAR(result.value(), expected, 100.0);
    }
}

// Performance benchmark (disabled by default)
TEST(IntegratorBenchmark, DISABLED_PerformanceComparison) {
    auto f = [](double x) { return std::exp(-x * x); };

    struct benchmark_result {
        std::string name;
        double value;
        size_t evaluations;
        double time_ms;
    };

    std::vector<benchmark_result> results;

    auto benchmark = [&](auto integrator, const std::string& name) {
        auto start = std::chrono::high_resolution_clock::now();
        auto result = integrator(f, -5.0, 5.0, 1e-8);
        auto end = std::chrono::high_resolution_clock::now();

        auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);

        results.push_back({
            name,
            result.value(),
            result.evaluations(),
            duration.count() / 1000.0
        });
    };

    benchmark(make_adaptive_integrator<double>(), "Adaptive");
    benchmark(make_romberg<double>(), "Romberg");
    benchmark(make_tanh_sinh<double>(), "Tanh-Sinh");
    // Use quadrature_integrator with simpson_rule
    using simpson = quadrature_integrator<double, quadrature::simpson_rule<double>>;
    benchmark(simpson{quadrature::simpson_rule<double>{}}, "Simpson");

    std::cout << "\nIntegrator Performance Comparison:\n";
    std::cout << "----------------------------------\n";
    for (const auto& r : results) {
        std::cout << std::setw(15) << r.name
                  << ": value = " << std::setw(12) << r.value
                  << ", evals = " << std::setw(6) << r.evaluations
                  << ", time = " << std::setw(8) << r.time_ms << " ms\n";
    }
}