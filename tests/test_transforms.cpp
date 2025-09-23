#include <gtest/gtest.h>
#include <cmath>
#include <limits>
#include <numbers>
#include <random>
#include "../include/transforms/coordinate_transforms.hpp"
#include "../include/algebraic_integrators.hpp"

using namespace algebraic_integrators::transforms;

// Test fixture for transform tests
template <typename T>
class TransformTest : public ::testing::Test {
protected:
    static constexpr T eps = std::numeric_limits<T>::epsilon();
    static constexpr T tol = eps * T(100);

    // Helper to test transform invertibility
    void test_invertibility(const auto& transform, T x, T tolerance = tol) {
        T y = transform.forward(x);
        T x_recovered = transform.inverse(y);

        EXPECT_NEAR(x_recovered, x, tolerance)
            << "Transform not invertible at x = " << x
            << ", forward(x) = " << y
            << ", inverse(forward(x)) = " << x_recovered;
    }

    // Helper to test Jacobian using finite differences
    void test_jacobian_numerical(const auto& transform, T x, T h = T(1e-8)) {
        T jacobian_analytical = transform.jacobian(x);

        // Numerical Jacobian using central difference
        T f_plus = transform.forward(x + h);
        T f_minus = transform.forward(x - h);
        T jacobian_numerical = (f_plus - f_minus) / (T(2) * h);

        T rel_error = std::abs(jacobian_analytical - jacobian_numerical) /
                     std::max(std::abs(jacobian_analytical), T(1));

        EXPECT_LT(rel_error, T(1e-6))
            << "Jacobian mismatch at x = " << x
            << ", analytical = " << jacobian_analytical
            << ", numerical = " << jacobian_numerical;
    }
};

using FloatTypes = ::testing::Types<float, double, long double>;
TYPED_TEST_SUITE(TransformTest, FloatTypes);

// Test interval mapping
TYPED_TEST(TransformTest, IntervalMap) {
    using T = TypeParam;

    // Map [0, 10] to [-1, 1]
    interval_map<T> transform(T(0), T(10));

    // Test endpoints
    EXPECT_NEAR(transform.forward(T(-1)), T(0), this->tol);
    EXPECT_NEAR(transform.forward(T(1)), T(10), this->tol);

    // Test midpoint
    EXPECT_NEAR(transform.forward(T(0)), T(5), this->tol);

    // Test inverse
    EXPECT_NEAR(transform.inverse(T(0)), T(-1), this->tol);
    EXPECT_NEAR(transform.inverse(T(10)), T(1), this->tol);
    EXPECT_NEAR(transform.inverse(T(5)), T(0), this->tol);

    // Test Jacobian (should be constant = 5)
    EXPECT_NEAR(transform.jacobian(T(-0.5)), T(5), this->tol);
    EXPECT_NEAR(transform.jacobian(T(0)), T(5), this->tol);
    EXPECT_NEAR(transform.jacobian(T(0.5)), T(5), this->tol);

    // Test invertibility at random points
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<T> dist(T(-1), T(1));

    for (int i = 0; i < 100; ++i) {
        T x = dist(gen);
        this->test_invertibility(transform, x);
    }
}

TYPED_TEST(TransformTest, IntervalMapNegative) {
    using T = TypeParam;

    // Map [-5, -2] to [-1, 1]
    interval_map<T> transform(T(-5), T(-2));

    EXPECT_NEAR(transform.forward(T(-1)), T(-5), this->tol);
    EXPECT_NEAR(transform.forward(T(1)), T(-2), this->tol);

    // Jacobian should be 1.5
    EXPECT_NEAR(transform.jacobian(T(0)), T(1.5), this->tol);

    this->test_invertibility(transform, T(0));
}

TYPED_TEST(TransformTest, IntervalMapLargeRange) {
    using T = TypeParam;

    // Map very large interval [0, 1000000] to [-1, 1]
    interval_map<T> transform(T(0), T(1000000));

    // Test that it correctly maps
    EXPECT_NEAR(transform.forward(T(-1)), T(0), this->tol * T(1000000));
    EXPECT_NEAR(transform.forward(T(1)), T(1000000), this->tol * T(1000000));
    EXPECT_NEAR(transform.forward(T(0)), T(500000), this->tol * T(1000000));

    // Jacobian should be 500000
    EXPECT_NEAR(transform.jacobian(T(0)), T(500000), this->tol * T(500000));

    // Test invertibility
    this->test_invertibility(transform, T(0), this->tol * T(1000000));
    this->test_invertibility(transform, T(0.5), this->tol * T(1000000));
    this->test_invertibility(transform, T(-0.5), this->tol * T(1000000));
}

TYPED_TEST(TransformTest, IntervalMapSmallRange) {
    using T = TypeParam;

    // Map very small interval [0, 0.001] to [-1, 1]
    interval_map<T> transform(T(0), T(0.001));

    // Test that it correctly maps
    EXPECT_NEAR(transform.forward(T(-1)), T(0), this->tol);
    EXPECT_NEAR(transform.forward(T(1)), T(0.001), this->tol);
    EXPECT_NEAR(transform.forward(T(0)), T(0.0005), this->tol);

    // Jacobian should be 0.0005
    EXPECT_NEAR(transform.jacobian(T(0)), T(0.0005), this->tol);

    // Test invertibility
    this->test_invertibility(transform, T(0));
}

// Test using transforms with integration
TEST(TransformIntegration, IntervalTransform) {
    // Integrate x^2 from 2 to 5 using transform to [-1,1]
    auto f = [](double x) { return x * x; };

    interval_map<double> transform(2.0, 5.0);

    auto transformed_f = [&](double t) {
        double x = transform.forward(t);
        return f(x) * transform.jacobian(t);
    };

    // Integrate over [-1, 1] using the high-level API
    auto result = algebraic_integrators::integrate_adaptive(transformed_f, -1.0, 1.0, 1e-10);

    // Exact integral: x^3/3 from 2 to 5 = 125/3 - 8/3 = 117/3 = 39
    double exact = 39.0;
    EXPECT_NEAR(result.value(), exact, 1e-8);
}

// Test transform composition manually
TEST(TransformComposition, ManualComposition) {
    // Compose two interval maps manually
    interval_map<double> map1(0.0, 10.0);  // Maps [-1,1] to [0,10]
    interval_map<double> map2(0.0, 100.0); // Maps [-1,1] to [0,100]

    // Compose: [-1,1] -> [0,10] -> map to [0,100] space
    double x = 0.5;
    double y1 = map1.forward(x);  // Should give 7.5

    // To use map2, we need to first convert y1 to [-1,1] space
    // This would require mapping [0,10] to [-1,1] then applying map2
    double t = (y1 - 5.0) / 5.0;  // Map [0,10] to [-1,1]
    double y2 = map2.forward(t);  // Map to [0,100]

    EXPECT_GT(y1, 7.0);
    EXPECT_LT(y1, 8.0);
    EXPECT_GT(y2, 70.0);
    EXPECT_LT(y2, 80.0);
}

// Test edge cases
TEST(TransformEdgeCases, BoundaryBehavior) {
    // Test interval map with zero-width interval
    {
        interval_map<double> transform(5.0, 5.0);
        EXPECT_EQ(transform.forward(0.0), 5.0);
        // Jacobian for zero-width interval should be 0 or inf
        double jac = transform.jacobian(0.0);
        EXPECT_TRUE(jac == 0.0 || std::isinf(jac));
    }

    // Test interval map with very small interval
    {
        double eps = std::numeric_limits<double>::epsilon();
        interval_map<double> transform(0.0, eps);

        // Should still work but with very small Jacobian
        double jac = transform.jacobian(0.0);
        EXPECT_LT(jac, 1e-10);
    }
}

// Test numerical stability
TEST(TransformStability, NumericalStability) {
    // Test with very large intervals
    {
        double large = 1e15;
        interval_map<double> transform(-large, large);

        // Should handle large values
        double y = transform.forward(0.0);
        EXPECT_FALSE(std::isnan(y));
        EXPECT_FALSE(std::isinf(y));
        EXPECT_NEAR(y, 0.0, large * 1e-10);
    }

    // Test with very small intervals
    {
        double tiny = std::numeric_limits<double>::min() * 100;
        interval_map<double> transform(0.0, tiny);

        // Should handle tiny values
        double y = transform.forward(0.0);
        EXPECT_FALSE(std::isnan(y));
        EXPECT_NEAR(y, tiny / 2.0, tiny);
    }
}

// Test transform usage patterns
TEST(TransformUsage, CommonPatterns) {
    // Pattern 1: Normalize domain to [-1, 1] for quadrature
    {
        double a = 3.0, b = 7.0;
        interval_map<double> normalize(a, b);

        // Any quadrature on [-1, 1] can now be applied to [a, b]
        auto f = [](double x) { return x; }; // Linear function

        double sum = 0.0;
        // Simple 2-point quadrature
        double x1 = -1.0 / std::sqrt(3.0);
        double x2 = 1.0 / std::sqrt(3.0);
        double w = 1.0;

        sum += w * f(normalize.forward(x1)) * normalize.jacobian(x1);
        sum += w * f(normalize.forward(x2)) * normalize.jacobian(x2);

        // Exact integral of x from 3 to 7 is (49-9)/2 = 20
        EXPECT_NEAR(sum, 20.0, 1e-10);
    }

    // Pattern 2: Scale results
    {
        interval_map<double> scale(0.0, 2.0 * M_PI);

        // Can be used to integrate periodic functions
        auto f = [](double x) { return std::sin(x); };

        // Transform from [-1, 1] to [0, 2π]
        double t = 0.0; // Middle of [-1, 1]
        double x = scale.forward(t);
        EXPECT_NEAR(x, M_PI, 1e-10);
        EXPECT_NEAR(f(x), 0.0, 1e-10); // sin(π) = 0
    }
}

// Performance test (disabled by default)
TEST(TransformBenchmark, DISABLED_Performance) {
    const size_t n = 1000000;
    std::vector<double> x_values(n);

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<double> dist(-1.0, 1.0);

    for (auto& x : x_values) {
        x = dist(gen);
    }

    interval_map<double> transform(0.0, 10.0);

    auto start = std::chrono::high_resolution_clock::now();

    double sum = 0.0;
    for (double x : x_values) {
        sum += transform.forward(x) * transform.jacobian(x);
    }

    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);

    std::cout << "Interval Map performance: "
              << duration.count() << " us for " << n << " evaluations"
              << ", sum = " << sum << "\n";
}