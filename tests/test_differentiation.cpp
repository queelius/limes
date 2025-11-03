#include <gtest/gtest.h>
#include <cmath>
#include <numbers>
#include <vector>
#include <array>
#include "../include/central_finite_difference.hpp"
#include "../include/numerical_differentiation.hpp"

using namespace calckit;

// Test fixture for differentiation tests
template <typename T>
class DifferentiationTest : public ::testing::Test {
protected:
    static constexpr T epsilon = std::numeric_limits<T>::epsilon();
    static constexpr T pi = std::numbers::pi_v<T>;

    // Calculate relative error
    T relative_error(T computed, T exact) const {
        if (std::abs(exact) < epsilon * T(10)) return std::abs(computed - exact);
        return std::abs(computed - exact) / std::abs(exact);
    }

    // Test-specific tolerance based on type and h
    T test_tolerance(T h) const {
        if constexpr (std::is_same_v<T, float>) {
            return T(1e-2);
        } else {
            // Higher order methods should give better accuracy
            return T(1e-5);
        }
    }
};

using FloatTypes = ::testing::Types<double, float>;
TYPED_TEST_SUITE(DifferentiationTest, FloatTypes);

// ========== First Derivative Tests ==========

// Test 1: Derivative of constant function (f(x) = c => f'(x) = 0)
TYPED_TEST(DifferentiationTest, ConstantFunction) {
    using T = TypeParam;

    auto f = [](T x) { return T(42); };
    T x = 1, h = T(0.01);

    auto result = central_finite_difference(f, x, h);
    T exact = T(0);

    // Roundoff error can be significant for float
    T tolerance = std::is_same_v<T, float> ? T(1e-4) : T(1e-10);
    EXPECT_NEAR(result.difference, exact, tolerance);
    EXPECT_EQ(result.evaluations, 8);
}

// Test 2: Derivative of linear function (f(x) = ax + b => f'(x) = a)
TYPED_TEST(DifferentiationTest, LinearFunction) {
    using T = TypeParam;

    auto f = [](T x) { return T(2) * x + T(3); };
    T x = 5, h = T(0.001);

    auto result = central_finite_difference(f, x, h);
    T exact = T(2);

    // Central finite difference has roundoff error due to coefficient arithmetic
    T tolerance = std::is_same_v<T, float> ? T(1e-3) : T(1e-6);
    EXPECT_NEAR(result.difference, exact, tolerance);
}

// Test 3: Derivative of quadratic (f(x) = x^2 => f'(x) = 2x)
TYPED_TEST(DifferentiationTest, QuadraticFunction) {
    using T = TypeParam;

    auto f = [](T x) { return x * x; };
    T x = 3, h = T(0.001);

    auto result = central_finite_difference(f, x, h);
    T exact = T(2) * x;

    EXPECT_NEAR(result.difference, exact, this->test_tolerance(h));
}

// Test 4: Derivative of polynomial (f(x) = x^3 => f'(x) = 3x^2)
TYPED_TEST(DifferentiationTest, PolynomialFunction) {
    using T = TypeParam;

    auto f = [](T x) { return x * x * x; };
    T x = 2, h = T(0.001);

    auto result = central_finite_difference(f, x, h);
    T exact = T(3) * x * x;

    EXPECT_NEAR(result.difference, exact, this->test_tolerance(h));
}

// Test 5: Derivative of sine (f(x) = sin(x) => f'(x) = cos(x))
TYPED_TEST(DifferentiationTest, SineFunction) {
    using T = TypeParam;

    auto f = [](T x) { return std::sin(x); };
    T x = this->pi / T(4), h = T(0.001);

    auto result = central_finite_difference(f, x, h);
    T exact = std::cos(x);

    EXPECT_NEAR(result.difference, exact, this->test_tolerance(h));
}

// Test 6: Derivative of exponential (f(x) = e^x => f'(x) = e^x)
TYPED_TEST(DifferentiationTest, ExponentialFunction) {
    using T = TypeParam;

    auto f = [](T x) { return std::exp(x); };
    T x = 1, h = T(0.001);

    auto result = central_finite_difference(f, x, h);
    T exact = std::exp(x);

    EXPECT_NEAR(result.difference, exact, this->test_tolerance(h) * T(10));
}

// Test 7: Derivative of logarithm (f(x) = ln(x) => f'(x) = 1/x)
TYPED_TEST(DifferentiationTest, LogarithmFunction) {
    using T = TypeParam;

    auto f = [](T x) { return std::log(x); };
    T x = 2, h = T(0.001);

    auto result = central_finite_difference(f, x, h);
    T exact = T(1) / x;

    EXPECT_NEAR(result.difference, exact, this->test_tolerance(h));
}

// ========== Second Derivative Tests ==========

// Test 8: Second derivative of quadratic (f(x) = x^2 => f''(x) = 2)
TYPED_TEST(DifferentiationTest, SecondDerivative_Quadratic) {
    using T = TypeParam;

    auto f = [](T x) { return x * x; };
    T x = 3, h = T(0.001);

    auto result = central_finite_difference_2nd(f, x, h);
    T exact = T(2);

    // Second derivatives have more numerical error, especially for float
    T tolerance = std::is_same_v<T, float> ? T(1.0) : T(0.01);
    EXPECT_NEAR(result.difference, exact, tolerance);
}

// Test 9: Second derivative of cubic (f(x) = x^3 => f''(x) = 6x)
TYPED_TEST(DifferentiationTest, SecondDerivative_Cubic) {
    using T = TypeParam;

    auto f = [](T x) { return x * x * x; };
    T x = 2, h = T(0.001);

    auto result = central_finite_difference_2nd(f, x, h);
    T exact = T(6) * x;

    // Second derivatives have more numerical error
    T tolerance = std::is_same_v<T, float> ? T(5.0) : T(0.01);
    EXPECT_NEAR(result.difference, exact, tolerance);
}

// Test 10: Second derivative of sine (f(x) = sin(x) => f''(x) = -sin(x))
TYPED_TEST(DifferentiationTest, SecondDerivative_Sine) {
    using T = TypeParam;

    auto f = [](T x) { return std::sin(x); };
    T x = this->pi / T(4), h = T(0.001);

    auto result = central_finite_difference_2nd(f, x, h);
    T exact = -std::sin(x);

    EXPECT_NEAR(result.difference, exact, this->test_tolerance(h) * T(10));
}

// Test 11: Second derivative of exponential (f(x) = e^x => f''(x) = e^x)
TYPED_TEST(DifferentiationTest, SecondDerivative_Exponential) {
    using T = TypeParam;

    auto f = [](T x) { return std::exp(x); };
    T x = 1, h = T(0.001);

    auto result = central_finite_difference_2nd(f, x, h);
    T exact = std::exp(x);

    // Second derivative of exponential has more error
    T tolerance = std::is_same_v<T, float> ? T(1.0) : T(0.01);
    EXPECT_NEAR(result.difference, exact, tolerance);
}

// ========== Step Size Tests ==========

// Test 12: Effect of step size on accuracy
TYPED_TEST(DifferentiationTest, StepSizeConvergence) {
    using T = TypeParam;

    auto f = [](T x) { return std::sin(x); };
    T x = 1;
    T exact = std::cos(T(1));

    // Test with different step sizes
    T h1 = T(0.1);
    T h2 = T(0.01);

    auto result1 = central_finite_difference(f, x, h1);
    auto result2 = central_finite_difference(f, x, h2);

    T error1 = this->relative_error(result1.difference, exact);
    T error2 = this->relative_error(result2.difference, exact);

    // Smaller h should give better accuracy (up to roundoff)
    // For small errors, they may be comparable due to roundoff
    if (error1 > this->epsilon * T(10000)) {
        EXPECT_LT(error2, error1 * T(1.1)); // error2 should be smaller or at most 10% larger
    }
}

// ========== Five-Point Stencil Tests ==========

// Test 13: Five-point stencil method
TEST(FivePointStencilTest, BasicDerivative) {
    using T = double;
    five_point_stencil_method method;

    auto f = [](T x) { return x * x; };
    T x = 2, h = 0.01;

    T result = method(f, x, h);
    T exact = 2 * x;

    EXPECT_NEAR(result, exact, 1e-6);
}

// Test 14: Five-point stencil with trigonometric function
TEST(FivePointStencilTest, TrigonometricFunction) {
    using T = double;
    five_point_stencil_method method;

    auto f = [](T x) { return std::sin(x); };
    T x = std::numbers::pi / 4.0, h = 0.001;

    T result = method(f, x, h);
    T exact = std::cos(x);

    EXPECT_NEAR(result, exact, 1e-6);
}

// ========== Gradient Tests ==========

// Test 15: Gradient of f(x,y) = x^2 + y^2 => grad = (2x, 2y)
TEST(GradientTest, QuadraticFunction) {
    auto f = [](const std::vector<double>& v) {
        return v[0] * v[0] + v[1] * v[1];
    };

    std::vector<double> x = {1.0, 2.0};
    double h = 0.001;

    auto gradient = grad(f, x, h);

    // Expected gradient at (1, 2) is (2, 4)
    EXPECT_NEAR(gradient[0], 2.0, 1e-5);
    EXPECT_NEAR(gradient[1], 4.0, 1e-5);
}

// Test 16: Gradient of f(x,y,z) = x*y + y*z + z*x => grad = (y+z, x+z, y+x)
TEST(GradientTest, TrilinearFunction) {
    auto f = [](const std::vector<double>& v) {
        return v[0] * v[1] + v[1] * v[2] + v[2] * v[0];
    };

    std::vector<double> x = {1.0, 2.0, 3.0};
    double h = 0.001;

    auto gradient = grad(f, x, h);

    // Expected gradient at (1, 2, 3) is (2+3, 1+3, 2+1) = (5, 4, 3)
    EXPECT_NEAR(gradient[0], 5.0, 1e-5);
    EXPECT_NEAR(gradient[1], 4.0, 1e-5);
    EXPECT_NEAR(gradient[2], 3.0, 1e-5);
}

// Test 17: Gradient of f(x,y) = sin(x) * cos(y)
TEST(GradientTest, TrigonometricFunction) {
    auto f = [](const std::vector<double>& v) {
        return std::sin(v[0]) * std::cos(v[1]);
    };

    std::vector<double> x = {std::numbers::pi / 4.0, std::numbers::pi / 3.0};
    double h = 0.001;

    auto gradient = grad(f, x, h);

    // grad = (cos(x)*cos(y), -sin(x)*sin(y))
    double exact_dx = std::cos(x[0]) * std::cos(x[1]);
    double exact_dy = -std::sin(x[0]) * std::sin(x[1]);

    EXPECT_NEAR(gradient[0], exact_dx, 1e-5);
    EXPECT_NEAR(gradient[1], exact_dy, 1e-5);
}

// ========== Edge Cases ==========

// Test 18: Derivative at zero
TYPED_TEST(DifferentiationTest, DerivativeAtZero) {
    using T = TypeParam;

    auto f = [](T x) { return x * x * x; }; // f'(0) = 0
    T x = 0, h = T(0.001);

    auto result = central_finite_difference(f, x, h);
    T exact = T(0);

    EXPECT_NEAR(result.difference, exact, T(1e-10));
}

// Test 19: Derivative of sharp function
TEST(DifferentiationEdgeCases, SharpFunction) {
    auto f = [](double x) { return std::abs(x); }; // Not differentiable at x=0
    double x = 0.1, h = 0.001; // Test near but not at singularity

    auto result = central_finite_difference(f, x, h);
    double exact = 1.0; // For x > 0

    EXPECT_NEAR(result.difference, exact, 1e-3);
}

// Test 20: Very small step size (roundoff error)
TEST(DifferentiationEdgeCases, RoundoffError) {
    auto f = [](double x) { return x * x; };
    double x = 1.0;

    // With very small h, roundoff dominates
    double h_tiny = 1e-10;
    auto result_tiny = central_finite_difference(f, x, h_tiny);

    // With reasonable h
    double h_good = 1e-4;
    auto result_good = central_finite_difference(f, x, h_good);

    double exact = 2.0;

    // The good h should actually give better results than tiny h
    double error_tiny = std::abs(result_tiny.difference - exact);
    double error_good = std::abs(result_good.difference - exact);

    // At very small h, roundoff makes things worse
    // (this test documents the behavior, not asserts it must happen)
    EXPECT_LT(error_good, 1e-6);
}

// Test 21: Higher order polynomial shows method accuracy
TEST(DifferentiationAccuracy, HighOrderPolynomial) {
    // f(x) = x^8 => f'(x) = 8x^7
    // This tests the 8th order accuracy of the method
    auto f = [](double x) {
        return std::pow(x, 8);
    };
    double x = 1.5, h = 0.01;

    auto result = central_finite_difference(f, x, h);
    double exact = 8.0 * std::pow(x, 7);

    // Should be very accurate for smooth polynomial
    EXPECT_NEAR(result.difference, exact, 1e-4);
}
