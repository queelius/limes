#include <gtest/gtest.h>
#include <cmath>
#include <numbers>
#include <limits>
#include "../include/ode/euler_ode1.hpp"
#include "../include/ode/rk4_ode1.hpp"
#include "../include/ode/euler_ode2.hpp"
#include "../include/ode/rk4_ode2.hpp"
#include "../include/accumulators/accumulators.hpp"

using namespace alex::math;
using namespace calckit::accumulators;

// Test fixture for ODE solvers
template <typename T>
class ODESolverTest : public ::testing::Test {
protected:
    static constexpr T epsilon = std::numeric_limits<T>::epsilon();
    static constexpr T pi = std::numbers::pi_v<T>;

    // Calculate relative error
    T relative_error(T computed, T exact) const {
        if (std::abs(exact) < epsilon) return std::abs(computed - exact);
        return std::abs(computed - exact) / std::abs(exact);
    }

    // Test-specific tolerance based on type
    T test_tolerance() const {
        if constexpr (std::is_same_v<T, float>) {
            return T(1e-3);
        } else if constexpr (std::is_same_v<T, double>) {
            return T(1e-6);
        } else {
            return T(1e-8);
        }
    }
};

using FloatTypes = ::testing::Types<float, double>;
TYPED_TEST_SUITE(ODESolverTest, FloatTypes);

// ========== First-Order ODE Tests ==========

// Test 1: y' = 0, y(0) = 1 => y(t) = 1 (constant)
TYPED_TEST(ODESolverTest, ConstantSolution_Euler) {
    using T = TypeParam;
    euler_ode1<T> solver;

    auto f = [](T t, T y) { return T(0); };
    T t0 = 0, y0 = 1, t = 1;

    auto result = solver(f, t0, y0, t, T(0.1));
    T exact = T(1);

    EXPECT_NEAR(result.value, exact, this->epsilon * T(10));
}

// Test 2: y' = y, y(0) = 1 => y(t) = e^t (exponential growth)
TYPED_TEST(ODESolverTest, ExponentialGrowth_Euler) {
    using T = TypeParam;
    euler_ode1<T> solver;

    auto f = [](T t, T y) { return y; };
    T t0 = 0, y0 = 1, t = 1;

    auto result = solver(f, t0, y0, t, T(0.01));
    T exact = std::exp(T(1));

    // Euler has first-order error, so tolerance is larger
    EXPECT_NEAR(result.value, exact, T(0.05));
}

// Test 3: y' = -y, y(0) = 1 => y(t) = e^(-t) (exponential decay)
TYPED_TEST(ODESolverTest, ExponentialDecay_RK4) {
    using T = TypeParam;
    using Acc = simple_accumulator<T>;
    rk4_ode1<Acc> solver;

    auto f = [](T t, T y) { return -y; };
    T t0 = 0, y0 = 1, t = 1;

    auto result = solver(f, t0, y0, t, T(0.1));
    T exact = std::exp(T(-1));

    // RK4 is fourth-order accurate, much better than Euler
    EXPECT_NEAR(result.value, exact, this->test_tolerance());
}

// Test 4: y' = t, y(0) = 0 => y(t) = t^2/2 (quadratic)
TYPED_TEST(ODESolverTest, QuadraticSolution_RK4) {
    using T = TypeParam;
    using Acc = simple_accumulator<T>;
    rk4_ode1<Acc> solver;

    auto f = [](T t, T y) { return t; };
    T t0 = 0, y0 = 0, t = 2;

    auto result = solver(f, t0, y0, t, T(0.1));
    T exact = T(2) * T(2) / T(2); // t^2/2 at t=2

    EXPECT_NEAR(result.value, exact, this->test_tolerance());
}

// Test 5: y' = cos(t), y(0) = 0 => y(t) = sin(t) (trigonometric)
TYPED_TEST(ODESolverTest, TrigonometricSolution_RK4) {
    using T = TypeParam;
    using Acc = neumaier_accumulator<T>;
    rk4_ode1<Acc> solver;

    auto f = [](T t, T y) { return std::cos(t); };
    T t0 = 0, y0 = 0, t = this->pi / T(2);

    auto result = solver(f, t0, y0, t, T(0.01));
    T exact = std::sin(this->pi / T(2)); // sin(π/2) = 1

    EXPECT_NEAR(result.value, exact, this->test_tolerance());
}

// Test 6: y' = t * y, y(0) = 1 => y(t) = e^(t^2/2)
TYPED_TEST(ODESolverTest, NonlinearODE_RK4) {
    using T = TypeParam;
    using Acc = kahan_accumulator<T>;
    rk4_ode1<Acc> solver;

    auto f = [](T t, T y) { return t * y; };
    T t0 = 0, y0 = 1, t = 1;

    auto result = solver(f, t0, y0, t, T(0.01));
    T exact = std::exp(T(0.5)); // e^(1^2/2) = e^0.5

    EXPECT_NEAR(result.value, exact, this->test_tolerance() * T(10));
}

// Test 7: Stiff equation y' = -10y, y(0) = 1 => y(t) = e^(-10t)
TYPED_TEST(ODESolverTest, StiffEquation_RK4) {
    using T = TypeParam;
    using Acc = klein_accumulator<T>;
    rk4_ode1<Acc> solver;

    auto f = [](T t, T y) { return T(-10) * y; };
    T t0 = 0, y0 = 1, t = T(0.5);

    // Need smaller step size for stiff equations
    auto result = solver(f, t0, y0, t, T(0.01));
    T exact = std::exp(T(-5)); // e^(-10*0.5)

    EXPECT_NEAR(result.value, exact, T(0.01));
}

// Test 8: Verify iterations count
TYPED_TEST(ODESolverTest, IterationCount_Euler) {
    using T = TypeParam;
    euler_ode1<T> solver;

    auto f = [](T t, T y) { return T(1); };
    T t0 = 0, y0 = 0, t = 1, h = T(0.1);

    auto result = solver(f, t0, y0, t, h);

    // Should take 10 steps to go from 0 to 1 with h=0.1
    EXPECT_EQ(result.iterations, 10);
}

// Test 9: Different step sizes convergence for RK4
TYPED_TEST(ODESolverTest, ConvergenceOrder_RK4) {
    using T = TypeParam;
    using Acc = simple_accumulator<T>;
    rk4_ode1<Acc> solver;

    // Test y' = y, y(0) = 1 => y(1) = e
    auto f = [](T t, T y) { return y; };
    T t0 = 0, y0 = 1, t = 1;
    T exact = std::exp(T(1));

    // Compute with two different step sizes
    T h1 = T(0.1);
    T h2 = T(0.05);

    auto result1 = solver(f, t0, y0, t, h1);
    auto result2 = solver(f, t0, y0, t, h2);

    T error1 = std::abs(result1.value - exact);
    T error2 = std::abs(result2.value - exact);

    // For RK4 (4th order), halving h should reduce error by ~2^4 = 16
    // We check for at least 2^3 = 8 to account for rounding
    if (error1 > this->epsilon * T(100)) {
        T ratio = error1 / error2;
        EXPECT_GT(ratio, T(8));
    }
}

// Test 10: Initial value preservation
TYPED_TEST(ODESolverTest, InitialValuePreservation) {
    using T = TypeParam;
    using Acc = simple_accumulator<T>;
    rk4_ode1<Acc> solver;

    auto f = [](T t, T y) { return T(0); };
    T t0 = 5, y0 = 42, t = 5; // No time evolution

    auto result = solver(f, t0, y0, t, T(0.1));

    EXPECT_EQ(result.value, y0);
    EXPECT_EQ(result.iterations, 0);
}

// ========== Edge Cases ==========

TYPED_TEST(ODESolverTest, BackwardIntegration_Euler) {
    using T = TypeParam;
    euler_ode1<T> solver;

    auto f = [](T t, T y) { return T(1); };
    T t0 = 1, y0 = 0, t = 0;

    // Backward integration (t < t0) should handle gracefully
    auto result = solver(f, t0, y0, t, T(0.1));

    // Since t < t0, no iterations should occur
    EXPECT_EQ(result.iterations, 0);
    EXPECT_EQ(result.value, y0);
}

// Test with very small time interval
TYPED_TEST(ODESolverTest, SmallTimeInterval_RK4) {
    using T = TypeParam;
    using Acc = simple_accumulator<T>;
    rk4_ode1<Acc> solver;

    auto f = [](T t, T y) { return y; };
    T t0 = 0, y0 = 1, t = T(1e-6);

    auto result = solver(f, t0, y0, t, T(1e-7));
    T exact = std::exp(T(1e-6));

    EXPECT_NEAR(result.value, exact, this->epsilon * T(1000));
}

// Test with different accumulator types
TEST(ODESolverAccumulatorTest, DifferentAccumulators) {
    auto f = [](double t, double y) { return y; };
    double t0 = 0, y0 = 1, t = 1, h = 0.01;
    double exact = std::exp(1.0);

    // Test with different accumulators
    {
        rk4_ode1<simple_accumulator<double>> solver;
        auto result = solver(f, t0, y0, t, h);
        EXPECT_NEAR(result.value, exact, 1e-6);
    }

    {
        rk4_ode1<kahan_accumulator<double>> solver;
        auto result = solver(f, t0, y0, t, h);
        EXPECT_NEAR(result.value, exact, 1e-6);
    }

    {
        rk4_ode1<neumaier_accumulator<double>> solver;
        auto result = solver(f, t0, y0, t, h);
        EXPECT_NEAR(result.value, exact, 1e-6);
    }

    {
        rk4_ode1<klein_accumulator<double>> solver;
        auto result = solver(f, t0, y0, t, h);
        EXPECT_NEAR(result.value, exact, 1e-6);
    }
}

// Test Euler vs RK4 comparison
TEST(ODESolverComparison, EulerVsRK4) {
    // Test y' = y, y(0) = 1 => y(1) = e
    auto f = [](double t, double y) { return y; };
    double t0 = 0, y0 = 1, t = 1, h = 0.1;
    double exact = std::exp(1.0);

    euler_ode1<double> euler;
    rk4_ode1<simple_accumulator<double>> rk4;

    auto euler_result = euler(f, t0, y0, t, h);
    auto rk4_result = rk4(f, t0, y0, t, h);

    double euler_error = std::abs(euler_result.value - exact);
    double rk4_error = std::abs(rk4_result.value - exact);

    // RK4 should be significantly more accurate than Euler
    EXPECT_LT(rk4_error, euler_error * 0.01);
}
