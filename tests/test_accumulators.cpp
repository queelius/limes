#include <gtest/gtest.h>
#include <numeric>
#include <random>
#include <vector>
#include <cmath>
#include <limits>
#include <array>
#include "../include/accumulators/accumulators.hpp"

using namespace algebraic_integrators::accumulators;

// Test fixture for accumulator tests
template <typename T>
class AccumulatorTest : public ::testing::Test {
protected:
    static constexpr T epsilon = std::numeric_limits<T>::epsilon();

    // Generate test data with known sum
    std::vector<T> generate_uniform_data(size_t n, T value) {
        return std::vector<T>(n, value);
    }

    // Generate problematic data for testing compensation
    std::vector<T> generate_ill_conditioned_data(size_t n) {
        std::vector<T> data;
        data.push_back(T(1.0));
        for (size_t i = 1; i < n; ++i) {
            data.push_back(epsilon / T(2));
        }
        return data;
    }

    // Generate random data
    std::vector<T> generate_random_data(size_t n, T min_val, T max_val) {
        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_real_distribution<T> dist(min_val, max_val);

        std::vector<T> data(n);
        for (auto& val : data) {
            val = dist(gen);
        }
        return data;
    }

    // Calculate relative error
    T relative_error(T computed, T exact) {
        if (exact == T(0)) return std::abs(computed);
        return std::abs(computed - exact) / std::abs(exact);
    }
};

// Use multiple floating-point types
using FloatTypes = ::testing::Types<float, double, long double>;
TYPED_TEST_SUITE(AccumulatorTest, FloatTypes);

// Test simple accumulator
TYPED_TEST(AccumulatorTest, SimpleAccumulatorBasic) {
    using T = TypeParam;
    simple_accumulator<T> acc;

    EXPECT_EQ(acc(), T(0));

    acc += T(1.0);
    EXPECT_EQ(acc(), T(1.0));

    acc += T(2.0);
    EXPECT_EQ(acc(), T(3.0));

    acc -= T(1.0);
    EXPECT_EQ(acc(), T(2.0));
}

TYPED_TEST(AccumulatorTest, SimpleAccumulatorRoundingError) {
    using T = TypeParam;
    simple_accumulator<T> acc;

    // Add many small numbers
    const size_t n = 1000000;
    const T tiny = std::numeric_limits<T>::epsilon();

    for (size_t i = 0; i < n; ++i) {
        acc += tiny;
    }

    T result = acc();
    T exact = tiny * n;

    // Simple accumulator will have significant rounding error
    T rel_error = this->relative_error(result, exact);

    // We expect error, just verify it's reasonable
    EXPECT_LT(rel_error, T(1.0)); // Error should be less than 100%
}

// Test Kahan accumulator
TYPED_TEST(AccumulatorTest, KahanAccumulatorBasic) {
    using T = TypeParam;
    kahan_accumulator<T> acc;

    EXPECT_EQ(acc(), T(0));

    acc += T(1.0);
    EXPECT_NEAR(acc(), T(1.0), this->epsilon);

    acc += T(2.0);
    EXPECT_NEAR(acc(), T(3.0), this->epsilon * T(3));
}

TYPED_TEST(AccumulatorTest, KahanAccumulatorCompensation) {
    using T = TypeParam;
    kahan_accumulator<T> acc;

    // Test case where compensation is crucial
    auto data = this->generate_ill_conditioned_data(100000);
    T exact = T(1.0) + (data.size() - 1) * (this->epsilon / T(2));

    for (const auto& val : data) {
        acc += val;
    }

    T result = acc();
    T rel_error = this->relative_error(result, exact);

    // Kahan should have much better accuracy than simple
    EXPECT_LT(rel_error, T(1e-10));
}

// Test Neumaier accumulator
TYPED_TEST(AccumulatorTest, NeumaierAccumulatorBasic) {
    using T = TypeParam;
    neumaier_accumulator<T> acc;

    EXPECT_EQ(acc(), T(0));

    // Test with sorted values (where Neumaier excels)
    std::vector<T> values = {T(1e10), T(1.0), T(1e-10)};

    for (const auto& val : values) {
        acc += val;
    }

    T exact = T(1e10) + T(1.0) + T(1e-10);
    EXPECT_NEAR(acc(), exact, exact * this->epsilon * T(10));
}

TYPED_TEST(AccumulatorTest, NeumaierAccumulatorLargeSmallMix) {
    using T = TypeParam;
    neumaier_accumulator<T> acc;

    // Mix of large and small values
    acc += T(1e15);
    for (int i = 0; i < 1000; ++i) {
        acc += T(1.0);
    }
    acc += T(1e-15);

    T exact = T(1e15) + T(1000.0) + T(1e-15);
    T result = acc();

    // Neumaier should handle this well
    EXPECT_NEAR(result, exact, exact * this->epsilon * T(100));
}

// Test Klein accumulator
TYPED_TEST(AccumulatorTest, KleinAccumulatorBasic) {
    using T = TypeParam;
    klein_accumulator<T> acc;

    EXPECT_EQ(acc(), T(0));

    acc += T(1.0);
    EXPECT_EQ(acc(), T(1.0));

    acc += T(2.0);
    EXPECT_EQ(acc(), T(3.0));
}

TYPED_TEST(AccumulatorTest, KleinAccumulatorHighPrecision) {
    using T = TypeParam;
    klein_accumulator<T> acc;

    // Klein uses second-order error compensation
    const size_t n = 10000;
    std::vector<T> data = this->generate_random_data(n, T(-100), T(100));

    // Calculate exact sum using higher precision if available
    long double exact_sum = 0.0L;
    for (const auto& val : data) {
        exact_sum += static_cast<long double>(val);
    }

    for (const auto& val : data) {
        acc += val;
    }

    T result = acc();
    T rel_error = this->relative_error(result, static_cast<T>(exact_sum));

    // Klein should have excellent accuracy
    EXPECT_LT(rel_error, T(1e-12));
}

// Test Pairwise accumulator
TYPED_TEST(AccumulatorTest, PairwiseAccumulatorBasic) {
    using T = TypeParam;
    pairwise_accumulator<T> acc;

    EXPECT_EQ(acc(), T(0));

    // Add powers of 2 (exact in floating point)
    for (int i = 0; i < 10; ++i) {
        acc += std::pow(T(2), i);
    }

    T exact = std::pow(T(2), 10) - T(1); // Sum of geometric series
    EXPECT_NEAR(acc(), exact, exact * this->epsilon);
}

TYPED_TEST(AccumulatorTest, PairwiseAccumulatorLargeDataset) {
    using T = TypeParam;
    pairwise_accumulator<T> acc;

    // Large dataset test
    const size_t n = 100000;
    std::vector<T> data = this->generate_uniform_data(n, T(0.1));

    for (const auto& val : data) {
        acc += val;
    }

    T exact = T(n) * T(0.1);
    T result = acc();

    EXPECT_NEAR(result, exact, exact * this->epsilon * std::log2(n));
}

// Test pairwise accumulator finalize functionality
TYPED_TEST(AccumulatorTest, PairwiseAccumulatorSum) {
    using T = TypeParam;
    pairwise_accumulator<T> acc;

    // Add some values
    for (int i = 1; i <= 100; ++i) {
        acc += T(i);
    }

    T result = acc();

    T exact = T(100 * 101 / 2); // Sum of 1 to 100
    EXPECT_NEAR(result, exact, exact * this->epsilon * T(10));
}

// Comparative accuracy test
TYPED_TEST(AccumulatorTest, ComparativeAccuracy) {
    using T = TypeParam;

    // Generate challenging data
    std::vector<T> data;
    data.push_back(T(1e10));
    for (int i = 0; i < 10000; ++i) {
        data.push_back(T(1.0));
    }
    data.push_back(T(-1e10));

    T exact = T(10000.0);

    // Test all accumulators
    simple_accumulator<T> simple;
    kahan_accumulator<T> kahan;
    neumaier_accumulator<T> neumaier;
    klein_accumulator<T> klein;
    pairwise_accumulator<T> pairwise;

    for (const auto& val : data) {
        simple += val;
        kahan += val;
        neumaier += val;
        klein += val;
        pairwise += val;
    }


    // Calculate errors
    T simple_error = this->relative_error(simple(), exact);
    T kahan_error = this->relative_error(kahan(), exact);
    T neumaier_error = this->relative_error(neumaier(), exact);
    T klein_error = this->relative_error(klein(), exact);
    T pairwise_error = this->relative_error(pairwise(), exact);

    // Verify hierarchy of accuracy (generally expected)
    // Klein and Neumaier should be best, followed by Kahan and Pairwise, then Simple
    EXPECT_GE(simple_error, kahan_error * T(0.1));
    EXPECT_LE(klein_error, kahan_error);
    EXPECT_LE(neumaier_error, kahan_error);

    // All compensated methods should be significantly better than simple
    EXPECT_LT(kahan_error, simple_error * T(0.5));
    EXPECT_LT(neumaier_error, simple_error * T(0.5));
    EXPECT_LT(klein_error, simple_error * T(0.5));
}

// Test edge cases
TYPED_TEST(AccumulatorTest, EdgeCases) {
    using T = TypeParam;

    // Test with NaN
    {
        simple_accumulator<T> acc;
        acc += std::numeric_limits<T>::quiet_NaN();
        EXPECT_TRUE(std::isnan(acc()));
    }

    // Test with infinity
    {
        kahan_accumulator<T> acc;
        acc += std::numeric_limits<T>::infinity();
        EXPECT_TRUE(std::isinf(acc()));
        EXPECT_GT(acc(), T(0));
    }

    // Test overflow behavior
    {
        neumaier_accumulator<T> acc;
        T max_val = std::numeric_limits<T>::max();
        acc += max_val / T(2);
        acc += max_val / T(2);
        acc += T(1.0);
        // Result should be close to max_val + 1
        EXPECT_GT(acc(), max_val * T(0.99));
    }

    // Test underflow behavior
    {
        klein_accumulator<T> acc;
        T min_normal = std::numeric_limits<T>::min();
        for (int i = 0; i < 1000; ++i) {
            acc += min_normal;
        }
        EXPECT_GT(acc(), T(0));
        EXPECT_NEAR(acc(), min_normal * T(1000), min_normal * T(1000) * this->epsilon * T(100));
    }
}

// Test accumulator reset functionality
TEST(AccumulatorResetTest, ResetFunctionality) {
    // Test that accumulators can be reused after getting result
    kahan_accumulator<double> acc;

    acc += 1.0;
    acc += 2.0;
    EXPECT_EQ(acc(), 3.0);

    // Reset by creating new instance
    acc = kahan_accumulator<double>();
    EXPECT_EQ(acc(), 0.0);

    acc += 5.0;
    EXPECT_EQ(acc(), 5.0);
}

// Performance comparison test (not a unit test, but useful for benchmarking)
TEST(AccumulatorPerformance, DISABLED_PerformanceBenchmark) {
    const size_t n = 10000000;
    std::vector<double> data(n);

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<double> dist(-1.0, 1.0);

    for (auto& val : data) {
        val = dist(gen);
    }

    // Benchmark different accumulators
    auto benchmark = [&data](auto& acc, const std::string& name) {
        auto start = std::chrono::high_resolution_clock::now();

        for (const auto& val : data) {
            acc += val;
        }

        // No finalize needed

        auto end = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);

        std::cout << name << ": " << duration.count() << " us, result = "
                  << acc() << std::endl;
    };

    simple_accumulator<double> simple;
    kahan_accumulator<double> kahan;
    neumaier_accumulator<double> neumaier;
    klein_accumulator<double> klein;
    pairwise_accumulator<double> pairwise;

    benchmark(simple, "Simple");
    benchmark(kahan, "Kahan");
    benchmark(neumaier, "Neumaier");
    benchmark(klein, "Klein");
    benchmark(pairwise, "Pairwise");
}

// Test accumulator with different numeric types
TEST(AccumulatorTypeTest, IntegerTypes) {
    // Accumulators should work with integer types too
    simple_accumulator<int> int_acc;
    int_acc += 5;
    int_acc += 10;
    EXPECT_EQ(int_acc(), 15);

    simple_accumulator<long long> ll_acc;
    ll_acc += 1000000000LL;
    ll_acc += 2000000000LL;
    EXPECT_EQ(ll_acc(), 3000000000LL);
}

// Test custom operations
TEST(AccumulatorCustomOps, MultiplyAccumulate) {
    // Test multiply-accumulate pattern
    kahan_accumulator<double> acc;

    // Simulate dot product
    std::vector<double> a = {1.0, 2.0, 3.0};
    std::vector<double> b = {4.0, 5.0, 6.0};

    for (size_t i = 0; i < a.size(); ++i) {
        acc += a[i] * b[i];
    }

    double expected = 1.0*4.0 + 2.0*5.0 + 3.0*6.0;
    EXPECT_NEAR(acc(), expected, 1e-14);
}