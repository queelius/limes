#include <gtest/gtest.h>
#include <cmath>
#include <numeric>
#include <functional>
#include <array>
#include <limits>
#include "../include/quadrature/quadrature_rules.hpp"

using namespace calckit::quadrature;

// Helper function for testing polynomial integration
template <typename T>
T integrate_polynomial(const auto& rule, int degree) {
    // Integrate x^n from -1 to 1
    auto f = [degree](T x) -> T {
        return std::pow(x, degree);
    };

    T sum = T(0);
    for (size_t i = 0; i < rule.size(); ++i) {
        sum += rule.weight(i) * f(rule.abscissa(i));
    }
    return sum;
}

// Calculate exact integral of x^n from -1 to 1
template <typename T>
T exact_polynomial_integral(int n) {
    if (n % 2 == 1) return T(0); // Odd powers integrate to 0
    return T(2) / T(n + 1); // Even powers
}

// Test fixture for quadrature rules
template <typename T>
class QuadratureTest : public ::testing::Test {
protected:
    static constexpr T tol = std::numeric_limits<T>::epsilon() * T(100);

    void test_polynomial_exactness(const auto& rule, int max_exact_degree) {
        for (int degree = 0; degree <= max_exact_degree; ++degree) {
            T computed = integrate_polynomial<T>(rule, degree);
            T exact = exact_polynomial_integral<T>(degree);

            EXPECT_NEAR(computed, exact, tol)
                << "Failed for polynomial degree " << degree
                << ", computed = " << computed << ", exact = " << exact;
        }

        // Should not be exact for higher degrees
        if (max_exact_degree < 20) {
            int test_degree = max_exact_degree + 1;
            T computed = integrate_polynomial<T>(rule, test_degree);
            T exact = exact_polynomial_integral<T>(test_degree);

            // Allow for some accuracy but not machine precision
            T error = std::abs(computed - exact);
            EXPECT_GT(error, tol * T(10))
                << "Unexpectedly exact for degree " << test_degree;
        }
    }

    void test_weight_sum(const auto& rule) {
        T sum = T(0);
        for (size_t i = 0; i < rule.size(); ++i) {
            sum += rule.weight(i);
        }
        // Weights should sum to 2 for interval [-1, 1]
        EXPECT_NEAR(sum, T(2), tol);
    }

    void test_node_bounds(const auto& rule) {
        for (size_t i = 0; i < rule.size(); ++i) {
            T node = rule.abscissa(i);
            EXPECT_GE(node, T(-1));
            EXPECT_LE(node, T(1));
        }
    }

    void test_symmetry(const auto& rule) {
        size_t n = rule.size();
        for (size_t i = 0; i < n / 2; ++i) {
            // Nodes should be symmetric
            EXPECT_NEAR(rule.abscissa(i), -rule.abscissa(n - 1 - i), tol);
            // Weights should be symmetric
            EXPECT_NEAR(rule.weight(i), rule.weight(n - 1 - i), tol);
        }
    }
};

using FloatTypes = ::testing::Types<float, double, long double>;
TYPED_TEST_SUITE(QuadratureTest, FloatTypes);

// Test Gauss-Legendre rules
TYPED_TEST(QuadratureTest, GaussLegendre2) {
    using T = TypeParam;
    gauss_legendre<T, 2> rule;

    EXPECT_EQ(rule.size(), 2);
    this->test_weight_sum(rule);
    this->test_node_bounds(rule);
    this->test_symmetry(rule);
    this->test_polynomial_exactness(rule, 3); // 2n-1 = 3
}

TYPED_TEST(QuadratureTest, GaussLegendre3) {
    using T = TypeParam;
    gauss_legendre<T, 3> rule;

    EXPECT_EQ(rule.size(), 3);
    this->test_weight_sum(rule);
    this->test_node_bounds(rule);
    this->test_symmetry(rule);
    this->test_polynomial_exactness(rule, 5); // 2n-1 = 5
}

// Gauss-Legendre 4 not provided in header, skip this test

TYPED_TEST(QuadratureTest, GaussLegendre5) {
    using T = TypeParam;
    gauss_legendre<T, 5> rule;

    EXPECT_EQ(rule.size(), 5);
    this->test_weight_sum(rule);
    this->test_node_bounds(rule);
    this->test_symmetry(rule);
    this->test_polynomial_exactness(rule, 9); // 2n-1 = 9
}

// Test Gauss-Kronrod rules
TYPED_TEST(QuadratureTest, GaussKronrod15) {
    using T = TypeParam;
    gauss_kronrod_15<T> rule;

    EXPECT_EQ(rule.size(), 15);
    this->test_weight_sum(rule);
    this->test_node_bounds(rule);
    this->test_symmetry(rule);

    // G-K 15 should integrate polynomials up to degree 22 exactly
    this->test_polynomial_exactness(rule, 22);
}

// Test Gauss-Kronrod embedded rules
TYPED_TEST(QuadratureTest, GaussKronrodEmbedded) {
    using T = TypeParam;
    gauss_kronrod_15<T> rule;

    // Test embedded Gauss rule
    EXPECT_EQ(rule.gauss_size, 7);

    // Test weight sums
    T gauss_weight_sum = T(0);
    for (size_t i = 0; i < rule.gauss_size; ++i) {
        gauss_weight_sum += rule.gauss_weights[i];
    }
    EXPECT_NEAR(gauss_weight_sum, T(2), this->tol);

    T kronrod_weight_sum = T(0);
    for (size_t i = 0; i < rule.size(); ++i) {
        kronrod_weight_sum += rule.weight(i);
    }
    EXPECT_NEAR(kronrod_weight_sum, T(2), this->tol);
}

// Test Clenshaw-Curtis rule
TYPED_TEST(QuadratureTest, ClenshawCurtis) {
    using T = TypeParam;

    // Test fixed size Clenshaw-Curtis
    clenshaw_curtis<T, 9> rule;

    EXPECT_EQ(rule.size(), 9);
    this->test_weight_sum(rule);
    this->test_node_bounds(rule);
    this->test_symmetry(rule);

    // Clenshaw-Curtis integrates polynomials up to degree n-1 exactly
    this->test_polynomial_exactness(rule, 8);
}

// Test Simpson's rule
TYPED_TEST(QuadratureTest, SimpsonsRule) {
    using T = TypeParam;
    simpson_rule<T> rule;

    EXPECT_EQ(rule.size(), 3);

    // Check specific nodes and weights
    EXPECT_NEAR(rule.abscissa(0), T(-1), this->tol);
    EXPECT_NEAR(rule.abscissa(1), T(0), this->tol);
    EXPECT_NEAR(rule.abscissa(2), T(1), this->tol);

    EXPECT_NEAR(rule.weight(0), T(1)/T(3), this->tol);
    EXPECT_NEAR(rule.weight(1), T(4)/T(3), this->tol);
    EXPECT_NEAR(rule.weight(2), T(1)/T(3), this->tol);

    this->test_weight_sum(rule);
    this->test_polynomial_exactness(rule, 3); // Simpson's rule is exact for cubics
}

// Test Trapezoidal rule
TYPED_TEST(QuadratureTest, TrapezoidalRule) {
    using T = TypeParam;
    trapezoidal_rule<T> rule;

    EXPECT_EQ(rule.size(), 2);

    // Check specific nodes and weights
    EXPECT_NEAR(rule.abscissa(0), T(-1), this->tol);
    EXPECT_NEAR(rule.abscissa(1), T(1), this->tol);

    EXPECT_NEAR(rule.weight(0), T(1), this->tol);
    EXPECT_NEAR(rule.weight(1), T(1), this->tol);

    this->test_weight_sum(rule);
    this->test_polynomial_exactness(rule, 1); // Trapezoidal rule is exact for linear
}

// Test Midpoint rule
TYPED_TEST(QuadratureTest, MidpointRule) {
    using T = TypeParam;
    midpoint_rule<T> rule;

    EXPECT_EQ(rule.size(), 1);

    // Check specific node and weight
    EXPECT_NEAR(rule.abscissa(0), T(0), this->tol);
    EXPECT_NEAR(rule.weight(0), T(2), this->tol);

    this->test_weight_sum(rule);
    this->test_polynomial_exactness(rule, 1); // Midpoint rule is exact for linear
}

// Newton-Cotes rules not in header, skip this test

// Chebyshev quadrature not in header, skip this test

// Test Tanh-Sinh quadrature nodes
TYPED_TEST(QuadratureTest, TanhSinhQuadrature) {
    using T = TypeParam;

    tanh_sinh_nodes<T> nodes;

    // Test a few nodes at level 3
    size_t level = 3;
    for (size_t i = 0; i < 5; ++i) {
        T x = nodes.abscissa(level, i);
        T w = nodes.weight(level, i);

        // Nodes should be in [-1, 1]
        EXPECT_GE(x, T(-1));
        EXPECT_LE(x, T(1));

        // Weights should be positive
        EXPECT_GT(w, T(0));
    }
}

// Test integration of specific functions
TEST(QuadratureIntegration, ExponentialFunction) {
    // Integrate exp(x) from -1 to 1
    auto f = [](double x) { return std::exp(x); };
    double exact = std::exp(1.0) - std::exp(-1.0);

    // Test with various rules
    {
        gauss_legendre<double, 5> rule;
        double sum = 0.0;
        for (size_t i = 0; i < rule.size(); ++i) {
            sum += rule.weight(i) * f(rule.abscissa(i));
        }
        EXPECT_NEAR(sum, exact, 1e-10);
    }

    {
        gauss_kronrod_15<double> rule;
        double sum = 0.0;
        for (size_t i = 0; i < rule.size(); ++i) {
            sum += rule.weight(i) * f(rule.abscissa(i));
        }
        EXPECT_NEAR(sum, exact, 1e-14);
    }
}

TEST(QuadratureIntegration, TrigonometricFunction) {
    // Integrate sin(pi*x) from -1 to 1 (should be 0)
    auto f = [](double x) { return std::sin(M_PI * x); };
    double exact = 0.0;

    gauss_legendre<double, 3> rule;
    double sum = 0.0;
    for (size_t i = 0; i < rule.size(); ++i) {
        sum += rule.weight(i) * f(rule.abscissa(i));
    }
    EXPECT_NEAR(sum, exact, 1e-14);
}

TEST(QuadratureIntegration, RationalFunction) {
    // Integrate 1/(1+x^2) from -1 to 1
    auto f = [](double x) { return 1.0 / (1.0 + x * x); };
    double exact = M_PI / 2.0;

    gauss_legendre<double, 5> rule;
    double sum = 0.0;
    for (size_t i = 0; i < rule.size(); ++i) {
        sum += rule.weight(i) * f(rule.abscissa(i));
    }
    EXPECT_NEAR(sum, exact, 1e-12);
}

// Test error estimation with embedded rules
TEST(QuadratureErrorEstimation, GaussKronrodError) {
    // Function with known integral
    auto f = [](double x) { return std::exp(-x * x); };

    gauss_kronrod_15<double> rule;

    // Compute with embedded Gauss rule
    double gauss_result = 0.0;
    for (size_t i = 0; i < rule.gauss_size; ++i) {
        size_t idx = rule.gauss_indices[i];
        gauss_result += rule.gauss_weights[i] * f(rule.abscissa(idx));
    }

    // Compute with Kronrod rule
    double kronrod_result = 0.0;
    for (size_t i = 0; i < rule.size(); ++i) {
        kronrod_result += rule.weight(i) * f(rule.abscissa(i));
    }

    // Error estimate
    double error_estimate = std::abs(kronrod_result - gauss_result);

    // Both should be accurate, so error should be small
    EXPECT_LT(error_estimate, 1e-10);

    // Kronrod should be more accurate
    double exact = 1.493648265624854; // Approximate value
    double gauss_error = std::abs(gauss_result - exact);
    double kronrod_error = std::abs(kronrod_result - exact);
    EXPECT_LT(kronrod_error, gauss_error);
}

// Test adaptive quadrature scenarios
TEST(QuadratureAdaptive, SingularFunction) {
    // Function with singularity at endpoint: sqrt(1-x^2)
    auto f = [](double x) -> double {
        if (std::abs(x) >= 1.0) return 0.0;
        return std::sqrt(1.0 - x * x);
    };

    // Exact integral is pi/2 (quarter circle area * 2)
    double exact = M_PI / 2.0;

    // Use tanh-sinh nodes
    tanh_sinh_nodes<double> nodes;
    double sum = 0.0;
    size_t level = 5;
    double h = 1.0 / std::pow(2.0, level);

    for (size_t i = 0; i < 50; ++i) {
        double x = nodes.abscissa(level, i);
        if (std::abs(x) < 1.0) {
            sum += nodes.weight(level, i) * f(x);
        }
        if (i > 0) {
            x = -nodes.abscissa(level, i);
            if (std::abs(x) < 1.0) {
                sum += nodes.weight(level, i) * f(x);
            }
        }
    }

    EXPECT_NEAR(sum, exact, 0.1); // Looser tolerance for this test
}

// Test quadrature rule factory/selection
TEST(QuadratureSelection, RuleSelection) {
    // Test that we can create rules dynamically
    auto create_rule = [](int order) {
        if (order <= 5) {
            return std::make_unique<gauss_legendre<double, 5>>();
        } else {
            return std::make_unique<gauss_legendre<double, 5>>();
        }
    };

    auto rule = create_rule(3);
    EXPECT_EQ(rule->size(), 5);
}

// Test custom quadrature rule
TEST(CustomQuadrature, UserDefinedRule) {
    // Create a custom 3-point rule
    struct custom_rule {
        static constexpr size_t size() { return 3; }

        double abscissa(size_t i) const {
            static const double nodes[3] = {-0.5, 0.0, 0.5};
            return nodes[i];
        }

        double weight(size_t i) const {
            static const double weights[3] = {2.0/3.0, 2.0/3.0, 2.0/3.0};
            return weights[i];
        }
    };

    custom_rule rule;

    // Test weight sum
    double sum = 0.0;
    for (size_t i = 0; i < rule.size(); ++i) {
        sum += rule.weight(i);
    }
    EXPECT_NEAR(sum, 2.0, 1e-14);
}

// Benchmark different quadrature rules (disabled by default)
TEST(QuadratureBenchmark, DISABLED_CompareRules) {
    // Complex oscillatory function
    auto f = [](double x) {
        return std::sin(10.0 * x) * std::exp(-x * x);
    };

    struct result {
        std::string name;
        double value;
        size_t points;
        double time_us;
    };

    std::vector<result> results;

    auto benchmark = [&](auto& rule, const std::string& name) {
        auto start = std::chrono::high_resolution_clock::now();

        double sum = 0.0;
        for (size_t i = 0; i < rule.size(); ++i) {
            sum += rule.weight(i) * f(rule.abscissa(i));
        }

        auto end = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);

        results.push_back({name, sum, rule.size(),
                          static_cast<double>(duration.count())});
    };

    gauss_legendre<double, 3> gl3;
    gauss_legendre<double, 5> gl5;
    gauss_kronrod_15<double> gk15;
    clenshaw_curtis<double, 17> cc17;
    simpson_rule<double> simp;
    trapezoidal_rule<double> trap;

    benchmark(gl3, "Gauss-Legendre 3");
    benchmark(gl5, "Gauss-Legendre 5");
    benchmark(gk15, "Gauss-Kronrod 15");
    benchmark(cc17, "Clenshaw-Curtis 17");
    benchmark(simp, "Simpson's Rule");
    benchmark(trap, "Trapezoidal Rule");

    std::cout << "\nQuadrature Rule Comparison:\n";
    std::cout << "----------------------------\n";
    for (const auto& r : results) {
        std::cout << std::setw(20) << r.name << ": "
                  << "value = " << std::setw(12) << r.value
                  << ", points = " << std::setw(3) << r.points
                  << ", time = " << std::setw(6) << r.time_us << " us\n";
    }
}