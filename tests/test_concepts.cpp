#include <gtest/gtest.h>
#include <type_traits>
#include <complex>
#include <functional>
#include <limes/algorithms/concepts/concepts.hpp>
#include <limes/algorithms/accumulators/accumulators.hpp>
#include <limes/algorithms/core/result.hpp>

using namespace limes::algorithms::concepts;

// Basic concept tests
TEST(ConceptsTest, BasicConcepts) {
    // Test that basic types satisfy Field concept
    EXPECT_TRUE(Field<double>);
    EXPECT_TRUE(Field<float>);
    EXPECT_TRUE(Field<long double>);

    // Complex numbers should also work
    EXPECT_TRUE(Field<std::complex<double>>);
}

// Test function concepts
TEST(ConceptsTest, FunctionConcepts) {
    // Lambda function
    auto lambda = [](double x) { return x * x; };
    EXPECT_TRUE((UnivariateFunction<decltype(lambda), double>));

    // Function pointer
    double (*func_ptr)(double) = [](double x) { return std::sin(x); };
    EXPECT_TRUE((UnivariateFunction<decltype(func_ptr), double>));

    // std::function
    std::function<double(double)> std_func = [](double x) { return x; };
    EXPECT_TRUE((UnivariateFunction<decltype(std_func), double>));
}

// Test that concepts work in practice
TEST(ConceptsRuntime, BasicUsage) {
    // Test Field operations
    if constexpr (Field<double>) {
        double a = 3.0;
        double b = 2.0;
        double sum = a + b;
        double diff = a - b;
        double prod = a * b;
        double quot = a / b;

        EXPECT_EQ(sum, 5.0);
        EXPECT_EQ(diff, 1.0);
        EXPECT_EQ(prod, 6.0);
        EXPECT_EQ(quot, 1.5);
    }

    // Test UnivariateFunction
    if constexpr (UnivariateFunction<double(*)(double), double>) {
        double (*f)(double) = [](double x) { return 2.0 * x; };
        double result = f(3.0);
        EXPECT_EQ(result, 6.0);
    }
}

// Test concept-constrained templates
TEST(ConceptsTemplates, ConstrainedTemplates) {
    // Simple field-constrained function
    auto add_fields = []<typename T>(T a, T b) -> T
        requires Field<T> {
        return a + b;
    };

    double result = add_fields(3.0, 4.0);
    EXPECT_EQ(result, 7.0);

    // Function-constrained template
    auto evaluate_at_zero = []<typename F>(F f) -> double
        requires UnivariateFunction<F, double> {
        return f(0.0);
    };

    auto constant_func = [](double) { return 5.0; };
    double value = evaluate_at_zero(constant_func);
    EXPECT_EQ(value, 5.0);
}

// Test type traits with concepts
TEST(ConceptsTraits, TypeChecking) {
    // Check Field types
    EXPECT_TRUE(Field<double>);
    EXPECT_TRUE(Field<float>);

    // Check function types
    struct Functor {
        double operator()(double x) const { return x; }
    };

    EXPECT_TRUE((UnivariateFunction<Functor, double>));

    // Lambda
    auto lambda = [](double x) { return x * x; };
    EXPECT_TRUE((UnivariateFunction<decltype(lambda), double>));
}

// Test practical integration usage
TEST(ConceptsIntegration, PracticalUsage) {
    // Use concepts to ensure proper types
    if constexpr (Field<double>) {
        // This code only compiles if double is a Field
        double x = 1.0;
        double y = 2.0;
        double z = x + y;
        EXPECT_EQ(z, 3.0);
    }

    // Test with function
    auto square = [](double x) { return x * x; };
    if constexpr (UnivariateFunction<decltype(square), double>) {
        double result = square(4.0);
        EXPECT_EQ(result, 16.0);
    }
}

// Verify concepts work with the library types
TEST(ConceptsLibrary, LibraryTypes) {
    // Test that library accumulators satisfy the concept
    using SimpleAcc = limes::algorithms::accumulators::simple_accumulator<double>;
    EXPECT_TRUE((Accumulator<SimpleAcc, double>));

    // Test that integration result works with Field concept
    EXPECT_TRUE(Field<double>);
    limes::algorithms::integration_result<double> result(1.0, 1e-6, 100);
    EXPECT_EQ(result.value(), 1.0);
}
