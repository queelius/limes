# Testing

## Running Tests

### All Tests

```bash
cd build
ctest --output-on-failure
```

### Individual Test Suites

```bash
./tests/test_accumulators
./tests/test_integrators
./tests/test_quadrature
./tests/test_ode_solvers
./tests/test_differentiation
./tests/test_antiderivative
./tests/test_transforms
./tests/test_parallel
```

### With Verbose Output

```bash
ctest -V
```

## Test Coverage

### Generate Coverage Report

```bash
mkdir build && cd build
cmake -DCMAKE_BUILD_TYPE=Debug ..
make
make coverage
```

View the report: `build/coverage/html/index.html`

### Current Coverage

The library has 118+ unit tests covering:

- **Accumulators** (42 tests): All accumulator types with float/double precision
- **Quadrature Rules** (tests): Gauss-Legendre, Gauss-Kronrod, Simpson
- **Integrators** (tests): Adaptive integration, transforms, edge cases
- **ODE Solvers** (26 tests): Euler and RK4 for 1st and 2nd order ODEs
- **Differentiation** (34 tests): First/second derivatives, gradients
- **Antiderivatives** (16 tests): Fundamental theorem, initial value problems
- **Transforms** (tests): Infinite intervals, singularities
- **Parallel** (tests): OpenMP parallelization

## Test Framework

Tests use **Google Test** framework with parameterized tests for type testing:

```cpp
template <typename T>
class IntegratorTest : public ::testing::Test {};

using FloatTypes = ::testing::Types<double, float>;
TYPED_TEST_SUITE(IntegratorTest, FloatTypes);

TYPED_TEST(IntegratorTest, BasicIntegration) {
    using T = TypeParam;
    // Test works for both float and double
}
```

## Writing Tests

### Basic Test

```cpp
TEST(IntegrationTest, SimplePolynomial) {
    auto f = [](double x) { return x * x; };
    auto result = integrate_adaptive(f, 0.0, 1.0, 1e-10);

    EXPECT_NEAR(result.value, 1.0/3.0, 1e-8);
}
```

### Parameterized Test

```cpp
TYPED_TEST(AccumulatorTest, LinearSum) {
    using T = TypeParam;
    using Acc = kahan_accumulator<T>;

    Acc acc;
    for (int i = 0; i < 100; i++) {
        acc.add(T(0.1));
    }

    T tolerance = std::is_same_v<T, float> ? T(1e-4) : T(1e-10);
    EXPECT_NEAR(acc.result(), T(10), tolerance);
}
```

## Continuous Integration

Tests run automatically on:
- Push to main branch
- Pull requests
- Tagged releases

## Benchmarks

Performance benchmarks are in `benchmarks/`:

```bash
./benchmarks/bench_integration
./benchmarks/bench_accumulators
./benchmarks/bench_ode
```

Results show:
- Function evaluation counts
- Wall-clock time
- Accuracy vs. performance tradeoffs
