# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

Header-only C++ library for numerical integration, differentiation, and ODE solving using modern C++20 with concept-based generic programming. The library emphasizes composability, allowing mix-and-match of quadrature rules, accumulators, and transforms.

## Build System

### CMake (Recommended)
```bash
# Configure
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release

# Build
cmake --build build

# Build and run all tests
cmake --build build --target test

# Run specific test executable
./build/tests/test_integrators

# Build examples
cmake --build build --target demo

# Coverage report (Debug mode only)
cmake -S . -B build -DCMAKE_BUILD_TYPE=Debug
cmake --build build --target coverage
# Report: build/coverage/html/index.html
```

### Direct Compilation
```bash
# Build ODE test
g++ -O3 -std=c++2a -Wall -o ode_test include/ode/ode_test.cpp -I.

# Build examples
g++ -O3 -std=c++2a -Wall -o demo examples/demo.cpp -I.
g++ -O3 -std=c++2a -Wall -o simple_demo examples/simple_demo.cpp -I.

# Legacy makefile support
make integration_test
```

### Compilation Requirements
- C++20 standard (`-std=c++2a`)
- Include path: `-I.` (root directory)
- Recommended flags: `-O3 -march=native -Wall`
- Optional: `-mavx2 -mfma` for SIMD, `-fopenmp` for parallelization

## Architecture

### Composable Design Philosophy

The library uses an **orthogonal component architecture** where integration algorithms are composed from independent, pluggable parts. See ARCHITECTURE.md for detailed design rationale.

Key composition layers:
1. **Accumulators**: Control numerical precision (Simple, Kahan, Neumaier, Klein, Pairwise)
2. **Quadrature Rules**: Define integration nodes/weights (Gauss-Legendre, Gauss-Kronrod, Tanh-Sinh, etc.)
3. **Integrators**: Combine rules + accumulators into algorithms
4. **Transforms**: Handle special domains (infinite intervals, singularities)
5. **Parallel Execution**: Work-stealing parallelization layer

### Namespace Organization
- `calckit` (alias: `ai`): Main namespace for integration
- `alex::math`: ODE solvers (legacy namespace from earlier design)

### Key Header Files

**Main Entry Point:**
- `include/calckit.hpp`: High-level interface with `integrate<T>::adaptive()`, builder pattern, convenience functions

**Core Infrastructure:**
- `include/concepts/integrator_concepts.hpp`: C++20 concepts for all template parameters
- `include/core/integration_result.hpp`: Result type with value, error, convergence info

**Component Libraries:**
- `include/accumulators/accumulators.hpp`: Precision control strategies
- `include/quadrature/quadrature_rules.hpp`: Node/weight generation
- `include/integrators/univariate_integrator.hpp`: Core integration algorithms
- `include/transforms/coordinate_transforms.hpp`: Change of variables
- `include/parallel/parallel_integration.hpp`: Parallel execution policies

**Legacy Integrators (older design):**
- `include/numerical_integrators/*.hpp`: Original single-file integrators

**Differentiation:**
- `include/numerical_differentiation.hpp`: Central finite difference, gradients
- `include/central_finite_difference.hpp`: Core differentiation algorithms

**Antiderivatives:**
- `include/antiderivative.hpp`: Symbolic antiderivative computation

**ODE Solvers:**
- `include/ode/euler_ode1.hpp`, `include/ode/rk4_ode1.hpp`: First-order ODEs
- `include/ode/euler_ode2.hpp`, `include/ode/rk4_ode2.hpp`: Second-order ODEs

### Usage Patterns

**Simple Integration:**
```cpp
auto result = integrate_adaptive(f, a, b, 1e-8);
```

**Custom Composition:**
```cpp
using acc = accumulators::klein_accumulator<double>;
using rule = quadrature::gauss_kronrod_15<double>;
quadrature_integrator<double, rule, acc> integrator{rule{}, acc{}};
auto result = integrator(f, a, b, tol);
```

**Builder Pattern:**
```cpp
auto integrator = make_integrator<double>()
    .with_quadrature("gauss-legendre")
    .with_accumulator("neumaier")
    .with_parallel(true)
    .with_tolerance(1e-12);
auto result = integrator.integrate(f, a, b);
```

## Testing

### Test Structure
- Google Test framework for unit tests (118+ tests total)
- Each component has dedicated test file:
  - `test_accumulators.cpp` (42 tests): Type-parameterized float/double tests
  - `test_quadrature.cpp`: Quadrature rule validation
  - `test_integrators.cpp`: Integration algorithm tests
  - `test_ode_solvers.cpp` (26 tests): ODE solver accuracy and convergence
  - `test_differentiation.cpp` (34 tests): Derivative and gradient tests
  - `test_antiderivative.cpp` (16 tests): Antiderivative computation tests
  - `test_transforms.cpp`: Coordinate transformation tests
  - `test_parallel.cpp`: Parallel execution tests
  - `test_concepts.cpp`: C++20 concept validation
  - `test_integration_result.cpp`: Result type tests
- `basic_tests.cpp`: Simple validation without GTest dependency
- Type-parameterized tests use type-aware tolerances (1e-5 for float, 1e-10 for double)

### Running Tests
```bash
# All tests
cd build && ctest --output-on-failure

# Specific test suite
./build/tests/test_integrators
./build/tests/test_accumulators
./build/tests/test_ode_solvers
./build/tests/test_differentiation
./build/tests/test_antiderivative
./build/tests/test_quadrature
./build/tests/test_transforms
./build/tests/test_parallel
./build/tests/test_concepts
./build/tests/test_integration_result

# Run with verbose output
ctest -V

# Run tests matching pattern
ctest -R "accumulator|integrator"

# With coverage (Debug build)
cmake --build build --target coverage
```

### Coverage Requirements
When adding features, ensure test coverage using the coverage target. Coverage reports exclude test files and external dependencies.

## Design Constraints

### Type Safety
All integrators define `value_type` and `result_type`. Use concepts to validate template parameters at compile time.

### Header-Only Implementation
All functionality in headers for zero-cost abstraction through inlining. No .cpp files in library core (only tests/examples).

### Performance Features
- SIMD: AVX2 support via `-mavx2` (controlled by `ENABLE_AVX2` CMake option)
- Parallel: OpenMP support via `ENABLE_OPENMP` option
- Sanitizers: Address/undefined behavior sanitizers enabled in Debug builds

## Common Development Patterns

### Adding New Quadrature Rules
Implement in `include/quadrature/quadrature_rules.hpp`, must satisfy `QuadratureRule` concept with nodes, weights, and order.

### Adding New Accumulators
Implement in `include/accumulators/accumulators.hpp`, must provide `operator+=()`, `operator()()`, and `value_type`.

### Adding New Integrators
Compose existing rules + accumulators in `include/integrators/`, or create new algorithm following concept requirements.

### ODE Solver Pattern
ODE solvers use interval-based stepping with `alex::math` namespace. See existing Euler/RK4 implementations for patterns.

### Type-Aware Testing
When writing tests for floating-point code, use type-aware tolerances:
```cpp
TYPED_TEST(MyTest, SomeTest) {
    using T = TypeParam;
    T tolerance = std::is_same_v<T, float> ? T(1e-5) : T(1e-10);
    EXPECT_NEAR(result, expected, tolerance);
}
```

### Differentiation and Antiderivatives
- Use `central_finite_difference()` for derivatives (8th-order accurate)
- Use `grad()` for gradients of multivariable functions
- Use `antideriv()` to create antiderivative functors
- Step size `h` typically 1e-3 to 1e-5 for double precision
