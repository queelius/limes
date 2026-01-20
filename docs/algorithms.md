# limes::algorithms - Numerical Integration Algorithms

The `limes::algorithms` namespace provides composable building blocks for numerical integration: **accumulators** for precision control, **quadrature rules** for node/weight generation, and **integrators** that combine them.

## Overview

```
┌─────────────────────────────────────────────────────┐
│  Integrators: quadrature_integrator, romberg, etc.  │
├──────────────────┬──────────────────────────────────┤
│  Quadrature Rules│  Accumulators                    │
│  (where to sample)│  (how to sum)                   │
├──────────────────┴──────────────────────────────────┤
│  Concepts: Field, UnivariateFunction, etc.          │
└─────────────────────────────────────────────────────┘
```

Each layer is independent and composable via C++20 concepts.

---

## Accumulators

**Namespace:** `limes::algorithms::accumulators`

Accumulators control how floating-point values are summed. Different strategies trade off speed vs. precision.

### Interface

All accumulators satisfy the `Accumulator<T>` concept:

```cpp
acc += value;      // Add a value
T sum = acc();     // Get current sum
acc.reset();       // Reset to zero
```

### Available Accumulators

| Accumulator | Description | When to Use |
|-------------|-------------|-------------|
| `simple_accumulator<T>` | Direct summation | Fast, low precision needs |
| `kahan_accumulator<T>` | Kahan-Babuska compensation | General purpose |
| `neumaier_accumulator<T>` | Improved Kahan | Mixed magnitude values |
| `klein_accumulator<T>` | Second-order compensation | Extreme precision |
| `pairwise_accumulator<T, ChunkSize>` | Divide-and-conquer | Many terms, cache-friendly |
| `any_accumulator<T>` | Type-erased wrapper | Runtime polymorphism |

### Details

**`simple_accumulator<T>`**

Standard floating-point addition. Fast but accumulates rounding errors.

```cpp
simple_accumulator<double> acc;
for (double x : values) acc += x;
double sum = acc();
```

**`kahan_accumulator<T>`**

Compensated summation that tracks and corrects rounding errors. Provides ~2x the precision of simple summation.

```cpp
kahan_accumulator<double> acc;
acc += 1e16;
acc += 1.0;      // Would be lost with simple summation
acc += -1e16;
// acc() ≈ 1.0, not 0.0
```

Extra methods:
- `correction()` — current error compensation term

**`neumaier_accumulator<T>`**

Improved Kahan that handles the case where the running sum is smaller than the value being added.

```cpp
neumaier_accumulator<double> acc;
acc += 1.0;
acc += 1e16;     // Neumaier handles this correctly
acc += -1e16;
// acc() ≈ 1.0
```

**`klein_accumulator<T>`**

Second-order error compensation. Tracks corrections to the corrections.

```cpp
klein_accumulator<double> acc;
// For extreme precision requirements
```

Extra methods:
- `correction()` — first-order correction
- `second_correction()` — second-order correction

**`pairwise_accumulator<T, ChunkSize>`**

Accumulates into chunks, then sums chunks pairwise. Good cache locality and O(log n) error growth.

```cpp
pairwise_accumulator<double, 256> acc;  // 256-element chunks
for (double x : large_dataset) acc += x;
```

---

## Quadrature Rules

**Namespace:** `limes::algorithms::quadrature`

Quadrature rules define where to sample a function and how to weight each sample. All rules work on the reference interval [-1, 1] and are scaled by integrators.

### Interface

All rules satisfy the `QuadratureRule<T>` concept:

```cpp
rule.size()        // Number of nodes
rule.abscissa(i)   // i-th sample point ∈ [-1, 1]
rule.weight(i)     // i-th weight
```

### Available Rules

| Rule | Nodes | Exactness | Best For |
|------|-------|-----------|----------|
| `gauss_legendre<T, N>` | N | Polynomials up to degree 2N-1 | Smooth functions |
| `gauss_kronrod_15<T>` | 15 | Embedded 7-point for error | Adaptive integration |
| `clenshaw_curtis<T, N>` | N | Uses Chebyshev nodes | Functions with endpoint issues |
| `simpson_rule<T>` | 3 | Degree 3 polynomials | Simple, pedagogical |
| `trapezoidal_rule<T>` | 2 | Degree 1 polynomials | Simple, periodic functions |
| `midpoint_rule<T>` | 1 | Degree 1 polynomials | Simplest rule |
| `tanh_sinh_nodes<T>` | Variable | Double exponential | Endpoint singularities |

### Details

**`gauss_legendre<T, N>`**

Optimal for smooth functions. N points integrate polynomials up to degree 2N-1 exactly.

```cpp
gauss_legendre<double, 5> rule;  // 5-point rule
for (size_t i = 0; i < rule.size(); ++i) {
    double x = rule.abscissa(i);  // Sample point
    double w = rule.weight(i);    // Weight
}
```

Specialized implementations for N = 2, 3, 5 with precomputed high-precision nodes.

**`gauss_kronrod_15<T>`**

15-point rule with an embedded 7-point Gauss rule. The difference between them estimates error.

```cpp
gauss_kronrod_15<double> rule;

// Full 15-point result
double kronrod_sum = ...;

// Embedded 7-point result (uses indices 1,3,5,7,9,11,13)
double gauss_sum = ...;

// Error estimate
double error = std::abs(kronrod_sum - gauss_sum);
```

Extra members:
- `gauss_weights` — weights for embedded 7-point rule
- `gauss_indices` — which nodes belong to embedded rule

**`clenshaw_curtis<T, N>`**

Uses Chebyshev nodes (cosine-spaced). Good for functions that are well-approximated by polynomials. Includes endpoints.

```cpp
clenshaw_curtis<double, 17> rule;  // 17-point rule
```

**`simpson_rule<T>`**

Classic Simpson's 1/3 rule. 3 points at {-1, 0, 1} with weights {1/3, 4/3, 1/3}.

**`trapezoidal_rule<T>`**

2 points at endpoints. Exact for linear functions. Surprisingly good for periodic functions over a full period.

**`midpoint_rule<T>`**

Single point at center. Simplest possible rule.

**`tanh_sinh_nodes<T>`**

Double exponential (tanh-sinh) transformation. Generates nodes that cluster near endpoints, excellent for integrating functions with endpoint singularities.

```cpp
tanh_sinh_nodes<double> nodes;
double x = nodes.abscissa(level, index);
double w = nodes.weight(level, index);
```

Used internally by `tanh_sinh_integrator`.

---

## Integrators

**Namespace:** `limes::algorithms`

Integrators combine quadrature rules with accumulators to compute definite integrals.

### Result Type

All integrators return `integration_result<T>`:

```cpp
auto result = integrator(f, a, b, tol);

result.value()       // Computed integral value
result.error()       // Estimated error
result.converged()   // Whether tolerance was achieved
result.iterations()  // Number of iterations/subdivisions
result.evaluations() // Number of function evaluations
```

Results support arithmetic for combining sub-integrals:
```cpp
auto total = result1 + result2;    // Values and errors add
auto scaled = result * 2.0;        // Scale value and error
```

### Available Integrators

| Integrator | Method | Best For |
|------------|--------|----------|
| `quadrature_integrator` | Single rule or adaptive | General purpose |
| `romberg_integrator` | Richardson extrapolation | Smooth functions |
| `tanh_sinh_integrator` | Double exponential | Endpoint singularities, infinite intervals |

### Details

**`quadrature_integrator<T, Rule, Acc>`**

Combines any quadrature rule with any accumulator. Supports both single-pass and adaptive integration.

```cpp
// Compose rule and accumulator
using rule = quadrature::gauss_legendre<double, 5>;
using acc = accumulators::kahan_accumulator<double>;
quadrature_integrator<double, rule, acc> integrator{rule{}, acc{}};

// Single-pass (no error control)
auto r1 = integrator(f, 0.0, 1.0);

// Adaptive (subdivides until tolerance met)
auto r2 = integrator(f, 0.0, 1.0, 1e-10);
```

Adaptive mode uses recursive bisection with Richardson extrapolation.

**`romberg_integrator<T, Acc>`**

Romberg integration: repeatedly refines trapezoidal rule and extrapolates. Very efficient for smooth functions.

```cpp
romberg_integrator<double> integrator;
auto result = integrator(f, 0.0, 1.0, 1e-12);
```

Uses Kahan accumulator by default. Converges rapidly (often 5-10 iterations) for analytic functions.

**`tanh_sinh_integrator<T, Acc>`**

Double exponential quadrature. Transforms the integral so that nodes cluster near endpoints, handling singularities gracefully.

```cpp
tanh_sinh_integrator<double> integrator;

// Finite interval with endpoint singularity
auto r1 = integrator([](double x) { return 1/std::sqrt(x); }, 0.0, 1.0);

// Semi-infinite interval
auto r2 = integrator([](double x) { return std::exp(-x); }, 0.0, INFINITY);

// Doubly-infinite interval
auto r3 = integrator([](double x) { return std::exp(-x*x); }, -INFINITY, INFINITY);
```

Automatically detects and handles infinite bounds.

### Type Alias

```cpp
// Recommended adaptive configuration (Gauss-Kronrod + Neumaier)
adaptive_integrator<double> integrator;
```

---

## Usage Examples

### Basic Integration

```cpp
#include <limes/limes.hpp>

using namespace limes::algorithms;

auto f = [](double x) { return x * x; };
adaptive_integrator<double> integrator;
auto result = integrator(f, 0.0, 1.0, 1e-10);

std::cout << "∫x² dx from 0 to 1 = " << result.value() << "\n";  // ≈ 0.333333
```

### Custom Composition

```cpp
using namespace limes::algorithms;

// High-order Gauss rule with Klein accumulator for extreme precision
using rule = quadrature::gauss_legendre<double, 5>;
using acc = accumulators::klein_accumulator<double>;
quadrature_integrator<double, rule, acc> integrator{rule{}, acc{}};

auto result = integrator(f, a, b, 1e-14);
```

### Handling Singularities

```cpp
using namespace limes::algorithms;

// 1/√x has a singularity at x=0
auto f = [](double x) { return 1.0 / std::sqrt(x); };

tanh_sinh_integrator<double> integrator;
auto result = integrator(f, 0.0, 1.0, 1e-10);
// result.value() ≈ 2.0
```

### Infinite Intervals

```cpp
using namespace limes::algorithms;

// Gaussian integral: ∫exp(-x²) from -∞ to ∞ = √π
auto f = [](double x) { return std::exp(-x * x); };

tanh_sinh_integrator<double> integrator;
auto result = integrator(f, -INFINITY, INFINITY, 1e-10);
// result.value() ≈ 1.7724538509
```

---

## Choosing an Integrator

| Situation | Recommended |
|-----------|-------------|
| Smooth function, moderate precision | `adaptive_integrator<T>` |
| Smooth function, high precision | `romberg_integrator<T>` |
| Endpoint singularity (e.g., 1/√x) | `tanh_sinh_integrator<T>` |
| Infinite interval | `tanh_sinh_integrator<T>` |
| Known polynomial degree | `gauss_legendre<T, N>` with N > degree/2 |
| Periodic function over full period | `trapezoidal_rule` (surprisingly good) |

## Choosing an Accumulator

| Situation | Recommended |
|-----------|-------------|
| Speed critical, low precision OK | `simple_accumulator` |
| General purpose | `kahan_accumulator` |
| Mixed magnitude values | `neumaier_accumulator` |
| Extreme precision | `klein_accumulator` |
| Many terms (>1000) | `pairwise_accumulator` |
