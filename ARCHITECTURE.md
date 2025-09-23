# Algebraic Integrators - Architecture Overview

## Design Philosophy

This library embodies the Unix philosophy of "do one thing well" and emphasizes:
- **Simplicity**: Each component has a single, clear responsibility
- **Composability**: Components can be freely combined to create powerful integrators
- **Generic Programming**: Templates and concepts provide zero-cost abstractions
- **Performance**: Support for parallel execution and SIMD optimizations
- **Flexibility**: Type erasure and generic interfaces for maximum adaptability

## Core Architecture

### 1. Concepts (`include/concepts/`)

C++20 concepts define clear contracts for all template parameters:

- **Numeric Types**:
  - `Field`: Types supporting +, -, *, / operations
  - `RealField`: Field types with transcendental functions

- **Function Types**:
  - `UnivariateFunction`: Single-variable functions f(x) → y
  - `MultivariateFunction`: Multi-variable functions f(x₁, ..., xₙ) → y

- **Component Types**:
  - `Accumulator`: Types that can accumulate values with controlled precision
  - `QuadratureRule`: Types providing integration nodes and weights
  - `IntegrationResult`: Types holding integration results with metadata
  - `CoordinateTransform`: Types for change of variables

### 2. Accumulators (`include/accumulators/`)

Pluggable accumulation strategies for precision control:

```cpp
// Simple summation (fast, less accurate)
auto acc1 = make_simple<double>();

// Kahan-Babuška compensated summation
auto acc2 = make_kahan<double>();

// Neumaier summation (improved Kahan)
auto acc3 = make_neumaier<double>();

// Klein summation (second-order compensation)
auto acc4 = make_klein<double>();

// Pairwise summation (cache-friendly)
auto acc5 = make_pairwise<double>();
```

### 3. Quadrature Rules (`include/quadrature/`)

Orthogonal quadrature rule components:

- **Gauss-Legendre**: Optimal for smooth functions
- **Gauss-Kronrod**: Adaptive integration with embedded error estimation
- **Clenshaw-Curtis**: FFT-friendly, nested rules
- **Tanh-Sinh**: Double exponential for endpoint singularities
- **Simpson/Trapezoidal/Midpoint**: Classical rules

### 4. Integration Algorithms (`include/integrators/`)

Core integration algorithms that compose rules and accumulators:

```cpp
// Basic quadrature with any rule and accumulator
quadrature_integrator<T, Rule, Accumulator> integrator{rule, acc};

// Tanh-sinh for infinite intervals
tanh_sinh_integrator<T, Accumulator> de_integrator;

// Romberg with Richardson extrapolation
romberg_integrator<T, Accumulator> romberg;
```

### 5. Coordinate Transforms (`include/transforms/`)

Transforms for handling special integration domains:

- **Linear**: Affine transformations
- **Interval mapping**: Maps [a, b] ↔ [-1, 1]
- **Sigmoid/Tanh**: For infinite intervals
- **Log/Sinh**: For heavy-tailed integrands
- **IMT**: For endpoint singularities
- **Double exponential**: Tanh-sinh transform

### 6. Parallel Execution (`include/parallel/`)

Parallel integration with work-stealing and SIMD:

```cpp
// Parallel adaptive integration
parallel_integrator<T, BaseIntegrator> par_int{base, policy};

// Parallel Monte Carlo
parallel_monte_carlo<T> mc{n_samples, policy};
```

## Usage Patterns

### Simple Integration

```cpp
// Automatic method selection
auto result = integrate_adaptive(f, a, b, tolerance);
```

### Custom Configuration

```cpp
// Compose specific components
using acc = accumulators::klein_accumulator<double>;
using rule = quadrature::gauss_kronrod_15<double>;
quadrature_integrator<double, rule, acc> integrator{rule{}, acc{}};
auto result = integrator(f, a, b, tol);
```

### Builder Pattern

```cpp
auto integrator = make_integrator<double>()
    .with_quadrature("gauss-legendre")
    .with_accumulator("neumaier")
    .with_parallel(true)
    .with_threads(8)
    .with_tolerance(1e-12);

auto result = integrator.integrate(f, a, b);
```

### Transform-Based Integration

```cpp
// Handle singularity at x=0
auto transform = transforms::make_imt<double>(0.5);
auto result = integrate<double>::with_transform(f, 0, 1, transform);
```

## Performance Features

### SIMD Optimization

The library automatically uses AVX2/AVX-512 instructions when available for:
- Batch function evaluation
- Vector accumulation
- Parallel reduction

### Parallel Execution

- Work-stealing scheduler for adaptive subdivision
- Parallel Monte Carlo with independent RNG streams
- Configurable thread pools and chunk sizes

### Cache Optimization

- Pairwise summation for cache-friendly accumulation
- Blocked algorithms for large-scale integration
- Memory-aligned data structures for SIMD

## Extensibility

### Adding New Quadrature Rules

```cpp
template<concepts::Field T, std::size_t N>
struct my_rule : quadrature_rule<T, N> {
    constexpr my_rule() {
        // Initialize weights and abscissas
    }
};
```

### Adding New Accumulators

```cpp
template<concepts::Field T>
class my_accumulator {
public:
    using value_type = T;

    my_accumulator& operator+=(T value) {
        // Custom accumulation logic
        return *this;
    }

    T operator()() const { return /* accumulated value */; }
};
```

### Adding New Transforms

```cpp
template<concepts::Field T>
class my_transform : coordinate_transform<T> {
public:
    T forward(T x) const override { /* x → y */ }
    T jacobian(T x) const override { /* dy/dx */ }
    T inverse(T y) const override { /* y → x */ }
};
```

## Best Practices

1. **Start Simple**: Use `integrate_adaptive()` for most cases
2. **Profile First**: Only optimize when measurements show need
3. **Choose Accumulator Wisely**:
   - Simple for speed when precision isn't critical
   - Kahan/Neumaier for general high precision
   - Klein for extreme precision requirements
4. **Parallelize Appropriately**: Only for expensive integrands or large domains
5. **Use Transforms**: For singularities and infinite intervals

## Design Rationale

### Why Header-Only?

- Zero-cost abstractions through aggressive inlining
- No ABI compatibility issues
- Easy integration (just include headers)
- Template instantiation flexibility

### Why Concepts?

- Clear, documented requirements for template parameters
- Better error messages
- Enables composition without inheritance
- Supports both static and dynamic polymorphism

### Why Orthogonal Components?

- Mix and match algorithms, rules, and accumulators
- Test components in isolation
- Reuse components across different contexts
- Clear separation of concerns

## Future Directions

- Multivariate integration with sparse grids
- Oscillatory integrals with specialized methods
- GPU acceleration via CUDA/SYCL
- Automatic differentiation integration
- Interval arithmetic for guaranteed bounds