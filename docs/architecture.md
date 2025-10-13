# Architecture

## Design Philosophy

The library follows a **composable architecture** where components can be mixed and matched:

- **Accumulators** control numerical precision
- **Quadrature Rules** define integration methods
- **Integrators** combine rules with accumulators
- **Transforms** handle special cases (infinite intervals, singularities)

## Core Components

### 1. Accumulators

Control how floating-point summation is performed:

- `simple_accumulator` - Standard summation
- `kahan_accumulator` - Compensated summation
- `neumaier_accumulator` - Improved Kahan
- `klein_accumulator` - Second-order compensation
- `pairwise_accumulator` - Divide-and-conquer summation

### 2. Quadrature Rules

Integration methods:

- `gauss_legendre_rule` - Gaussian quadrature
- `gauss_kronrod_rule` - Adaptive Gauss-Kronrod
- `tanh_sinh_rule` - Double exponential (best for singularities)
- `simpson_rule` - Simpson's method

### 3. Integrators

Combine accumulators with rules:

```cpp
// Custom composition
using Acc = kahan_accumulator<double>;
using Rule = gauss_kronrod_rule<Acc, 15, 31>;
gauss_kronrod_integrator<Acc, 15, 31> integrator;
```

### 4. ODE Solvers

- `euler_ode1` - First-order Euler method
- `rk4_ode1` - Runge-Kutta 4th order
- Second-order variants for systems

### 5. Differentiation

- `central_finite_difference` - 8th order accurate
- `central_finite_difference_2nd` - Second derivatives
- `grad` - Multivariable gradients

## Type Safety

The library uses C++20 concepts for compile-time interface checking:

```cpp
template <typename T>
concept accumulator = requires(T a, typename T::value_type v) {
    { a.add(v) } -> std::same_as<void>;
    { a.result() } -> std::same_as<typename T::value_type>;
};
```

## Performance

- **Header-only** - Zero runtime overhead
- **Template specialization** - Compiler optimizations
- **SIMD support** - Optional AVX2 acceleration
- **Parallelization** - Optional OpenMP support
