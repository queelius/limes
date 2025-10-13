# Integration API

## High-Level Functions

### `integrate_adaptive`

Adaptive integration with automatic error control.

```cpp
template <typename F, typename T>
auto integrate_adaptive(F f, T a, T b, T tolerance);
```

**Parameters:**
- `f` - Function to integrate
- `a`, `b` - Integration bounds
- `tolerance` - Desired accuracy

**Returns:** `integration_result<T>` with `.value` and `.error`

**Example:**
```cpp
auto f = [](double x) { return std::sin(x); };
auto result = integrate_adaptive(f, 0.0, M_PI, 1e-10);
// result.value ≈ 2.0
```

## Integrator Classes

### `gauss_kronrod_integrator`

High-accuracy adaptive integration.

```cpp
template <typename Acc, int N, int M>
class gauss_kronrod_integrator;
```

**Example:**
```cpp
using Acc = kahan_accumulator<double>;
gauss_kronrod_integrator<Acc, 15, 31> integrator;
auto result = integrator([](double x) { return x*x; }, 0.0, 1.0, 1e-8);
```

### `tanh_sinh_integrator`

Best for functions with endpoint singularities.

```cpp
using Acc = neumaier_accumulator<double>;
tanh_sinh_integrator<Acc> integrator;
auto result = integrator([](double x) { return 1/sqrt(x); }, 0.0, 1.0, 1e-8);
```

### `simpson_univariate_integrator`

Simple Simpson's rule integration.

```cpp
simpson_univariate_integrator<Acc> integrator;
auto result = integrator(f, a, b, n_intervals);
```

## Transforms

### Infinite Intervals

```cpp
// Integrate from 0 to ∞
auto result = integrate_transform(
    [](double x) { return std::exp(-x); },
    0.0, std::numeric_limits<double>::infinity(),
    1e-8
);
```

### Singularities

```cpp
// Handle singularity at x=0
auto result = integrate_tanh_sinh(
    [](double x) { return 1/std::sqrt(x); },
    0.0, 1.0, 1e-8
);
```

## Result Type

```cpp
template <typename T>
struct integration_result {
    T value;           // Computed integral
    T error;           // Error estimate
    size_t evaluations; // Function evaluations
};
```
