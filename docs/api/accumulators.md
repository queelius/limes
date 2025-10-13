# Accumulators API

Accumulators control how floating-point summation is performed, trading speed for accuracy.

## Available Accumulators

### `simple_accumulator<T>`

Standard floating-point summation.

```cpp
simple_accumulator<double> acc;
acc.add(1.0);
acc.add(2.0);
double sum = acc.result(); // 3.0
```

**Pros:** Fastest
**Cons:** Accumulates roundoff error
**Use when:** Speed is critical, high precision not required

### `kahan_accumulator<T>`

Compensated summation (Kahan algorithm).

```cpp
kahan_accumulator<double> acc;
for (int i = 0; i < 1000000; i++) {
    acc.add(0.1);
}
// More accurate than simple accumulator
```

**Pros:** Much better precision than simple
**Cons:** Slightly slower
**Use when:** Summing many small values

### `neumaier_accumulator<T>`

Improved Kahan algorithm.

```cpp
neumaier_accumulator<double> acc;
```

**Pros:** More robust than Kahan for mixed-magnitude values
**Cons:** Slightly more expensive
**Use when:** Values vary widely in magnitude

### `klein_accumulator<T>`

Second-order compensated summation.

```cpp
klein_accumulator<double> acc;
```

**Pros:** Best accuracy of compensated algorithms
**Cons:** More expensive than Kahan/Neumaier
**Use when:** Maximum precision needed, performance acceptable

### `pairwise_accumulator<T>`

Divide-and-conquer summation.

```cpp
pairwise_accumulator<double> acc;
```

**Pros:** Better balance of speed and accuracy
**Cons:** Uses more memory for internal tree
**Use when:** Good general-purpose choice

## Accumulator Interface

All accumulators satisfy the `accumulator` concept:

```cpp
template <typename T>
concept accumulator = requires(T a, typename T::value_type v) {
    typename T::value_type;
    { a.add(v) } -> std::same_as<void>;
    { a.result() } -> std::same_as<typename T::value_type>;
};
```

## Choosing an Accumulator

```cpp
// Quick guide:
simple_accumulator      // Speed-critical, low precision OK
kahan_accumulator       // Good default choice
neumaier_accumulator    // Mixed magnitudes
klein_accumulator       // Maximum precision
pairwise_accumulator    // Balance of speed/accuracy
```

## Usage in Integrators

```cpp
// Compose with any integrator
using Acc = kahan_accumulator<double>;
gauss_kronrod_integrator<Acc, 15, 31> integrator;

// Or use with ODE solvers
rk4_ode1<klein_accumulator<double>> solver;

// Simpson's rule with neumaier
simpson_univariate_integrator<neumaier_accumulator<float>> integrator;
```

## Performance Comparison

| Accumulator | Speed | Accuracy | Memory |
|-------------|-------|----------|--------|
| Simple      | ⭐⭐⭐⭐⭐ | ⭐ | ⭐⭐⭐⭐⭐ |
| Kahan       | ⭐⭐⭐⭐ | ⭐⭐⭐⭐ | ⭐⭐⭐⭐⭐ |
| Neumaier    | ⭐⭐⭐⭐ | ⭐⭐⭐⭐ | ⭐⭐⭐⭐⭐ |
| Klein       | ⭐⭐⭐ | ⭐⭐⭐⭐⭐ | ⭐⭐⭐⭐ |
| Pairwise    | ⭐⭐⭐⭐ | ⭐⭐⭐⭐ | ⭐⭐⭐ |
