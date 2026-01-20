# Motivation {#motivation}

## The Question That Started It All

> *Can we build something for integration analogous to what autograd does for differentiation?*

This question launched limes. Automatic differentiation (autograd) revolutionized machine learning by computing gradients automatically through expression graphs. Could we do the same for integrals?

## The Honest Answer: No, But...

The honest answer is **no**—not in the same way.

**Differentiation is local.** The chain rule decomposes perfectly: to differentiate `f(g(x))`, you only need `f'(g(x))` and `g'(x)` at a single point. The computation is algebraic.

**Integration is global.** There's no "chain rule" for integrals. To compute an integral, you fundamentally need to visit infinitely many points (approximated by quadrature). The computation is analytic.

This asymmetry means we can't build an "autointegrate" that automatically computes integrals by understanding expression structure the way autograd computes gradients.

## What We CAN Build

Instead of mimicking autograd's *mechanism*, we take inspiration from its *design philosophy*: **expressions as first-class composable objects**.

The insight is:

> **Integrals are algebraic objects that compose according to mathematical laws.**

These laws are well-known:
- **Linearity**: `∫(af + bg) = a∫f + b∫g`
- **Domain additivity**: `∫[a,c] = ∫[a,b] + ∫[b,c]`
- **Fubini's theorem**: `∫∫f dxdy = ∫∫f dydx` (when valid)
- **Separability**: `∫∫f(x)g(y) dxdy = (∫f dx)(∫g dy)`

limes makes these laws **explicit and exploitable** in the API.

## The Black-Box Problem

Most numerical integration libraries treat functions as black boxes:

```cpp
// Traditional approach
double f(double x) { return std::sin(x * x); }
double result = integrate(f, 0.0, 1.0);
```

The integrator has no idea what `f` does. It can only:
1. Evaluate `f` at points
2. Use generic adaptive strategies

This works, but it's leaving information on the table.

## The Expression Approach

limes builds expressions from primitives, making structure visible:

```cpp
// limes approach
auto x = arg<0>;
auto f = sin(x * x);
auto I = integral(f).over<0>(0.0, 1.0);
```

Now the system *knows* `f` is `sin(x²)`. This enables:

### 1. Symbolic Differentiation

Since `f = sin(x²)` is built from known primitives, we can differentiate symbolically:

```cpp
auto df = derivative(f).wrt<0>();  // 2x·cos(x²)
```

The chain rule is applied at compile time. No numerical differentiation, no epsilon parameters, no round-off error accumulation.

### 2. Separability Detection

Consider this double integral:

```cpp
auto f = exp(-x*x) * exp(-y*y);  // Gaussian in x times Gaussian in y
auto I = integral(f).over<0>(-10, 10).over<1>(-10, 10);
```

A black-box integrator would use nested quadrature: for each `y` point, integrate over `x`. That's O(n²) evaluations.

But limes can detect that `f(x,y) = g(x)·h(y)` where `g(x) = exp(-x²)` and `h(y) = exp(-y²)`. The integral factors:

```cpp
∫∫ f(x,y) dx dy = (∫ g(x) dx)(∫ h(y) dy)
```

Now we need only O(n) evaluations for each 1D integral.

### 3. Domain Intelligence

Integrals carry domain information that can be manipulated:

```cpp
auto I = integral(f).over<0>(0.0, 1.0);

// Split at a discontinuity
auto [left, right] = I.split(0.5);

// Change of variables to remove a singularity
auto T = I.transform(
    [](double t) { return t*t; },     // x = t²
    [](double t) { return 2*t; },     // |dx/dt| = 2t
    0.0, 1.0
);
```

### 4. Method Selection

Different integrands need different methods. With expressions, we can (eventually) make intelligent choices:

```cpp
// Smooth periodic function → trapezoidal rule excels
// Endpoint singularity → tanh-sinh transformation
// High dimension → Monte Carlo
// Polynomial → Gauss with appropriate order
```

## Design Principles

limes follows these principles:

### 1. Separation of What and How

The *expression* (what to compute) is separate from the *method* (how to compute it):

```cpp
auto I = integral(f).over<0>(0.0, 1.0);  // What

I.eval();              // How: default method
I.eval(gauss<7>());    // How: 7-point Gauss
I.eval(adaptive(1e-12));  // How: adaptive
```

### 2. Composability Over Features

Simple pieces that compose well beat monolithic features:

```cpp
// Methods compose
auto method = make_adaptive(gauss<5>(), 1e-10);  // Adaptive Gauss

// Integrals compose
auto IJ = I * J;  // Product of independent integrals

// Expressions compose
auto df = derivative(f).wrt<0>();
auto I_df = integral(df).over<0>(a, b);  // ∫f'(x)dx = f(b) - f(a)
```

### 3. Concepts Over Inheritance

Extension through concepts enables generic programming:

```cpp
template<typename M>
concept IntegrationMethod = requires(M m, std::function<T(T)> f, T a, T b) {
    { m(f, a, b) } -> std::convertible_to<integration_result<T>>;
};
```

Any type satisfying the concept works with the system.

### 4. Inspectability

Expressions are trees you can examine:

```cpp
std::cout << I.to_string();
// Output: (integral (sin (* x x)) [0 0.0 1.0])
```

This makes debugging possible and the system pedagogically valuable.

## Stepanov's Influence

limes draws from Alex Stepanov's approach to generic programming:

1. **Identify algebraic structure** — Integrals form a vector space over ℝ
2. **Express it in the type system** — Integral expressions are types with algebraic operations
3. **Program generically against concepts** — Methods and accumulators are pluggable

The goal is not just to compute integrals, but to express the *mathematics of integration* in C++.

## What limes Is Not

- **Not a computer algebra system** — We don't solve integrals symbolically in general
- **Not a replacement for numerical libraries** — The algorithms layer wraps standard techniques
- **Not the fastest possible implementation** — Clarity and composability come first

## What limes Is

- A **pedagogical library** demonstrating expression-based numeric computing
- A **composable API** where mathematical laws are first-class citizens
- A **testbed** for exploring how much structure we can exploit
- A **proof of concept** that expressions > black boxes for many use cases

## Summary

The key insight of limes:

> While we cannot automate integration like autograd automates differentiation, we CAN build a system where both differentiation and integration are **first-class algebraic objects** with rich structure that the computer can exploit.

Expressions beat black boxes when structure matters.
