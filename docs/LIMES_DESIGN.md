# limes: Composable Calculus Expressions

*limes* (Latin: "limit", "boundary") — a C++20 library for composable calculus expressions.

Limits are the foundation of calculus: derivatives and integrals are both defined through limits. This library makes calculus operations first-class composable objects with symbolic computation where possible and robust numerical methods as fallback.

## 1. Motivation

### The Original Question

Can we build something for integration analogous to what autograd does for differentiation?

### The Honest Answer

No—not in the same way. Differentiation is **local** (chain rule decomposes perfectly), while integration is **global** (no analogous decomposition exists). We cannot build an "autointegrate" that automatically computes integrals by understanding expression structure the way autograd computes gradients.

### What We CAN Build

Instead of mimicking autograd's mechanism, we take inspiration from its **design philosophy**: expressions as first-class composable objects with deferred evaluation. The insight is:

> **Integrals are algebraic objects that compose according to mathematical laws.**

This library makes those laws explicit in the API.

## 2. Core Vision

### Integrals as First-Class Expressions

An `Expr<T>` is a **lazy**, **composable**, **inspectable** representation of an integration problem. It captures:
- The integrand (what to integrate)
- The domain (where to integrate)
- The structure (how integrals compose)

Evaluation is separate from construction—build complex expressions, then evaluate with chosen methods.

### Pedagogical Goals

Users should learn:
1. **Integration concepts**: bounds, domains, iteration, Fubini's theorem
2. **Numerical methods**: quadrature rules, error estimation, convergence
3. **Generic programming**: composable abstractions with concepts and policies

### Design Principles

1. **Expressiveness**: Complex integration problems should be writable clearly
2. **Composability**: Expressions combine via algebraic laws (linearity, additivity, nesting)
3. **Inspectability**: Structure is visible—print expressions, debug construction
4. **Separation**: The *what* (expression) is separate from the *how* (numerical method)

## 3. API Overview

### Building Integrands

Integrands are built from expression primitives, not passed as black-box functions:

```cpp
using namespace limes;

// Variables (positional)
auto x = arg<0>;
auto y = arg<1>;

// Build expressions from primitives
auto f = sin(x) * exp(-y);           // sin(x) · e^(-y)
auto g = x*x + y*y;                  // x² + y²
auto h = sqrt(1 - x*x);              // √(1 - x²)
```

### Constructing Integrals

```cpp
// Basic integral: ∫₀¹ x² dx
auto I = integral(x*x).over(0, 0.0, 1.0);

// Double integral: ∫₀¹∫₀ˣ sin(x)·e^(-y) dy dx
auto J = integral(sin(x) * exp(-y))
    .over(1, 0.0, x)                   // integrate y from 0 to x
    .over(0, 0.0, 1.0);                // integrate x from 0 to 1

// Marginal: integrate out y
auto joint = exp(-(x*x + y*y) / 2);
auto marginal = integral(joint).over(1, -10.0, 10.0);
```

### Composing Expressions

```cpp
// Linearity: ∫(af + bg) = a∫f + b∫g
auto sum = 2.0 * I + 3.0 * J;

// Domain splitting: ∫₀¹ = ∫₀^{0.5} + ∫_{0.5}^1
auto [left, right] = I.split(0, 0.5);

// Change of variables
auto transformed = I.transform(0, phi, phi_jacobian);
```

### Differentiation

Symbolic differentiation via chain rule:

```cpp
auto f = sin(x*x);                    // sin(x²)
auto df = derivative(f, 0);           // d/dx[sin(x²)] = 2x·cos(x²)

auto g = exp(x) * sin(y);
auto grad_g = gradient(g);            // [∂g/∂x, ∂g/∂y]

// Fundamental Theorem: ∫(df/dx)dx = f
auto F = integral(df).over(0, 0.0, x);  // Should equal sin(x²) - sin(0)
```

### Evaluating Expressions

```cpp
// Evaluate with default method
auto result = I.evaluate();

// Evaluate with specific method (template parameter)
auto result = I.evaluate<AdaptiveGaussKronrod>();

// Result contains value, error estimate, convergence info
double value = result.value();
double error = result.error();
bool ok = result.converged();
```

### Inspecting Expressions

```cpp
// Print structure (S-expression format)
std::cout << I.to_string() << "\n";
// Output: (Integral f [0 0.0 1.0])

std::cout << J.to_string() << "\n";
// Output: (Integral (Integral f [1 0.0 (Bound x)]) [0 0.0 1.0])
```

## 4. Core Types

### `Expr<T>`

The central type—a lazy integral expression.

Properties:
- **Value type**: `T` (typically `double` or `float`)
- **Arity**: Number of free (unintegrated) variables, tracked at runtime
- **Structure**: Tree of expression nodes

Key operations:
- `.over(dim, lo, hi)` — bind dimension `dim` to integration with bounds `[lo, hi]`
- `.evaluate<Method>()` — compute numerical result
- `.to_string()` — structural representation

### `Result<T>`

The result of evaluation.

Properties:
- `.value()` — the computed integral value
- `.error()` — estimated error (method-dependent)
- `.converged()` — whether tolerance was achieved
- `.iterations()` — number of function evaluations or subdivisions

### Bounds

Bounds can be:
- **Constant**: `0.0`, `1.0`, `inf`
- **Callable**: `[](auto x){ return x*x; }` — depends on outer variables

## 5. Composition Model

### Algebraic Laws

The API exposes integration's algebraic structure:

| Law | Expression | API |
|-----|------------|-----|
| Linearity | ∫(af+bg) = a∫f + b∫g | `a*I + b*J` |
| Additivity | ∫ₐᵇ + ∫ᵇᶜ = ∫ₐᶜ | `left + right` after `.split()` |
| Separability | ∫∫f(x)g(y) dxdy = (∫f dx)(∫g dy) | Future: auto-detected with expression nodes |
| Fubini | ∫∫f dxdy = ∫∫f dydx (when valid) | `.swap(0, 1)` |

### Variable Model

Variables are **positional** (indices, not names):
- `arg<0>` is the first variable
- `arg<1>` is the second, etc.

When `.over(dim, lo, hi)` is called:
- Dimension `dim` becomes the innermost integration variable
- Expression arity decreases by 1
- Remaining variables shift down if needed

### Nesting Semantics

Chained `.over()` calls create **iterated integrals**, evaluated inside-out:

```cpp
integral(f).over(1, a, b).over(0, c, d)
```

Means: for each x in [c,d], compute ∫ₐᵇ f(x,y) dy, then integrate that over x.

## 6. Expression Nodes

All integrands are built from expression nodes—there are no opaque black-box functions. This enables structural analysis, symbolic antiderivatives, and pattern detection (like separability).

### Core Nodes

| Node | Meaning |
|------|---------|
| `Const` | A constant value |
| `Var` | A positional variable reference `arg<i>` |
| `BinaryOp` | +, -, *, /, ^ (power) |
| `UnaryOp` | Negation, etc. |
| `Integral` | Integration node with integrand + domain |
| `Bound` | Bound expression (constant or callable) |

### Primitive Functions (v1)

| Node | Meaning | Derivative | Antiderivative |
|------|---------|------------|----------------|
| `Exp` | e^x | d/dx[e^x] = e^x | ∫e^x = e^x |
| `Log` | ln(x) | d/dx[ln(x)] = 1/x | ∫ln(x) = x·ln(x) - x |
| `Sin` | sin(x) | d/dx[sin(x)] = cos(x) | ∫sin(x) = -cos(x) |
| `Cos` | cos(x) | d/dx[cos(x)] = -sin(x) | ∫cos(x) = sin(x) |
| `Sqrt` | √x | d/dx[√x] = 1/(2√x) | ∫√x = (2/3)x^(3/2) |
| `Abs` | |x| | d/dx[|x|] = sign(x) | ∫|x| = (x·|x|)/2 |

Each primitive knows both its derivative and antiderivative. The chain rule composes derivatives automatically; when an integral has a known antiderivative, numerical integration is avoided.

### Extensibility via Concepts

Users can define custom nodes by modeling the `ExprNode` concept:

```cpp
template<typename N, typename T>
concept ExprNode = requires(N node, std::span<T> args) {
    { node.evaluate(args) } -> std::convertible_to<T>;
    { node.arity() } -> std::convertible_to<std::size_t>;
    { node.to_string() } -> std::convertible_to<std::string>;
};

// Optional: nodes that know their antiderivative
template<typename N, typename T, std::size_t Dim>
concept IntegrableNode = ExprNode<N, T> && requires(N node) {
    { node.template antiderivative<Dim>() };  // returns an expression
};
```

Custom nodes that model these concepts work seamlessly with the system.

## 7. Domain Operations

### Split

Divide a domain at a point:

```cpp
auto I = integral(f).over(0, 0.0, 1.0);
auto [left, right] = I.split(0, 0.5);
// left: ∫₀^{0.5} f dx
// right: ∫_{0.5}^1 f dx
// Invariant: left + right ≈ I
```

### Transform (Change of Variables)

Apply a substitution x = φ(t):

```cpp
// Transform x ∈ [0, ∞) to t ∈ [0, 1) via x = t/(1-t)
auto transformed = I.transform(0,
    [](auto t){ return t / (1 - t); },        // φ(t)
    [](auto t){ return 1 / ((1-t)*(1-t)); }   // |φ'(t)|
);
```

## 8. Evaluation

### Method Selection

Numerical methods are specified as template parameters:

```cpp
auto result = expr.evaluate<AdaptiveGaussKronrod>();
auto result = expr.evaluate<Romberg>();
auto result = expr.evaluate<MonteCarlo>();
```

The default `evaluate()` uses a sensible adaptive method.

### Error Estimation

Methods that support error estimation populate `result.error()`. Methods without inherent error estimation may return an estimate of 0 or use heuristics.

### Inside-Out Evaluation

For nested integrals, evaluation proceeds inside-out:
1. Inner integral becomes a function of outer variables
2. Outer integration evaluates that function at quadrature points
3. Repeat outward

## 9. Namespace Organization

`limes` is the top-level namespace. The numerical algorithms (currently in `calckit`) become a subnamespace:

```
limes::                       # Top-level: expression API
    integral()                  # Expression construction
    Expr<T>                  # Core expression type
    Result<T>                # Evaluation result

limes::algorithms::           # Numerical backend
    accumulators::              # Precision control (Kahan, Neumaier, etc.)
    quadrature::                # Node/weight generation
    integrators::               # Direct numerical methods
```

The expression layer **uses** the algorithms internally. Users can:
- Use `limes::integral()` for composable expressions (recommended)
- Use `limes::algorithms::` directly for low-level control

## 10. Scope and Phasing

### v1: Core Expression System

**Expression-only design** — all integrands built from nodes, no black-box functions.

- **Nodes**: Const, Var, BinaryOp, UnaryOp, Integral, Bound
- **Primitives**: exp, log, sin, cos, sqrt, abs (with derivatives and antiderivatives)
- **Differentiation**: `derivative(expr, dim)`, `gradient(expr)` — symbolic via chain rule
- **Integration**: `integral(expr).over(dim, lo, hi)` — symbolic when possible, numerical otherwise
- **Multivariate**: N-dimensional, chained `.over()` for iterated integrals
- **Dependent bounds**: Lambda bounds for non-rectangular regions
- **Arithmetic composition**: +, -, scalar *, /
- **Concepts**: `ExprNode`, `DifferentiableNode`, `IntegrableNode` for extensibility
- **Evaluation**: Adaptive Gauss-Kronrod with `Result<T>` (value, error, convergence)
- **Inspection**: S-expression printing

### v1.1: Domain Operations

- **Split**: Divide domain at a point
- **Transform**: Change of variables with Jacobian
- **Separability detection**: Auto-factor ∫∫f(x)g(y) → (∫f)(∫g)

### v2+: Extended Features

**Extended nodes:**
- `Sum`, `Product` — finite Σ and Π
- `Derivative` — for mixed integration/differentiation
- `Conditional` — piecewise functions
- `Limit`, `Convolution`

**Numerical methods:**
- Additional methods (tanh-sinh, Monte Carlo, sparse grids)
- Automatic method selection based on expression structure

**Symbolic integration techniques:**
- U-substitution — algorithmic pattern matching: recognize ∫f(g(x))·g'(x)dx
- Integration by parts — ∫u·dv = uv - ∫v·du
- Partial fractions — for rational functions
- Trigonometric identities — sin²+cos²=1, etc.

**Advanced features:**
- Parallelization of multi-dimensional integrals
- Integration with automatic differentiation
- Caching/memoization of repeated subexpressions

## 11. Example: Bayesian Marginal

A motivating use case—computing a marginal distribution:

```cpp
using namespace limes;

// Variables
auto x = arg<0>;
auto y = arg<1>;

// Joint PDF: p(x, y) = exp(-(x² + y²)/2) / 2π
auto joint_pdf = exp(-(x*x + y*y) / 2) / (2 * pi);

// Marginal: p(x) = ∫ p(x,y) dy
auto marginal_expr = integral(joint_pdf).over(1, -10.0, 10.0);

// Evaluate for specific x values
for (double xval : {-2.0, -1.0, 0.0, 1.0, 2.0}) {
    auto p_x = marginal_expr.bind(0, xval).evaluate();
    std::cout << "p(" << xval << ") = " << p_x.value() << "\n";
}

// Integrate again for normalization check
auto total = integral(marginal_expr).over(0, -10.0, 10.0);
std::cout << "∫p(x)dx = " << total.evaluate().value() << "\n";  // ≈ 1
```

### Why Expression-Based?

Because `joint_pdf` is built from primitives (`exp`, `*`, `/`), the system can:
1. Detect that exp(-(x² + y²)/2) = exp(-x²/2) · exp(-y²/2) — **separable!**
2. Factor the integral: (∫exp(-x²/2)dx)(∫exp(-y²/2)dy)
3. Use symbolic antiderivatives where possible
4. Choose optimal numerical methods for what remains

## 12. Summary

limes provides:

1. **A novel API** for calculus as composable expressions — both differentiation and integration
2. **Symbolic computation** where possible — chain rule for derivatives, known antiderivatives for integrals
3. **Numerical fallback** — robust quadrature when symbolic methods don't apply
4. **Pedagogical clarity** through inspectable expression trees and explicit algebraic laws
5. **Extensibility** via concepts — define custom nodes that "just work"

The key insight: while we can't fully automate integration like autograd automates differentiation, we CAN build a unified expression system where both operations are first-class citizens with rich algebraic structure and Stepanov-style generic programming.
