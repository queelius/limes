# limes: Composable Calculus Expressions {#mainpage}

**limes** (Latin: *limit*, *boundary*) is a C++20 header-only library for composable calculus expressions with symbolic differentiation and numerical integration.

## Quick Example

```cpp
#include <limes/limes.hpp>
using namespace limes::expr;

auto x = arg<0>;
auto f = sin(x * x);                       // f(x) = sin(x^2)
auto df = derivative(f).wrt<0>();          // df/dx = 2x*cos(x^2)

auto I = integral(f).over<0>(0.0, 1.0);    // Integral from 0 to 1
auto result = I.eval();                    // Approx 0.3103
```

## Why limes?

Most numerical libraries treat functions as black boxes. limes treats mathematical expressions as **algebraic objects** that compose according to mathematical laws:

- **Derivatives** are computed symbolically via chain rule at compile time
- **Integrals** compose via linearity, Fubini's theorem, and separability
- **Methods** are first-class objects you can mix and match

This design enables optimizations impossible with black-box functions:

```cpp
// limes detects these are independent and evaluates them separately
auto I = integral(sin(x)).over<0>(0.0, pi);   // depends on x
auto J = integral(exp(y)).over<1>(0.0, 1.0);  // depends on y
auto IJ = I * J;                               // ProductIntegral: (I)(J), not nested
```

## Architecture

| Layer | Namespace | Purpose |
|-------|-----------|---------|
| **Expression** | `limes::expr` | User-facing API for composable calculus |
| **Methods** | `limes::methods` | Integration method objects |
| **Algorithms** | `limes::algorithms` | Low-level numerical backend |

## Documentation

### Guides

- @subpage motivation "Motivation" - The problem limes solves and design philosophy
- @subpage tutorial "Tutorial" - Step-by-step introduction
- @subpage examples "Examples" - Complete worked examples

### Reference

- @ref expr_nodes "Expression Nodes" - Building blocks: Var, Const, Binary, Unary
- @ref expr_ops "Expression Operations" - Differentiation, integration, analysis
- @ref methods "Integration Methods" - gauss, monte_carlo, adaptive
- **[Algorithms](algorithms.md)** - Low-level integrators and accumulators

### Design

- **[Design Document](LIMES_DESIGN.md)** - Full API design and philosophy

## Key Features

### Symbolic Differentiation

Chain rule applied at compile time:

```cpp
auto f = exp(sin(x));
auto df = derivative(f).wrt<0>();      // cos(x)*exp(sin(x))
auto d2f = derivative(f).wrt<0, 0>();  // Second derivative
auto grad = derivative(f).gradient();  // All partials as tuple
```

See @ref limes::expr::derivative and @ref limes::expr::DerivativeBuilder.

### Numerical Integration

Fluent builder with method selection:

```cpp
auto I = integral(x*x).over<0>(0.0, 1.0);

I.eval();                    // Default (adaptive)
I.eval(gauss<7>());          // 7-point Gauss-Legendre
I.eval(monte_carlo(10000));  // Monte Carlo
I.eval(adaptive(1e-12));     // Adaptive with tolerance
```

See @ref limes::expr::Integral and @ref limes::expr::IntegralBuilder.

### Box Integration

Monte Carlo over N-dimensional rectangular regions:

```cpp
auto I = integral(x*y*z).over_box({{0,1}, {0,1}, {0,1}});
auto result = I.eval(monte_carlo(100000));
```

See @ref limes::expr::BoxIntegral.

### Separable Composition

Multiply independent integrals for automatic factorization:

```cpp
auto I = integral(sin(x)).over<0>(0, pi);   // depends on x
auto J = integral(exp(y)).over<1>(0, 1);    // depends on y
auto IJ = I * J;                             // ProductIntegral
```

See @ref limes::expr::ProductIntegral.

### Integration Methods

First-class method objects with builder configuration:

```cpp
using namespace limes::methods;

gauss<7>()                           // 7-point Gauss-Legendre
monte_carlo(10000).with_seed(42)     // Reproducible Monte Carlo
adaptive().with_tolerance(1e-12)     // Custom tolerance
make_adaptive(gauss<5>(), 1e-10)     // Adaptive wrapper
```

See @ref limes::methods.

## Getting Started

Include the main header:

```cpp
#include <limes/limes.hpp>
```

Or include specific modules:

```cpp
#include <limes/expr/expr.hpp>           // Expression layer
#include <limes/methods/methods.hpp>     // Method objects
#include <limes/algorithms/algorithms.hpp>  // Low-level algorithms
```

## Installation

limes is header-only. Copy `include/limes` to your project or use CMake:

```cmake
find_package(limes REQUIRED)
target_link_libraries(your_target PRIVATE limes::limes)
```

## Requirements

- C++20 compiler (GCC 10+, Clang 12+, MSVC 19.29+)
- CMake 3.16+
