# limes — Composable Calculus Expressions

*limes* (Latin: "limit", "boundary") — a C++20 header-only library for composable calculus expressions with symbolic differentiation and numerical integration.

[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)
[![C++20](https://img.shields.io/badge/C%2B%2B-20-blue.svg)](https://en.cppreference.com/w/cpp/20)

## What is limes?

limes makes calculus operations **first-class composable objects**. Build mathematical expressions, differentiate them symbolically via chain rule, and integrate them numerically with your choice of method—all with a clean, fluent API.

```cpp
#include <limes/limes.hpp>
using namespace limes::expr;

auto x = arg<0>;
auto f = sin(x * x);                           // f(x) = sin(x²)
auto df = derivative(f).wrt<0>();              // df/dx = 2x·cos(x²)

auto I = integral(f).over<0>(0.0, 1.0);        // ∫₀¹ sin(x²) dx
auto result = I.eval();                        // ≈ 0.3103
```

## Features

### Expression Layer (`limes::expr`)
- **Symbolic differentiation** via chain rule — derivatives computed at compile time
- **Fluent builder API** — `derivative(f).wrt<0>()`, `integral(f).over<0>(a, b)`
- **Named variables** — `var(0, "x")` for readable `to_string()` output
- **N-dimensional integration** — nested integrals with dependent bounds
- **Box integration** — Monte Carlo over rectangular regions
- **Separable composition** — `I * J` for independent integrals

### Algorithm Layer (`limes::algorithms`)
- **Multiple quadrature rules** — Gauss-Legendre, Gauss-Kronrod, Clenshaw-Curtis, Tanh-Sinh
- **Adaptive integration** — automatic interval subdivision with error estimation
- **High-precision accumulators** — Kahan, Neumaier, Klein compensated summation
- **Integration methods as objects** — `gauss<7>()`, `monte_carlo(10000)`, `adaptive()`

## Quick Start

### Basic Integration

```cpp
using namespace limes::expr;

auto x = arg<0>;
auto I = integral(x * x).over<0>(0.0, 1.0);    // ∫₀¹ x² dx
auto result = I.eval();                         // = 1/3
```

### Symbolic Differentiation

```cpp
auto x = arg<0>;
auto f = exp(sin(x));                           // e^(sin(x))
auto df = derivative(f).wrt<0>();               // cos(x)·e^(sin(x))

// Evaluate at x = 0
std::array<double, 1> args{0.0};
double slope = df.eval(args);                   // = 1.0
```

### Named Variables

```cpp
auto [x, y] = vars_xy();                        // Named variables for debugging
auto f = sin(x) * cos(y);
std::cout << f.to_string();                     // "(* (sin x) (cos y))"
```

### Integration Methods

```cpp
using namespace limes::methods;

auto I = integral(sin(x)).over<0>(0.0, pi);

// Choose your method
auto r1 = I.eval();                             // Default adaptive
auto r2 = I.eval(gauss<15>());                  // 15-point Gauss-Legendre
auto r3 = I.eval(adaptive(1e-12));              // Adaptive with tolerance
```

### Multi-dimensional Integration

```cpp
auto x = arg<0>;
auto y = arg<1>;

// Double integral: ∫₀¹∫₀ˣ xy dy dx
auto I = integral(x * y)
    .over<1>(0.0, x)                            // y from 0 to x
    .over<0>(0.0, 1.0);                         // x from 0 to 1

auto result = I.eval();                         // = 1/8
```

### Box Integration (Monte Carlo)

```cpp
auto x = arg<0>;
auto y = arg<1>;

// Volume under paraboloid over unit square
auto I = integral(x*x + y*y)
    .over_box({0.0, 1.0}, {0.0, 1.0});

auto result = I.eval(monte_carlo(100000));      // ≈ 2/3

// With constraint (unit circle)
auto disk = integral(One<double>{})
    .over_box({-1.0, 1.0}, {-1.0, 1.0})
    .where([](double x, double y) { return x*x + y*y <= 1.0; });

auto area = disk.eval(monte_carlo(100000));     // ≈ π
```

### Separable Integral Composition

```cpp
auto x = arg<0>;
auto y = arg<1>;

// Independent integrals
auto I = integral(sin(x)).over<0>(0.0, pi);     // = 2
auto J = integral(exp(y)).over<1>(0.0, 1.0);   // = e - 1

// Multiply: evaluates each independently, then multiplies
auto IJ = I * J;                                // ProductIntegral
auto result = IJ.eval();                        // = 2(e - 1) ≈ 3.44
```

## Installation

limes is header-only. Copy `include/limes` to your project or use CMake:

```cmake
find_package(limes REQUIRED)
target_link_libraries(your_target PRIVATE limes::limes)
```

## Architecture

```
limes/
├── limes.hpp                    # Main entry point
├── expr/                        # Expression layer (user-facing API)
│   ├── expr.hpp                 # Aggregated expression header
│   ├── nodes/                   # Expression nodes (Const, Var, Binary, Unary, ...)
│   ├── derivative.hpp           # Symbolic differentiation
│   ├── derivative_builder.hpp   # derivative(f).wrt<D>() API
│   ├── integral.hpp             # Integration expressions
│   ├── box_integral.hpp         # N-dimensional box integration
│   ├── product_integral.hpp     # Separable integral composition
│   └── analysis.hpp             # Variable set analysis
├── methods/                     # Integration method objects
│   └── methods.hpp              # gauss<N>(), monte_carlo(), adaptive()
└── algorithms/                  # Numerical backend
    ├── concepts/                # Field, Accumulator, QuadratureRule
    ├── accumulators/            # Kahan, Neumaier, Klein
    ├── quadrature/              # Gauss-Legendre, Kronrod, Clenshaw-Curtis
    └── integrators/             # Adaptive, Romberg, Tanh-Sinh
```

## Namespace Structure

- `limes` (alias: `li`) — Root namespace
- `limes::expr` — Expression layer (composable calculus)
- `limes::methods` — Integration method objects
- `limes::algorithms` — Low-level numerical backend
- `limes::algorithms::quadrature` — Quadrature rules
- `limes::algorithms::accumulators` — Compensated summation

## Design Philosophy

limes is inspired by Stepanov's approach to generic programming:

1. **Expressions as algebraic objects** — Integrals compose according to mathematical laws (linearity, Fubini, separability)
2. **Separation of concerns** — The *what* (expression) is separate from the *how* (numerical method)
3. **Concepts over inheritance** — Clear, minimal requirements enable generic programming
4. **Compile-time where possible** — Symbolic differentiation and expression simplification at compile time

## Building

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build
ctest --test-dir build --output-on-failure
```

## Requirements

- C++20 compiler (GCC 10+, Clang 12+, MSVC 19.29+)
- CMake 3.16+
- Google Test (for tests, fetched automatically)

## Documentation

### Online Documentation

Full documentation is available at: **[GitHub Pages](https://queelius.github.io/limes/)**

### Local Build

Build the API documentation locally with Doxygen:

```bash
cmake -S . -B build -DBUILD_DOCS=ON
cmake --build build --target docs
# Open build/docs/html/index.html in your browser
```

### Guides

| Guide | Description |
|-------|-------------|
| [Motivation](docs/motivation.md) | Why limes? The problem it solves |
| [Tutorial](docs/tutorial.md) | Step-by-step introduction |
| [Examples](docs/examples.md) | Complete worked examples |
| [Algorithms](docs/algorithms.md) | Numerical method details |
| [Design](docs/LIMES_DESIGN.md) | Full API design and philosophy |
| [About](docs/about.md) | Author and project info |

## License

MIT
