# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Vision

**limes** is a pedagogical C++20 header-only library demonstrating generic programming principles for numerical calculus. Inspired by Stepanov's approach: define minimal algebraic concepts, implement algorithms generically against those concepts, and compose solutions from simple, orthogonal parts.

The goal is **clarity over feature count**. Each component should be a teachable example of good API design.

## Build Commands

```bash
# Configure and build
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build

# Run all tests
ctest --test-dir build --output-on-failure

# Run specific test suite
./build/tests/test_expr
./build/tests/test_integrators

# Run tests matching a pattern
ctest --test-dir build -R "BoxIntegral"
ctest --test-dir build -R "ProductIntegral"

# Coverage (Debug mode)
cmake -S . -B build -DCMAKE_BUILD_TYPE=Debug
cmake --build build --target coverage

# Build and run demo
cmake --build build --target demo && ./build/examples/demo
```

## Architecture: Two-Layer Design

```
┌─────────────────────────────────────────────────────────────────┐
│  limes::expr — Expression Layer (user-facing API)              │
│  Composable calculus expressions with symbolic differentiation  │
│  and numerical integration                                      │
├─────────────────────────────────────────────────────────────────┤
│  limes::methods — Integration Method Objects                    │
│  gauss<N>(), monte_carlo(), adaptive()                         │
├─────────────────────────────────────────────────────────────────┤
│  limes::algorithms — Numerical Backend                          │
│  Integrators, quadrature rules, accumulators                    │
└─────────────────────────────────────────────────────────────────┘
```

### Expression Layer (`limes::expr`)

Expression templates for zero-overhead calculus expressions:

```cpp
using namespace limes::expr;

auto x = arg<0>;                              // Variable (positional)
auto f = sin(x*x) + exp(-x);                  // Composable expression

// v3.0 API: fluent builders
auto df = derivative(f).wrt<0>();             // Symbolic differentiation
auto I = integral(f).over<0>(0.0, 1.0);       // Integration expression
auto result = I.eval();                       // Returns integration_result<T>

// Methods as objects
auto r = I.eval(methods::gauss<15>());        // Specific method
auto r = I.eval(methods::monte_carlo(10000)); // Monte Carlo
```

**Key types:**
- `Const<T>`, `Var<N, T>`, `NamedVar<T>` — Leaf nodes
- `Binary<Op, L, R>`, `Unary<Op, E>` — Composite nodes
- `UnaryFunc<Tag, E>` — Primitives: `exp`, `log`, `sin`, `cos`, `sqrt`, `abs`
- `Integral<E, Dim, Lo, Hi>` — 1D integration node
- `BoxIntegral<E, Dims>` — N-D box integration (Monte Carlo)
- `ProductIntegral<I1, I2>` — Separable integral composition

### Methods Layer (`limes::methods`)

Integration methods as first-class objects:

```cpp
using namespace limes::methods;

gauss_legendre<15, double>{}      // Or use gauss<15>()
monte_carlo<double>{10000}        // Or monte_carlo(n).with_seed(42)
adaptive<double>{1e-10}           // Adaptive with tolerance
```

### Algorithms Layer (`limes::algorithms`)

Low-level numerical methods with pluggable components:

```cpp
using namespace limes::algorithms;

adaptive_integrator<double> integrator;
auto result = integrator(f, 0.0, 1.0, 1e-10);
```

**Key abstractions:**
- `Accumulator` concept — Precision strategies (Kahan, Neumaier, Klein)
- `QuadratureRule` concept — Node/weight generators
- `Integrator` concept — Combines rules and accumulators

## Namespace Structure

- `limes` (alias: `li`) — Root namespace
- `limes::expr` — Expression layer (user-facing)
- `limes::methods` — Integration method objects
- `limes::algorithms` — Numerical backend
- `limes::algorithms::concepts` — C++20 concepts
- `limes::algorithms::accumulators` — Precision control
- `limes::algorithms::quadrature` — Quadrature rules

## Key Design Patterns

### Variable Set Analysis (`analysis.hpp`)

Compile-time bitset tracking which dimensions an expression depends on:

```cpp
variable_set_v<Var<0, double>>           // = 0b001
variable_set_v<Binary<Mul, Var<0>, Var<1>>>  // = 0b011
```

Used for: separability detection, independence verification, bound dependency analysis.

### Expression Simplification

Operators use `if constexpr` for compile-time simplification:
- `Zero * anything = Zero`
- `One * x = x`
- `Const + Const = Const` (constant folding)

### is_integral Type Trait

Forward-declared in `nodes/binary.hpp`, specialized in `integral.hpp`. This allows `operator*` in `binary.hpp` to exclude `Integral` types (which use `ProductIntegral` instead).

## Extension Points

**New expression primitive:** Add a tag and specialize `UnaryFunc<Tag, E>` with `eval()`, `derivative<Dim>()`, and `to_string()`.

**New accumulator:** Implement `operator+=(T)`, `operator()() const`, and default constructor.

**New quadrature rule:** Implement `size()`, `weight(i)`, `abscissa(i)`.

**New integration method:** Model `IntegrationMethod` concept in `methods/concepts.hpp`.

## Testing

Google Test with type-parameterized tests for float/double/long double.

```bash
# Run expression layer tests
ctest --test-dir build -R "Expr"

# Run specific test groups
ctest --test-dir build -R "BoxIntegral"
ctest --test-dir build -R "ProductIntegral"
ctest --test-dir build -R "DerivativeBuilder"

# Verbose output
ctest --test-dir build -V
```

Use type-aware tolerances:
```cpp
T tol = std::is_same_v<T, float> ? T(1e-5) : T(1e-10);
```

## API Conventions

- **eval() vs evaluate()**: Use `eval()`. The `evaluate()` method is deprecated.
- **Builder pattern**: `derivative(f).wrt<0>()`, `integral(f).over<0>(a, b)`
- **Method objects**: Pass to `eval()`: `I.eval(methods::gauss<7>())`
- **Named variables**: `var(0, "x")` or `named<0>("x")` for debugging
