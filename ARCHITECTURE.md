# limes - Architecture Overview

## Design Philosophy

limes is a pedagogical C++20 header-only library for composable calculus expressions. Inspired by Stepanov's approach to generic programming:

- **Minimal concepts**: Define clear algebraic contracts via C++20 concepts
- **Generic algorithms**: Implement against concepts, not concrete types
- **Composability**: Expressions, methods, and accumulators combine freely
- **Clarity over features**: Each component is a teachable example of good API design

## Three-Layer Architecture

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

### Layer 1: Expression Layer (`include/limes/expr/`)

Expression templates providing zero-overhead composable calculus:

**Node types** (`expr/nodes/`):
- `Const<T>`, `Zero<T>`, `One<T>` — Constant leaf nodes
- `Var<N, T>` — Positional variable reference (compile-time index)
- `NamedVar<T>`, `StaticNamedVar<N, T>` — Variables with debug names
- `Binary<Op, L, R>` — Binary operators (+, -, *, /) with compile-time simplification
- `Unary<Op, E>` — Unary negation
- `UnaryFunc<Tag, E>` — Math primitives (sin, cos, exp, log, sqrt, abs)
- `BinaryFunc<Tag, L, R>` — Two-argument functions (pow, atan2, min, max)
- `Pow<N, E>` — Integer power with optimized derivative
- `Bound<E, Dim>` — Variable binding (partial application)
- `Conditional<C, T, F>` — Conditional expressions, `sign()`, `heaviside()`
- `FiniteSum<E, Idx>`, `FiniteProduct<E, Idx>` — Summation/product over index ranges

**Integration nodes** (`expr/`):
- `Integral<E, Dim, Lo, Hi>` — 1D definite integral with dependent bounds
- `BoxIntegral<E, Dims>` — N-D rectangular integration (Monte Carlo)
- `ConstrainedBoxIntegral<E, Dims, Constraint>` — Constrained region integration
- `ProductIntegral<I1, I2>` — Separable integral composition

**Analysis and differentiation** (`expr/`):
- `analysis.hpp` — Compile-time variable set tracking (bitset), separability detection
- `derivative.hpp` — Chain rule differentiation for all node types
- `derivative_builder.hpp` — Fluent API: `derivative(f).wrt<0>()`
- `antiderivatives.hpp` — Known closed-form antiderivatives

### Layer 2: Methods Layer (`include/limes/methods/`)

Integration methods as first-class objects modelling the `IntegrationMethod` concept:

- `gauss_legendre<N, T>` — N-point Gauss-Legendre quadrature
- `monte_carlo<T>` — Monte Carlo with configurable samples and seed
- `adaptive<T>` — Adaptive subdivision with tolerance control

Factory functions: `gauss<N>()`, `monte_carlo(n)`, `adaptive(tol)`

### Layer 3: Algorithms Layer (`include/limes/algorithms/`)

Low-level numerical backend with pluggable components:

**Concepts** (`algorithms/concepts/`):
- `Field` — Types supporting +, -, *, / operations
- `Accumulator` — Types accumulating values with controlled precision
- `QuadratureRule` — Types providing integration nodes and weights
- `UnivariateFunction` — Single-variable callable

**Accumulators** (`algorithms/accumulators/`):
- `simple_accumulator<T>` — Direct summation
- `kahan_accumulator<T>` — Kahan compensated summation
- `neumaier_accumulator<T>` — Improved Kahan (handles large+small)
- `klein_accumulator<T>` — Second-order compensation
- `pairwise_accumulator<T>` — Cache-friendly pairwise summation

**Quadrature rules** (`algorithms/quadrature/`):
- `gauss_legendre<T, N>` — Gauss-Legendre for N = 2, 3, 5, 7, 15
- `gauss_kronrod<T, N>` — Embedded error estimation (N = 7, 15)
- `clenshaw_curtis<T, N>` — Nested FFT-friendly rules
- `tanh_sinh_nodes<T, N>` — Double exponential for endpoint singularities

**Integrators** (`algorithms/integrators/`):
- `quadrature_integrator<T, Rule, Acc>` — Basic quadrature integration
- `adaptive_integrator<T, Acc>` — Adaptive interval subdivision
- `romberg_integrator<T, Acc>` — Richardson extrapolation
- `tanh_sinh_integrator<T, Acc>` — Double exponential for semi/infinite intervals

## Key Design Patterns

### Expression Templates with Compile-Time Simplification

Operators use `if constexpr` to simplify at compile time:
- `Zero * anything = Zero`
- `One * x = x`
- `Const + Const = Const` (constant folding)

### Variable Set Analysis

Compile-time bitset tracking which dimensions an expression depends on:
```cpp
variable_set_v<Var<0, double>>                     // = 0b001
variable_set_v<Binary<Mul, Var<0>, Var<1>>>        // = 0b011
```
Used for separability detection, independence verification, and bound dependency analysis.

### is_integral Type Trait

Forward-declared in `nodes/binary.hpp`, specialized in `integral.hpp`. This allows `operator*` in `binary.hpp` to route `Integral * Integral` to `ProductIntegral` instead of `Binary<Mul>`.

### Detail Helpers

Shared implementation helpers reduce duplication:
- `detail::tag_name<Tag>()` — Maps function tags to string names
- `detail::make_integrand_fn<Dim>()` — Builds integrand lambdas for eval
- `detail::mc_box_sampler<T, Dims>` — Monte Carlo RNG and sampling
- `detail::make_extended_args()` — Argument tuple construction for sum/product

## Namespace Structure

- `limes` (alias: `li`) — Root namespace
- `limes::expr` — Expression layer (user-facing)
- `limes::methods` — Integration method objects
- `limes::algorithms` — Numerical backend
- `limes::algorithms::concepts` — C++20 concepts
- `limes::algorithms::accumulators` — Precision control
- `limes::algorithms::quadrature` — Quadrature rules

## Extension Points

**New expression primitive:** Add a tag struct and specialize `UnaryFunc<Tag, E>` with `eval()`, `derivative<Dim>()`, and `to_string()`.

**New accumulator:** Implement `operator+=(T)`, `operator()() const`, and default constructor.

**New quadrature rule:** Implement `size()`, `weight(i)`, `abscissa(i)`.

**New integration method:** Model `IntegrationMethod` concept (provide `operator()` taking a callable and bounds).

## File Organization

```
include/limes/
├── limes.hpp                         # Main entry point
├── fwd.hpp                           # Forward declarations
├── expr/
│   ├── expr.hpp                      # Aggregated expression header
│   ├── nodes/                        # 11 expression node types
│   ├── analysis.hpp                  # Variable set analysis
│   ├── antiderivatives.hpp           # Known antiderivatives
│   ├── concepts.hpp                  # Expression concepts
│   ├── derivative.hpp                # Differentiation rules
│   ├── derivative_builder.hpp        # Builder API
│   ├── integral.hpp                  # 1D integration
│   ├── box_integral.hpp              # N-D box integration
│   ├── product_integral.hpp          # Separable integrals
│   └── to_string.hpp                 # Expression formatting
├── methods/
│   ├── concepts.hpp                  # IntegrationMethod concept
│   └── methods.hpp                   # Method implementations
└── algorithms/
    ├── algorithms.hpp                # Aggregated header
    ├── concepts/concepts.hpp         # Core concepts
    ├── core/result.hpp               # integration_result<T>
    ├── accumulators/accumulators.hpp  # 5 accumulator types
    ├── quadrature/quadrature.hpp     # 4 quadrature families
    └── integrators/integrators.hpp   # 4 integrator types
```
