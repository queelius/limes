# Tutorial {#tutorial}

This tutorial walks you through limes step by step.

## Installation

limes is header-only. Clone the repository and point your compiler at the `include` directory:

```bash
git clone https://github.com/queelius/limes.git
cd limes
```

### Using CMake

```cmake
# In your CMakeLists.txt
add_subdirectory(limes)
target_link_libraries(your_target PRIVATE limes::limes)
```

Or if installed system-wide:

```cmake
find_package(limes REQUIRED)
target_link_libraries(your_target PRIVATE limes::limes)
```

### Header Only

Just add the include path:

```bash
g++ -std=c++20 -I/path/to/limes/include your_program.cpp
```

## Step 1: Your First Expression

Include the library and use the expression namespace:

```cpp
#include <limes/limes.hpp>
#include <iostream>

int main() {
    using namespace limes::expr;

    // Create a variable
    auto x = arg<0>;

    // Build an expression: f(x) = x²
    auto f = x * x;

    // Evaluate at x = 3.0
    std::array<double, 1> args{3.0};
    double result = f.eval(args);

    std::cout << "f(3) = " << result << "\n";  // Output: f(3) = 9
}
```

### What's Happening?

- `arg<0>` creates a variable representing the first argument
- `x * x` builds a `Binary<Mul, Var, Var>` expression at compile time
- `f.eval(args)` evaluates the expression tree with the given arguments

The expression is **not** a function pointer—it's a type that carries structure.

## Step 2: More Complex Expressions

Combine primitives to build complex expressions:

```cpp
using namespace limes::expr;

auto x = arg<0>;

// Trigonometric
auto f = sin(x);           // sin(x)
auto g = cos(x * x);       // cos(x²)

// Exponential and logarithmic
auto h = exp(-x * x);      // e^(-x²)
auto j = log(1 + x);       // ln(1 + x)

// Combinations
auto k = sin(x) * exp(-x);              // sin(x)·e^(-x)
auto m = sqrt(1 - x*x);                 // √(1 - x²)
auto n = exp(sin(x)) + log(cos(x));     // e^(sin(x)) + ln(cos(x))
```

Available primitives: `sin`, `cos`, `tan`, `exp`, `log`, `sqrt`, `abs`, `pow`.

## Step 3: Symbolic Differentiation

Differentiate expressions symbolically:

```cpp
using namespace limes::expr;

auto x = arg<0>;
auto f = sin(x * x);  // sin(x²)

// Differentiate with respect to x (dimension 0)
auto df = derivative(f).wrt<0>();  // 2x·cos(x²)

// Evaluate the derivative at x = 1
std::array<double, 1> args{1.0};
double slope = df.eval(args);  // ≈ 1.0806

std::cout << df.to_string() << "\n";
// Output: (* (* 2 x) (cos (* x x)))
```

### Chain Rule at Compile Time

The derivative is computed symbolically using the chain rule. For `sin(x²)`:

1. Outer function: `sin(u)` → derivative is `cos(u)`
2. Inner function: `u = x²` → derivative is `2x`
3. Chain rule: `d/dx[sin(x²)] = cos(x²) · 2x`

All of this happens at compile time through template metaprogramming.

## Step 4: Basic Integration

Integrate expressions numerically:

```cpp
using namespace limes::expr;

auto x = arg<0>;
auto f = x * x;  // f(x) = x²

// Integral: ∫₀¹ x² dx
auto I = integral(f).over<0>(0.0, 1.0);

// Evaluate
auto result = I.eval();

std::cout << "Value: " << result.value() << "\n";  // ≈ 0.333333
std::cout << "Error: " << result.error() << "\n";  // Small error estimate
```

### The Fluent Builder

The pattern `integral(f).over<Dim>(lo, hi)` is a fluent builder:

1. `integral(f)` creates an `IntegralBuilder` wrapping the integrand
2. `.over<0>(0.0, 1.0)` specifies: integrate dimension 0 from 0 to 1
3. The result is an `Integral` expression node

## Step 5: Integration Methods

Choose different numerical methods:

```cpp
using namespace limes::expr;
using namespace limes::methods;

auto x = arg<0>;
auto I = integral(sin(x)).over<0>(0.0, 3.14159);

// Default adaptive integration
auto r1 = I.eval();

// 7-point Gauss-Legendre quadrature
auto r2 = I.eval(gauss<7>());

// 15-point Gauss-Legendre
auto r3 = I.eval(gauss<15>());

// Monte Carlo with 100,000 samples
auto r4 = I.eval(monte_carlo_method(100000));

// Adaptive with custom tolerance
auto r5 = I.eval(adaptive_method(1e-12));

// Simpson's rule with 100 subdivisions
auto r6 = I.eval(simpson_method<100>());
```

### Method Objects

Methods are **objects**, not just type tags. You can configure them:

```cpp
// Monte Carlo with seed for reproducibility
auto mc = monte_carlo_method(100000).with_seed(42);
auto result = I.eval(mc);

// Adaptive with custom tolerance and max subdivisions
auto adaptive = adaptive_method()
    .with_tolerance(1e-14)
    .with_max_subdivisions(5000);
```

## Step 6: Multivariate Expressions

Work with multiple variables:

```cpp
using namespace limes::expr;

auto x = arg<0>;
auto y = arg<1>;

// f(x, y) = x² + y²
auto f = x*x + y*y;

// Evaluate at (x, y) = (3, 4)
std::array<double, 2> args{3.0, 4.0};
double result = f.eval(args);  // 25.0

// Partial derivatives
auto df_dx = derivative(f).wrt<0>();  // 2x
auto df_dy = derivative(f).wrt<1>();  // 2y

// Gradient (all first partials)
auto [fx, fy] = derivative(f).gradient();
```

## Step 7: Double Integrals

Chain `.over()` calls for iterated integrals:

```cpp
using namespace limes::expr;

auto x = arg<0>;
auto y = arg<1>;

// ∫₀¹ ∫₀¹ xy dx dy = 1/4
auto I = integral(x * y)
    .over<0>(0.0, 1.0)   // Integrate x from 0 to 1
    .over<1>(0.0, 1.0);  // Integrate y from 0 to 1

auto result = I.eval();
std::cout << result.value() << "\n";  // ≈ 0.25
```

### Dependent Bounds

The inner bound can depend on outer variables:

```cpp
// ∫₀¹ ∫₀ˣ xy dy dx (triangle region)
auto I = integral(x * y)
    .over<1>(0.0, x)     // y from 0 to x (depends on x!)
    .over<0>(0.0, 1.0);  // x from 0 to 1

auto result = I.eval();  // = 1/8
```

## Step 8: Box Integration

For N-dimensional rectangular regions, use Monte Carlo:

```cpp
using namespace limes::expr;
using namespace limes::methods;

auto x = arg<0>;
auto y = arg<1>;
auto z = arg<2>;

// Volume integral over unit cube
auto I = integral(x * y * z)
    .over_box({{0.0, 1.0}, {0.0, 1.0}, {0.0, 1.0}});

auto result = I.eval(monte_carlo_method(100000));
std::cout << result.value() << "\n";  // ≈ 0.125
```

### Constrained Regions

Use `.where()` for non-rectangular regions:

```cpp
// Area of unit circle using Monte Carlo
auto I = integral(One<double>{})
    .over_box({{-1.0, 1.0}, {-1.0, 1.0}})
    .where([](double x, double y) {
        return x*x + y*y <= 1.0;  // Inside unit circle
    });

auto result = I.eval(100000);  // ≈ π
```

## Step 9: Separable Integrals

When integrals are over independent variables, multiply them:

```cpp
using namespace limes::expr;

auto x = arg<0>;
auto y = arg<1>;

// These are independent!
auto I = integral(sin(x)).over<0>(0.0, 3.14159);  // depends only on x
auto J = integral(exp(y)).over<1>(0.0, 1.0);      // depends only on y

// ProductIntegral: evaluates each separately
auto IJ = I * J;
auto result = IJ.eval();

// Equivalent to: (∫sin(x)dx) × (∫exp(y)dy)
// Much faster than nested integration!
```

The library verifies independence at compile time—if the integrals share variables, you get a clear error message.

## Step 10: Named Variables

For better debugging, use named variables:

```cpp
using namespace limes::expr;

auto x = var(0, "x");
auto y = var(1, "y");

auto f = sin(x) * cos(y);
std::cout << f.to_string() << "\n";
// Output: (* (sin x) (cos y))

// Compare with positional:
auto a = arg<0>;
auto b = arg<1>;
auto g = sin(a) * cos(b);
std::cout << g.to_string() << "\n";
// Output: (* (sin (arg 0)) (cos (arg 1)))
```

## Next Steps

Now that you understand the basics:

- See @ref examples for complete worked examples
- Explore @ref limes::methods for all integration methods
- Read @ref motivation for the design philosophy
- Check the API reference for detailed documentation

## Quick Reference

```cpp
// Variables
auto x = arg<0>;              // Positional
auto x = var(0, "x");         // Named

// Expressions
sin(x), cos(x), exp(x), log(x), sqrt(x), abs(x)
x + y, x - y, x * y, x / y
pow(x, 2), pow(x, y)

// Differentiation
derivative(f).wrt<0>()        // ∂f/∂x
derivative(f).wrt<0, 1>()     // ∂²f/∂x∂y
derivative(f).gradient()       // (∂f/∂x₀, ∂f/∂x₁, ...)

// Integration
integral(f).over<0>(a, b)     // ∫ₐᵇ f dx
integral(f).over_box(bounds)  // N-dimensional box
I * J                         // Product of independent integrals

// Evaluation
I.eval()                      // Default method
I.eval(gauss<7>())            // Specific method
result.value()                // Integral value
result.error()                // Error estimate
```
