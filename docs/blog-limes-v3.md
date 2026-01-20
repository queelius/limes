# Calculus as Algebra: Designing a Composable Integration Library in C++20

*Building mathematical expressions that know their own derivatives and can be integrated with any method.*

---

## The Problem with Numerical Integration APIs

Most numerical integration libraries have an API that looks something like this:

```cpp
double result = integrate(f, 0.0, 1.0);  // f is a function pointer or lambda
```

This works, but it treats the integrand as an opaque black box. The library can't inspect the function's structure, can't simplify `∫(2x + 3x) dx` to `∫5x dx`, can't recognize that `∫∫f(x)g(y) dxdy` factors into `(∫f dx)(∫g dy)`, and can't symbolically differentiate through the integral.

What if we could build integrals that understand their own algebraic structure?

## limes v3.0: Integrals as First-Class Objects

With the v3.0 release of [limes](https://github.com/...), we now have a C++20 library where calculus operations are composable algebraic objects:

```cpp
using namespace limes::expr;

auto x = arg<0>;
auto f = sin(x * x);                       // Build an expression

auto df = derivative(f).wrt<0>();          // Symbolic differentiation
auto I = integral(f).over<0>(0.0, 1.0);    // Integration expression

auto result = I.eval();                    // Numerical evaluation
```

The key insight: **expressions carry their structure**. The type `sin(x * x)` isn't just "some function"—it's literally `UnaryFunc<SinTag, Binary<Mul, Var<0>, Var<0>>>`. The derivative is computed at compile time by recursively applying the chain rule through that type structure.

## Six Features That Make This Work

### 1. Uniform Builder Pattern

Both differentiation and integration use fluent builders that mirror mathematical notation:

```cpp
// Derivative: d²f/dx∂y
auto d2f = derivative(f).wrt<0>().wrt<1>();
auto d2f = derivative(f).wrt<0, 1>();      // Same thing

// Integral: ∫₀¹∫₀ˣ f(x,y) dy dx
auto I = integral(f)
    .over<1>(0.0, x)                       // Inner: integrate y from 0 to x
    .over<0>(0.0, 1.0);                    // Outer: integrate x from 0 to 1
```

The dimension parameter (`<0>`, `<1>`) specifies which variable to differentiate/integrate with respect to. This handles multivariate expressions naturally.

### 2. Methods as First-Class Objects

Instead of hardcoding numerical methods, limes treats them as composable objects:

```cpp
using namespace limes::methods;

// Basic methods
I.eval(gauss<7>());                        // 7-point Gauss-Legendre
I.eval(simpson());                         // Simpson's rule
I.eval(monte_carlo(10000));                // Monte Carlo with 10k samples

// Composed methods
I.eval(adaptive(gauss<5>(), 1e-12));       // Adaptive subdivision

// Reproducible Monte Carlo
I.eval(monte_carlo(10000).with_seed(42));
```

Methods satisfy a concept, so users can define their own:

```cpp
template<typename M, typename T>
concept IntegrationMethod = requires(M m, auto f, T a, T b) {
    { m(f, a, b) } -> std::same_as<integration_result<T>>;
};
```

### 3. Named Variables for Debugging

Mathematical expressions can be hard to debug when everything is `x0`, `x1`, `x2`. Named variables solve this:

```cpp
auto x = var(0, "x");
auto y = var(1, "y");
auto f = sin(x) * cos(y);

std::cout << f.to_string();  // "(* (sin x) (cos y))"
```

Or use the convenience functions:

```cpp
auto [x, y, z] = vars_xyz();
```

### 4. N-Dimensional Box Integration

For rectangular regions, Monte Carlo integration is often the most practical approach:

```cpp
// Integrate over the 3D unit cube
auto I = integral(f).over_box(
    {0.0, 1.0},   // x bounds
    {0.0, 1.0},   // y bounds
    {0.0, 1.0}    // z bounds
);

auto result = I.eval(monte_carlo(100000));
```

For non-rectangular regions, add a constraint:

```cpp
// Integrate over the unit disk
auto disk = integral(One<double>{})
    .over_box({-1.0, 1.0}, {-1.0, 1.0})
    .where([](double x, double y) {
        return x*x + y*y <= 1.0;
    });

auto area = disk.eval(monte_carlo(100000));  // ≈ π
```

The `.where()` clause triggers rejection sampling—points outside the region are discarded.

### 5. Separable Integral Composition

Here's where the algebraic structure really pays off. For independent integrals:

```cpp
auto I = integral(sin(x)).over<0>(0.0, pi);    // Depends only on x
auto J = integral(exp(y)).over<1>(0.0, 1.0);   // Depends only on y

auto IJ = I * J;  // ProductIntegral
```

The library verifies at compile time that `I` and `J` are over independent variables using a `variable_set` bitset:

```cpp
static_assert(
    (variable_set_v<I::integrand_type> &
     variable_set_v<J::integrand_type>) == 0,
    "Integrals must be over independent variables"
);
```

When you evaluate `IJ`, it computes `I` and `J` separately (potentially in parallel in a future version) and multiplies the results. Error propagation follows the standard formula: `δ(ab) = |b|δa + |a|δb`.

### 6. Compile-Time Expression Simplification

The type system does constant folding and algebraic simplification:

```cpp
auto x = arg<0>;
auto f = x * One<double>{};    // Type: Var<0, double>, not Binary<Mul, ...>
auto g = x + Zero<double>{};   // Type: Var<0, double>
auto h = Zero<double>{} * f;   // Type: Zero<double>
```

This happens through `if constexpr` chains in the operator overloads:

```cpp
template<typename L, typename R>
constexpr auto operator*(L l, R r) {
    if constexpr (is_zero_v<L>) {
        return Zero<typename R::value_type>{};
    } else if constexpr (is_one_v<R>) {
        return l;
    } else {
        return Binary<Mul, L, R>{l, r};
    }
}
```

## Implementation Highlights

### Variable Set Analysis

Every expression type carries a compile-time bitset of which dimensions it depends on:

```cpp
template<typename E>
struct variable_set;

template<std::size_t N, typename T>
struct variable_set<Var<N, T>> {
    static constexpr uint64_t value = (1ULL << N);
};

template<typename Op, typename L, typename R>
struct variable_set<Binary<Op, L, R>> {
    static constexpr uint64_t value =
        variable_set<L>::value | variable_set<R>::value;
};

template<typename E, std::size_t Dim, typename Lo, typename Hi>
struct variable_set<Integral<E, Dim, Lo, Hi>> {
    // Remove integrated dimension, add bound dependencies
    static constexpr uint64_t value =
        (variable_set<E>::value & ~(1ULL << Dim)) |
        variable_set<Lo>::value | variable_set<Hi>::value;
};
```

This enables:
- Checking if expressions are separable
- Verifying integral independence
- Detecting if bounds depend on integration variables

### Symbolic Differentiation via Chain Rule

Each expression node knows how to differentiate itself:

```cpp
template<typename E>
struct UnaryFunc<SinTag, E> {
    E child;

    template<std::size_t Dim>
    constexpr auto derivative() const {
        // d/dx[sin(f)] = cos(f) * df/dx
        auto df = child.template derivative<Dim>();
        return cos(child) * df;
    }
};
```

The result is a new expression type. Constant expressions (`Zero<T>`, `One<T>`) short-circuit the chain rule—no runtime overhead for multiplying by 1 or adding 0.

### Bounds as Expressions

Integration bounds can be constants or expressions:

```cpp
// Constant bounds
integral(f).over<0>(0.0, 1.0);

// Expression bounds (depends on outer variable)
integral(f).over<1>(Zero<double>{}, x);   // y from 0 to x

// Lambda bounds
integral(f).over<1>(0.0, [](auto x) { return x*x; });
```

This enables iterated integrals over non-rectangular regions like triangles, cones, or arbitrary shapes defined by inequalities.

## What's Next?

The v3.0 API is stable, but several enhancements are planned:

1. **Parallel evaluation** — `ProductIntegral::eval()` could evaluate both children concurrently
2. **Automatic separability detection** — Recognize `∫∫f(x)g(y) dxdy` and factor automatically
3. **Symbolic antiderivatives** — For known primitives, avoid numerical integration entirely
4. **GPU acceleration** — Monte Carlo integration is embarrassingly parallel

## Try It

limes is header-only and requires only C++20:

```bash
git clone https://github.com/.../limes
cd limes
cmake -B build && cmake --build build
ctest --test-dir build
```

The test suite has 489 tests covering expression evaluation, differentiation, integration, box integrals, and separable composition.

---

## Appendix: The Design Philosophy

limes follows Stepanov's approach to generic programming:

1. **Find the algebraic structure** — Integrals satisfy linearity, additivity, Fubini's theorem. Make these laws explicit in the API.

2. **Define minimal concepts** — An `IntegrationMethod` just needs to compute a result from a function and bounds. An `ExprNode` just needs `eval()`, `derivative<Dim>()`, and `to_string()`.

3. **Let composition emerge** — Don't design a "separable integral" feature. Design expression types with variable tracking, then let `operator*` check independence as a natural consequence.

4. **Separate what from how** — The expression `integral(sin(x)).over<0>(0, pi)` describes *what* to compute. The method `gauss<15>()` describes *how* to compute it. They combine at evaluation time.

The result is a library where calculus operations are as composable as arithmetic—because mathematically, that's exactly what they are.
