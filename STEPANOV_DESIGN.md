# Stepanov Design Principles for Algebraic Integrators

## Alexander Stepanov's Philosophy

Core principles from *Elements of Programming* and STL design:

1. **Generic Programming** - Algorithms independent of data structures
2. **Concepts** - Explicit requirements on template parameters
3. **Regular Types** - Value semantics, copyable, comparable
4. **Algebraic Structures** - Identify and use mathematical structure (monoids, groups, etc.)
5. **Algorithm/Data Separation** - Orthogonal components that compose
6. **Zero Overhead** - Abstractions should be as efficient as hand-coded
7. **Simplicity** - Do one thing well, compose for complexity

## Current Design Alignment

### ✅ What We're Already Doing Right

#### 1. **Generic Programming**
```cpp
// Algorithm (integration) is generic over:
// - Value type (float, double, long double, custom)
// - Accumulator strategy
// - Quadrature rule
template <typename T, typename Rule, typename Accumulator>
class quadrature_integrator;
```

#### 2. **Concepts** (C++20)
```cpp
template <typename T>
concept Field = requires(T a, T b) {
    { a + b } -> std::convertible_to<T>;
    { a - b } -> std::convertible_to<T>;
    { a * b } -> std::convertible_to<T>;
    { a / b } -> std::convertible_to<T>;
};
```

#### 3. **Composability**
```cpp
// Accumulators compose with integrators
using Acc = kahan_accumulator<double>;
using Rule = gauss_kronrod_15<double>;
quadrature_integrator<double, Rule, Acc> integrator;
```

#### 4. **Regular Types**
```cpp
// integration_result is Regular:
// - Default constructible
// - Copyable, movable
// - Equality comparable
struct integration_result {
    T value;
    T error;
    size_t evaluations;
};
```

### 🟡 What Needs Improvement

## Stepanov Principle 1: **Identify Algebraic Structures**

### Accumulators as Monoids

Accumulators should be recognized as **commutative monoids**:
- **Set:** Values of type `T`
- **Operation:** `add(a, b)` or `a + b`
- **Identity:** `zero()`
- **Associative:** `(a + b) + c = a + (b + c)`
- **Commutative:** `a + b = b + a`

```cpp
template <typename T>
concept Monoid = requires(T a, T b) {
    { a + b } -> std::same_as<T>;           // Binary operation
    { T::identity() } -> std::same_as<T>;   // Identity element
    // Associativity and identity laws (axioms, not checkable)
};

template <typename T>
concept CommutativeMonoid = Monoid<T> && requires(T a, T b) {
    // Commutativity (axiom)
    requires std::same_as<decltype(a + b), decltype(b + a)>;
};

// Accumulators satisfy CommutativeMonoid
template <typename T>
class kahan_accumulator {
    T sum = T(0);
    T compensation = T(0);

public:
    static kahan_accumulator identity() { return {}; }  // Zero element

    kahan_accumulator& operator+=(const T& value) {
        T y = value - compensation;
        T t = sum + y;
        compensation = (t - sum) - y;
        sum = t;
        return *this;
    }

    friend kahan_accumulator operator+(kahan_accumulator a, const kahan_accumulator& b) {
        a += b.result();  // Combine accumulators
        return a;
    }

    T result() const { return sum; }
};
```

### Integration as a Fold/Reduce

Integration is a **generalized fold** over function evaluations:

```cpp
// Traditional fold/reduce
template <Monoid M, typename F, typename Range>
M fold(Range range, F transform, M init = M::identity()) {
    for (auto x : range) {
        init = init + transform(x);
    }
    return init;
}

// Integration is fold over quadrature nodes
template <Field T, QuadratureRule Rule, Monoid Acc>
T integrate(Function f, T a, T b, Rule rule, Acc acc = Acc::identity()) {
    auto [nodes, weights] = rule.generate(a, b);

    for (size_t i = 0; i < nodes.size(); ++i) {
        acc += f(nodes[i]) * weights[i];
    }

    return acc.result();
}
```

**Benefit:** Automatic parallelization via monoid associativity!

```cpp
// Parallel reduce because monoids are associative
template <Monoid M>
M parallel_reduce(std::vector<M> partials) {
    // Can combine in any order: ((a+b)+(c+d)) or (a+(b+(c+d)))
    return std::reduce(std::execution::par, partials.begin(), partials.end());
}
```

## Stepanov Principle 2: **Regular Types**

### Make All Types Regular

A type is **Regular** if it has:
- Default constructor
- Copy constructor, copy assignment
- Move constructor, move assignment (C++11)
- Equality comparison `==`
- Total ordering `<` (for use in containers)

```cpp
// Current: integration_result is almost Regular
template <typename T>
struct integration_result {
    T value{};
    T error{};
    size_t evaluations{};

    // Add equality
    friend bool operator==(const integration_result&, const integration_result&) = default;

    // Add ordering (lexicographic on value, then error)
    friend auto operator<=>(const integration_result& a, const integration_result& b) {
        if (auto cmp = a.value <=> b.value; cmp != 0) return cmp;
        return a.error <=> b.error;
    }
};

// Now can use in containers, algorithms
std::vector<integration_result<double>> results;
std::sort(results.begin(), results.end());  // Works!
auto min_error = std::min_element(results.begin(), results.end());
```

### Value Semantics

Prefer value types over reference/pointer semantics:

```cpp
// ❌ Bad: Reference semantics, unclear ownership
class integrator {
    Function* f;  // Who owns this?
    Accumulator& acc;  // Lifetime issues
};

// ✅ Good: Value semantics, clear ownership
template <typename Function, typename Accumulator>
class integrator {
    Function f;      // Owns the function (small, likely lambda)
    Accumulator acc; // Owns the accumulator state
};
```

## Stepanov Principle 3: **Algorithms Over Data Structures**

### Separate Algorithm from Storage

Current design mixes algorithm (how to integrate) with strategy (accumulator):

```cpp
// Current: Integrator is tightly coupled to accumulator type
template <typename T, typename Rule, typename Accumulator>
class quadrature_integrator {
    Rule rule;
    Accumulator acc;  // Specific accumulator type
};
```

Improve with **strategy pattern via concepts**:

```cpp
// Algorithm: General integration using any accumulator
template <Field T, QuadratureRule Rule>
class quadrature_algorithm {
    Rule rule;

public:
    // Generic over any Accumulator satisfying concept
    template <Accumulator Acc, Function F>
    auto operator()(F f, T a, T b, Acc acc = Acc::identity()) const {
        auto [nodes, weights] = rule.generate(a, b);

        for (size_t i = 0; i < nodes.size(); ++i) {
            acc += f(nodes[i]) * weights[i];
        }

        return integration_result{acc.result(), /* ... */};
    }
};
```

Now the algorithm doesn't care about accumulator implementation:

```cpp
auto gauss = quadrature_algorithm<double, gauss_legendre<7>>{};

// Use with different accumulators
auto r1 = gauss(f, 0, 1, simple_accumulator<double>{});
auto r2 = gauss(f, 0, 1, kahan_accumulator<double>{});
auto r3 = gauss(f, 0, 1, neumaier_accumulator<double>{});
```

## Stepanov Principle 4: **Concept-Driven Design**

### Define Precise Concepts

```cpp
namespace algebraic_integrators::concepts {

// Foundational mathematical structures
template <typename T>
concept Semigroup = requires(T a, T b) {
    { a + b } -> std::same_as<T>;
    // Axiom: Associative (a+b)+c = a+(b+c)
};

template <typename T>
concept Monoid = Semigroup<T> && requires {
    { T::identity() } -> std::same_as<T>;
    // Axiom: a + identity() = a
};

template <typename T>
concept Group = Monoid<T> && requires(T a) {
    { -a } -> std::same_as<T>;  // Inverse
    // Axiom: a + (-a) = identity()
};

template <typename T>
concept Field = Group<T> && requires(T a, T b) {
    { a * b } -> std::convertible_to<T>;
    { a / b } -> std::convertible_to<T>;
    // Axioms: Ring and field axioms
};

// Accumulator is a Monoid with reduction
template <typename A>
concept Accumulator = Monoid<A> && requires(A acc, typename A::value_type x) {
    typename A::value_type;
    { acc += x } -> std::same_as<A&>;
    { acc.result() } -> std::same_as<typename A::value_type>;
};

// Quadrature rule generates nodes and weights
template <typename R, typename T>
concept QuadratureRule = requires(R rule, T a, T b) {
    typename R::value_type;
    requires std::same_as<typename R::value_type, T>;
    { rule.nodes() } -> std::ranges::range;
    { rule.weights() } -> std::ranges::range;
    { rule.order() } -> std::convertible_to<int>;
};

// Function that can be integrated
template <typename F, typename T>
concept IntegrableFunction = requires(F f, T x) {
    { f(x) } -> std::convertible_to<T>;
};

} // namespace concepts
```

### Use Concepts for Compile-Time Checks

```cpp
template <concepts::Field T,
          concepts::QuadratureRule<T> Rule,
          concepts::Accumulator Acc>
    requires std::same_as<typename Acc::value_type, T>
auto integrate(concepts::IntegrableFunction<T> auto f, T a, T b) {
    // Implementation
}
```

Benefits:
- Clear requirements
- Better error messages
- Self-documenting
- Enables optimizations

## Stepanov Principle 5: **Transformations and Iterators**

### Functional Transformations

```cpp
// Transform is a functor: domain → codomain
template <typename Domain, typename Codomain>
concept Transform = requires(Transform t, Domain x) {
    { t(x) } -> std::convertible_to<Codomain>;
};

// Bijection: invertible transform
template <typename T, typename U>
concept Bijection = Transform<T, U> && requires(Bijection b, T x, U y) {
    { b.forward(x) } -> std::convertible_to<U>;
    { b.inverse(y) } -> std::convertible_to<T>;
    // Axiom: inverse(forward(x)) = x
};

// Diffeomorphism: smooth bijection (for integration)
template <typename T>
concept Diffeomorphism = Bijection<T, T> && requires(Diffeomorphism d, T x) {
    { d.jacobian(x) } -> std::convertible_to<T>;  // Derivative of forward map
};
```

Use for change of variables:

```cpp
// Generic integration with transformation
template <Field T, Diffeomorphism<T> Transform>
auto integrate_transformed(Function f, T a, T b, Transform transform) {
    auto g = [&](T t) {
        T x = transform.forward(t);
        return f(x) * transform.jacobian(t);
    };

    T ta = transform.inverse(a);
    T tb = transform.inverse(b);

    return integrate(g, ta, tb);
}
```

### Range-Based Design (C++20)

```cpp
// Quadrature as a range of (node, weight) pairs
template <Field T>
struct quadrature_point {
    T node;
    T weight;
};

template <typename R>
concept QuadratureRange = std::ranges::range<R> &&
    requires { typename std::ranges::range_value_t<R>::node; };

// Integration as fold over range
template <QuadratureRange Range, Function F>
auto integrate_range(F f, Range&& quad_points) {
    return std::ranges::fold_left(
        quad_points,
        0.0,
        [&](auto sum, auto qp) { return sum + f(qp.node) * qp.weight; }
    );
}
```

## Stepanov Principle 6: **Efficiency Through Abstraction**

### Zero-Cost Abstractions

```cpp
// Template instantiation means zero runtime overhead
template <typename Acc>
auto integrate_efficient(Function f, double a, double b) {
    Acc acc;  // Type known at compile time

    // Compiler can inline everything
    for (auto [node, weight] : quadrature_points) {
        acc += f(node) * weight;  // No virtual calls!
    }

    return acc.result();
}

// vs. Runtime polymorphism (slower)
auto integrate_virtual(Function f, double a, double b, AccumulatorBase* acc) {
    for (auto [node, weight] : quadrature_points) {
        acc->add(f(node) * weight);  // Virtual call overhead
    }
}
```

### Compile-Time Computation

```cpp
// Quadrature nodes/weights known at compile-time
template <int N>
struct gauss_legendre_rule {
    static constexpr std::array<double, N> nodes = /* computed at compile-time */;
    static constexpr std::array<double, N> weights = /* computed at compile-time */;

    // No runtime computation needed!
};
```

## Stepanov Principle 7: **Generic Algorithms**

### Lifting 1D to ND

Use generic programming to extend algorithms:

```cpp
// 1D integration
template <Field T>
T integrate_1d(Function auto f, T a, T b);

// Lift to 2D via Fubini's theorem (iterated integrals)
template <Field T>
T integrate_2d(auto f, std::array<T, 2> lower, std::array<T, 2> upper) {
    return integrate_1d([&](T x) {
        return integrate_1d([&](T y) {
            return f({x, y});
        }, lower[1], upper[1]);
    }, lower[0], upper[0]);
}

// Generic N-dimensional (recursive)
template <int Dim, Field T>
T integrate_nd(auto f, std::array<T, Dim> lower, std::array<T, Dim> upper) {
    if constexpr (Dim == 1) {
        return integrate_1d([&](T x) { return f({x}); }, lower[0], upper[0]);
    } else {
        return integrate_1d([&](T x) {
            // Integrate over remaining dimensions
            std::array<T, Dim-1> sub_lower, sub_upper;
            std::copy(lower.begin()+1, lower.end(), sub_lower.begin());
            std::copy(upper.begin()+1, upper.end(), sub_upper.begin());

            return integrate_nd<Dim-1>([&](auto sub_x) {
                std::array<T, Dim> full_x;
                full_x[0] = x;
                std::copy(sub_x.begin(), sub_x.end(), full_x.begin()+1);
                return f(full_x);
            }, sub_lower, sub_upper);
        }, lower[0], upper[0]);
    }
}
```

## Implementation Guidelines

### 1. **Concept-First Design**
- Define concepts before implementations
- Document axioms (even if not enforceable)
- Use concepts for all template parameters

### 2. **Regular Types Everywhere**
- Make all types Regular by default
- Document if type cannot be Regular (and why)
- Prefer value semantics

### 3. **Identify Algebraic Structure**
- Is it a semigroup? Monoid? Group?
- Does it form a ring? Field?
- Can operations be parallelized (via associativity)?

### 4. **Separate Algorithm from Strategy**
- Algorithms should be generic over strategies
- Use concepts for strategy requirements
- Avoid coupling to specific implementations

### 5. **Transformations as First-Class**
- Transforms are functions/functors
- Compose naturally
- Use for coordinate changes, domain mapping

### 6. **Test Algebraic Properties**
- Test associativity: `(a+b)+c == a+(b+c)`
- Test identity: `a + identity() == a`
- Test commutativity: `a+b == b+a`
- Use property-based testing

## References

1. **Elements of Programming** - Stepanov & McJones
   - Fundamental algorithms and data structures
   - Generic programming methodology

2. **From Mathematics to Generic Programming** - Stepanov & Rose
   - Historical development of generic programming
   - Algebraic structures in programming

3. **STL Design** - Stepanov's lectures
   - Iterator concepts
   - Algorithm/container separation

4. **C++20 Ranges** - Eric Niebler
   - Modern generic programming with ranges
   - Composable algorithms

## Summary: Stepanov Principles Applied

| Principle | Current Status | Improvement |
|-----------|----------------|-------------|
| Generic Programming | ✅ Strong | Extend to multivariate |
| Concepts | ✅ Good | Add algebraic structure concepts |
| Regular Types | 🟡 Partial | Complete for all types |
| Algebraic Structures | 🟡 Implicit | Make explicit (monoids, etc.) |
| Algorithm/Data Separation | 🟡 Good | Perfect separation |
| Zero Overhead | ✅ Excellent | Maintain |
| Transformations | ✅ Good | Extend to multivariate |

The library is already well-aligned with Stepanov's philosophy! The main improvements are:
1. Make algebraic structures (monoids) explicit
2. Complete Regular type requirements
3. Perfect algorithm/strategy separation
4. Extend generic algorithms to multivariate
