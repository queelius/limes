# Numerical Calculus Library - Design Vision

## Core Philosophy: Hybrid Symbolic-Numeric Computing

The library uses **generic programming** to preserve analytical relationships while providing numerical methods as a fallback. This creates a "best of both worlds" system.

## Analytical Preservation Through Type System

### Current Implementation (Already Works!)

```cpp
// Integrate f numerically to get F
auto F = antideriv(f, integrator);

// Differentiate F - returns original f analytically!
auto f_recovered = deriv(F);  // Returns f, not numerical approximation
```

**How it works:** `antideriv_of<F,N>` stores the original function `f` as a member, so `deriv()` just returns it. Zero computational cost, perfect accuracy.

### Extended Vision

#### 1. **Compose Operations with Type Information**

```cpp
// Each operation wraps the function and preserves metadata
auto f = /* original function */;
auto F = integrate(f);           // F knows it came from integrating f
auto f2 = differentiate(F);      // Returns f analytically (stored in F)
auto F2 = integrate(f2);         // Returns F analytically (same as F)

// But also:
auto g = sin_function{};
auto g_prime = differentiate(g); // Returns cos_function{} analytically!
```

#### 2. **Type-Based Dispatch**

```cpp
// Generic deriv() function
template <typename T>
auto deriv(T func) {
    if constexpr (is_antiderivative_of<T>) {
        return func.get_original_function();  // Analytical
    } else if constexpr (has_analytical_derivative<T>) {
        return func.derivative();             // Analytical (e.g., sin -> cos)
    } else {
        return numerical_derivative(func);    // Numerical fallback
    }
}
```

#### 3. **Multivariate Extension**

```cpp
// 2D integration
auto F = integrate_2d(f, region);  // Stores f and region
auto df_dx = partial_x(F);         // Returns partial of f analytically
auto df_dy = partial_y(F);         // Returns partial of f analytically
```

## Multivariate Integration Design

### Monte Carlo Integration

**When to use:**
- High dimensions (>3-4)
- Complex/non-rectangular domains
- Discontinuous functions
- Need probabilistic error bounds

**Implementation:**
```cpp
template <typename F, int Dim>
class monte_carlo_integrator {
    size_t samples;

    template <typename Region>
    auto operator()(F f, Region region, double tol) {
        // Stratified sampling
        // Quasi-random sequences (Sobol, Halton)
        // Variance reduction techniques
    }
};
```

**Features needed:**
- Quasi-random sequences (Sobol, Halton) for better convergence
- Stratified sampling
- Importance sampling for peaked functions
- Parallel sampling (embarrassingly parallel)
- Variance reduction

### Cubature Rules

**When to use:**
- Low dimensions (2-4)
- Smooth functions
- Rectangular/simple domains
- Need deterministic error bounds

**Implementation:**
```cpp
template <typename F, int Dim>
class cubature_integrator {
    // Tensor product of 1D rules
    quadrature_rule_1d rule;

    // Or specialized multidimensional rules
    // - Genz-Malik (adaptive)
    // - Clenshaw-Curtis (nested)
    // - Sparse grids (high-dim efficient)
};
```

**Strategies:**
1. **Tensor Product** (2-3D): Extend 1D Gauss-Kronrod to multiple dimensions
2. **Sparse Grids** (4-7D): Smolyak construction for curse of dimensionality
3. **Adaptive Subdivision** (any dim): Subdivide high-error regions

### Hybrid Approach

```cpp
auto integrate_adaptive(F f, Region region, double tol) {
    if (dimension <= 3 && is_smooth(f)) {
        return cubature(f, region, tol);      // Deterministic
    } else {
        return monte_carlo(f, region, tol);   // Stochastic
    }
}
```

## Proposed Type Hierarchy

```cpp
// Base: Any callable
template <typename F>
concept Function = requires(F f) {
    { f(/* args */) } -> std::convertible_to</* return type */>;
};

// Wraps a function with metadata about its origin
template <typename F, typename Metadata>
class analytical_function {
    F func;
    Metadata meta;  // How was this function created?
};

// Specifically for antiderivatives
template <typename F, typename Integrator>
class antiderivative : public analytical_function<F, integrator_metadata> {
    // deriv() returns the stored F
};

// Specifically for derivatives
template <typename F>
class derivative : public analytical_function<F, derivative_metadata> {
    // integrate() might return the original F if available
};

// Known analytical functions
struct sin_function {
    auto derivative() const { return cos_function{}; }
    auto antiderivative() const { return neg_cos_function{}; }
};
```

## Implementation Priority

### Phase 1: Multivariate Integration (High Priority)
1. **2D/3D Cubature** - Tensor product Gauss-Kronrod
2. **Monte Carlo** - Basic with Sobol sequences
3. **Adaptive subdivision** - Error-based region splitting
4. **Tests** - Compare against known analytical integrals

### Phase 2: Enhanced Analytical Preservation (Medium Priority)
5. **Concept-based dispatch** - `has_analytical_derivative<T>`
6. **Composition tracking** - Chain rule for composed functions
7. **Known functions library** - sin, cos, exp, log with analytical derivatives
8. **Expression templates** - Lazy evaluation with type preservation

### Phase 3: Advanced Features (Lower Priority)
9. **Sparse grids** - High-dimensional cubature
10. **Variance reduction** - Importance sampling for Monte Carlo
11. **Symbolic differentiation** - Basic AD or symbolic manipulation
12. **GPU acceleration** - CUDA/SYCL for Monte Carlo

## Examples of Analytical Preservation

### Example 1: Fundamental Theorem of Calculus
```cpp
auto f = [](double x) { return x * x; };
auto F = antideriv(f, integrator);

// Differentiate - gets f back analytically
auto f_prime = deriv(F);  // Returns f, not numerical approximation

// Verify
assert(f(5.0) == f_prime(5.0));  // Exact equality!
```

### Example 2: Chain Rule (Future)
```cpp
auto f = sin_function{};
auto g = [](double x) { return 2*x; };
auto composed = f ∘ g;  // sin(2x)

// Derivative uses chain rule analytically
auto d_composed = deriv(composed);  // 2*cos(2x), stored symbolically
```

### Example 3: Multivariate (Future)
```cpp
auto f = [](double x, double y) { return x*x + y*y; };
auto F = integrate_x(f, 0, 1);  // F(y) = 1/3 + y²x evaluated at x=1

// Partial derivative
auto dF_dy = partial_y(F);  // Returns 2y analytically
```

## Benefits of This Approach

1. **Performance**: Analytical operations are instant, no numerical error
2. **Accuracy**: Preserve exact relationships where possible
3. **Composability**: Operations compose naturally through type system
4. **Flexibility**: Numerical fallback always available
5. **Discoverability**: Type system documents analytical relationships
6. **Zero Overhead**: Template metaprogramming, no runtime cost

## Comparison to Other Libraries

| Library | Approach | Strengths | Weaknesses |
|---------|----------|-----------|------------|
| **SymPy** | Pure symbolic | Exact math | Slow, Python |
| **Eigen** | Pure numeric | Fast, battle-tested | No analytical preservation |
| **Boost.Math** | Pure numeric | Comprehensive | No symbolic features |
| **SymEngine** | Symbolic C++ | Fast symbolic | Complex API |
| **This Library** | Hybrid | Best of both worlds | New, needs validation |

## Open Questions

1. **How much symbolic capability?**
   - Start minimal (just preserving operations)
   - Or build full symbolic engine?

2. **Expression templates?**
   - Could use ET to delay evaluation and preserve structure
   - More complexity but more power

3. **Automatic differentiation?**
   - Forward/reverse mode AD for complex functions
   - Dual numbers for exact derivatives

4. **User-defined analytical functions?**
   - How do users specify analytical derivatives?
   - Traits? Specializations? Macros?

## Naming

Consider renaming to clarify hybrid approach:
- `numerical_calculus` ✓ (current proposal)
- `hybrid_calculus`
- `analytical_numerics`
- `calculus_kit`
- `symbolic_numeric_calculus`
