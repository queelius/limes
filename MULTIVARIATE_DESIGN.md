# Multivariate Integration Design

## Overview

Extend the composable architecture to support multivariate integration with both deterministic (cubature) and stochastic (Monte Carlo) methods.

## Design Goals

1. **Maintain composability** - Reuse accumulators, error estimation, parallel execution
2. **Flexible domains** - Rectangles, simplices, arbitrary regions
3. **Multiple methods** - Cubature for smooth functions, Monte Carlo for high dimensions
4. **Adaptive subdivision** - Focus computation on high-error regions
5. **Type-safe** - Use C++20 concepts for multidimensional functions

## API Design

### Basic Usage

```cpp
// 2D integration over rectangle [0,1] × [0,1]
auto f = [](std::array<double, 2> x) { return x[0] * x[1]; };
auto result = integrate_2d(f, {0, 0}, {1, 1}, 1e-8);

// 3D integration with Monte Carlo
auto g = [](std::array<double, 3> x) { return x[0]*x[0] + x[1]*x[1] + x[2]*x[2]; };
auto result = integrate_mc(g, {0, 0, 0}, {1, 1, 1}, 1e-6, 1000000);

// N-dimensional with automatic method selection
auto result = integrate_nd(f, lower_bounds, upper_bounds, tol);
```

### Region Types

```cpp
// Hyperrectangle (most common)
struct hyperrectangle {
    std::vector<double> lower;
    std::vector<double> upper;
    int dimension() const;
    double volume() const;
};

// Simplex (triangles, tetrahedra)
struct simplex {
    std::vector<std::array<double, N>> vertices;
    int dimension() const;
    double volume() const;
};

// Ball/Sphere
struct ball {
    std::array<double, N> center;
    double radius;
    int dimension() const;
    double volume() const;
};

// Arbitrary region (Monte Carlo rejection sampling)
template <int Dim>
struct arbitrary_region {
    hyperrectangle<Dim> bounding_box;
    std::function<bool(std::array<double, Dim>)> contains;
};
```

### Concepts

```cpp
namespace algebraic_integrators::concepts {

template <typename F, typename T, int Dim>
concept MultivariateFunction = requires(F f, std::array<T, Dim> x) {
    { f(x) } -> std::convertible_to<T>;
};

template <typename R, int Dim>
concept IntegrationRegion = requires(R region) {
    { region.dimension() } -> std::same_as<int>;
    { region.volume() } -> std::floating_point;
    // For sampling
    { region.sample_uniform() } -> std::convertible_to<std::array<double, Dim>>;
};

} // namespace algebraic_integrators::concepts
```

## Cubature Methods

### 1. Tensor Product Rules

Extend 1D quadrature to multiple dimensions:

```cpp
template <typename T, int Dim, typename Accumulator>
class tensor_product_integrator {
    quadrature_rule_1d<T> rule_1d;
    Accumulator acc;

public:
    template <MultivariateFunction<T, Dim> F>
    auto operator()(F f, hyperrectangle<Dim> region, T tol) {
        // Tensor product: n^d nodes for d dimensions with n nodes per dimension
        // Example for 2D with Gauss-Legendre-5:
        // (x_i, y_j) for all i,j in {1,...,5} = 25 points
        return integrate_tensor_product(f, region);
    }
};
```

**Pros:**
- Simple to implement
- Reuses 1D infrastructure
- Exact for polynomial integrands up to certain degree

**Cons:**
- Exponential growth in nodes: O(n^d)
- Only practical for d ≤ 3-4

**Implementation:**
```cpp
// For 2D with n-point 1D rule:
for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
        point = {nodes_1d[i], nodes_1d[j]};
        weight = weights_1d[i] * weights_1d[j];
        acc.add(f(transform(point, region)) * weight * jacobian);
    }
}
```

### 2. Sparse Grids (Smolyak Algorithm)

Reduce curse of dimensionality for smooth functions:

```cpp
template <typename T, int Dim, typename Accumulator>
class sparse_grid_integrator {
    int level;  // Controls accuracy vs. cost tradeoff

public:
    // O(n * log(n)^(d-1)) points instead of O(n^d)
    auto operator()(F f, hyperrectangle<Dim> region, T tol) {
        return smolyak_integration(f, region, level);
    }
};
```

**Pros:**
- Much fewer points than tensor product
- Good for d = 4-10
- Nested rules enable adaptive refinement

**Cons:**
- More complex implementation
- Requires nested 1D rules (Clenshaw-Curtis ideal)

### 3. Adaptive Cubature (Genz-Malik)

Recursively subdivide high-error regions:

```cpp
template <typename T, int Dim, typename Accumulator>
class adaptive_cubature {
    struct region_error {
        hyperrectangle<Dim> region;
        T value;
        T error;
        bool operator<(const region_error& other) const {
            return error < other.error;  // Min-heap by error
        }
    };

    std::priority_queue<region_error> work_queue;

public:
    auto operator()(F f, hyperrectangle<Dim> region, T tol) {
        // 1. Integrate entire region with embedded error estimate
        // 2. If error > tol, subdivide along longest dimension
        // 3. Add subregions to priority queue
        // 4. Process highest-error regions first
        // 5. Stop when total error < tol
    }
};
```

**Algorithm:**
```cpp
auto integrate_adaptive(F f, Region region, T global_tol) {
    // Initial estimate
    auto [value, error] = integrate_region(f, region);
    if (error < global_tol) return value;

    // Subdivide and process
    priority_queue<region_error> queue;
    queue.push({region, value, error});

    T total_value = value;
    T total_error = error;

    while (total_error > global_tol && !queue.empty()) {
        auto worst = queue.top(); queue.pop();

        // Subdivide worst region
        auto [left, right] = subdivide(worst.region);
        auto [v1, e1] = integrate_region(f, left);
        auto [v2, e2] = integrate_region(f, right);

        // Update totals
        total_value += (v1 + v2 - worst.value);
        total_error += (e1 + e2 - worst.error);

        queue.push({left, v1, e1});
        queue.push({right, v2, e2});
    }

    return {total_value, total_error};
}
```

## Monte Carlo Methods

### 1. Basic Monte Carlo

Simple random sampling:

```cpp
template <typename T, int Dim, typename Accumulator>
class monte_carlo_integrator {
    size_t num_samples;
    std::mt19937_64 rng;

public:
    auto operator()(F f, hyperrectangle<Dim> region, T tol) {
        T volume = region.volume();

        // Sample uniformly
        Accumulator sum, sum_sq;
        for (size_t i = 0; i < num_samples; ++i) {
            auto point = sample_uniform(region, rng);
            T value = f(point);
            sum.add(value);
            sum_sq.add(value * value);
        }

        T mean = sum() / num_samples;
        T variance = (sum_sq() / num_samples - mean * mean) / num_samples;
        T error = volume * std::sqrt(variance);

        return {volume * mean, error};
    }
};
```

**Error:** O(1/√N) - probabilistic, independent of dimension!

### 2. Quasi-Monte Carlo

Use low-discrepancy sequences instead of pseudo-random:

```cpp
template <typename T, int Dim>
class sobol_sequence {
    // Van der Corput sequences
    // Sobol sequence generation

public:
    std::array<T, Dim> next() {
        // Generate next point in [0,1]^Dim with low discrepancy
    }
};

template <typename T, int Dim, typename Accumulator>
class quasi_monte_carlo_integrator {
    sobol_sequence<T, Dim> sobol;

public:
    auto operator()(F f, Region region, T tol) {
        // Use Sobol sequence instead of random sampling
        // Better convergence: O(log(N)^d / N) instead of O(1/√N)
    }
};
```

**Pros:**
- Better convergence than random MC for smooth functions
- Still works in high dimensions
- Deterministic error estimation

**Sequences:**
- Sobol sequence (best overall)
- Halton sequence (simple, but degrades in high-d)
- Latin hypercube sampling (stratified)

### 3. Stratified Sampling

Divide domain into strata and sample proportionally:

```cpp
template <typename T, int Dim>
auto stratified_monte_carlo(F f, Region region, size_t samples) {
    // Divide each dimension into k strata
    int strata_per_dim = std::pow(samples, 1.0 / Dim);

    // Sample proportionally from each stratum
    for (auto stratum : stratify(region, strata_per_dim)) {
        // Sample within this stratum
        samples_per_stratum = total_samples / num_strata;
        // ...
    }
}
```

**Variance reduction:** Much better than pure random sampling.

### 4. Importance Sampling

Sample more where |f| is large:

```cpp
template <typename T, int Dim>
auto importance_sampling(F f, Region region, Distribution proposal) {
    // proposal(x) ≈ |f(x)| / ∫|f|
    // Estimate: ∫f = E[f(X)/proposal(X)] where X ~ proposal

    for (size_t i = 0; i < num_samples; ++i) {
        auto x = sample(proposal);
        T weight = f(x) / proposal(x);
        sum.add(weight);
    }
}
```

**Use case:** Functions with peaks or singularities.

## Unified Interface

```cpp
namespace algebraic_integrators {

template <int Dim, typename T = double>
class integrate_nd {
public:
    // Automatic method selection
    template <MultivariateFunction<T, Dim> F>
    static auto adaptive(F f, std::array<T, Dim> lower,
                        std::array<T, Dim> upper, T tol) {
        if constexpr (Dim <= 3) {
            // Use adaptive cubature for low dimensions
            return adaptive_cubature<T, Dim>{}(f, {lower, upper}, tol);
        } else if constexpr (Dim <= 7) {
            // Use sparse grids for medium dimensions
            return sparse_grid<T, Dim>{}(f, {lower, upper}, tol);
        } else {
            // Use quasi-Monte Carlo for high dimensions
            return quasi_monte_carlo<T, Dim>{}(f, {lower, upper}, tol);
        }
    }

    // Explicit method selection
    template <MultivariateFunction<T, Dim> F>
    static auto cubature(F f, Region region, T tol) {
        return adaptive_cubature<T, Dim>{}(f, region, tol);
    }

    template <MultivariateFunction<T, Dim> F>
    static auto monte_carlo(F f, Region region, T tol, size_t samples = 1000000) {
        return monte_carlo_integrator<T, Dim>{}(f, region, tol, samples);
    }

    template <MultivariateFunction<T, Dim> F>
    static auto quasi_monte_carlo(F f, Region region, T tol) {
        return quasi_monte_carlo_integrator<T, Dim>{}(f, region, tol);
    }
};

// Convenience functions
template <typename F, typename T>
auto integrate_2d(F f, std::array<T, 2> lower, std::array<T, 2> upper, T tol) {
    return integrate_nd<2, T>::adaptive(f, lower, upper, tol);
}

template <typename F, typename T>
auto integrate_3d(F f, std::array<T, 3> lower, std::array<T, 3> upper, T tol) {
    return integrate_nd<3, T>::adaptive(f, lower, upper, tol);
}

} // namespace algebraic_integrators
```

## Implementation Phases

### Phase 1: Basic 2D/3D Cubature (Week 1-2)
- [ ] Tensor product integrator for 2D
- [ ] Tensor product integrator for 3D
- [ ] Basic adaptive subdivision
- [ ] Tests against known integrals

### Phase 2: Monte Carlo (Week 3)
- [ ] Basic Monte Carlo integrator
- [ ] Sobol sequence generator
- [ ] Quasi-Monte Carlo integrator
- [ ] Comparison benchmarks

### Phase 3: Advanced Cubature (Week 4)
- [ ] Sparse grids (Smolyak)
- [ ] Genz-Malik adaptive algorithm
- [ ] Error estimation improvements

### Phase 4: Advanced Monte Carlo (Week 5)
- [ ] Stratified sampling
- [ ] Importance sampling framework
- [ ] Variance reduction techniques

### Phase 5: Integration & Polish (Week 6)
- [ ] Unified API
- [ ] Parallel execution support
- [ ] Documentation and examples
- [ ] Performance benchmarking

## Test Cases

### Analytical Test Functions

```cpp
// 2D Gaussian (exact integral known)
auto gaussian_2d = [](std::array<double, 2> x) {
    return std::exp(-(x[0]*x[0] + x[1]*x[1]));
};
// Exact: ∫∫ exp(-(x²+y²)) dx dy from -∞ to ∞ = π

// 2D product (separable)
auto separable = [](std::array<double, 2> x) {
    return std::sin(x[0]) * std::cos(x[1]);
};
// Exact: [0,π]×[0,π] = 2 * 0 = 0

// 3D sphere volume
auto sphere = [r](std::array<double, 3> x) {
    return (x[0]*x[0] + x[1]*x[1] + x[2]*x[2] <= r*r) ? 1.0 : 0.0;
};
// Exact: 4πr³/3

// Genz test functions (standard benchmark suite)
// - Oscillatory, product peak, corner peak, Gaussian, continuous, discontinuous
```

## Performance Targets

| Method | Dimension | Function Evals | Error | Notes |
|--------|-----------|----------------|-------|-------|
| Tensor Product | 2 | ~100-1000 | 1e-10 | Smooth functions |
| Adaptive Cubature | 2-3 | ~1000-10000 | 1e-8 | General purpose |
| Sparse Grid | 4-7 | ~10000 | 1e-6 | Smooth only |
| Monte Carlo | any | 1M | 1e-3 | Universal fallback |
| Quasi-MC | any | 100k | 1e-4 | Better than MC |

## References

1. **Cubature:**
   - Genz & Malik (1980) - Adaptive algorithms
   - Smolyak (1963) - Sparse grids
   - Clenshaw & Curtis (1960) - Nested rules

2. **Monte Carlo:**
   - Sobol (1967) - Sobol sequences
   - Halton (1960) - Halton sequences
   - Caflisch (1998) - Monte Carlo and quasi-Monte Carlo methods

3. **Software:**
   - Cuba library (Hahn) - Comprehensive cubature/MC
   - HIntLib - C++ cubature library
   - GNU Scientific Library - Basic Monte Carlo
