# Differentiation API

## First Derivatives

### `central_finite_difference`

8th-order accurate central difference method.

```cpp
template <typename F, typename T>
auto central_finite_difference(F f, T x, T h);
```

**Parameters:**
- `f` - Function to differentiate
- `x` - Point to evaluate derivative
- `h` - Step size (typically 1e-3 to 1e-5)

**Returns:** `differentiation_result<T>` with `.difference` and `.evaluations`

**Example:**
```cpp
auto f = [](double x) { return x * x; };
auto result = central_finite_difference(f, 2.0, 0.001);
// result.difference ≈ 4.0 (exact derivative at x=2)
```

## Second Derivatives

### `central_finite_difference_2nd`

Second-order derivative using central differences.

```cpp
auto f = [](double x) { return std::sin(x); };
auto result = central_finite_difference_2nd(f, M_PI/4, 0.001);
// result.difference ≈ -sin(π/4)
```

## Gradients

### `grad`

Compute gradient of multivariable function.

```cpp
template <typename F>
auto grad(F f, const std::vector<double>& x, double h = 1e-5);
```

**Example:**
```cpp
// f(x,y) = x² + y²
auto f = [](const std::vector<double>& v) {
    return v[0]*v[0] + v[1]*v[1];
};

std::vector<double> point = {1.0, 2.0};
auto gradient = grad(f, point);
// gradient = {2.0, 4.0}
```

## Functional Interface

### `deriv`

Create a derivative function:

```cpp
auto f = [](double x) { return x * x * x; };
auto f_prime = deriv(f);

double derivative_at_2 = f_prime(2.0); // ≈ 12.0
```

## Choosing Step Size

**General guidance:**
- Too large: Truncation error dominates
- Too small: Roundoff error dominates
- Sweet spot: h ≈ 1e-3 to 1e-5 for double precision

**Adaptive step size:**
```cpp
// Use smaller h for faster-changing functions
auto h = std::pow(std::numeric_limits<double>::epsilon(), 1.0/3.0);
```

## Result Type

```cpp
template <typename T>
struct differentiation_result {
    T difference;      // Computed derivative
    size_t evaluations; // Function evaluations
};
```
