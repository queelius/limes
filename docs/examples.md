# Examples

## Basic Integration

### Simple Integral

```cpp
#include "calckit.hpp"

int main() {
    // Integrate x² from 0 to 1
    auto f = [](double x) { return x * x; };
    auto result = integrate_adaptive(f, 0.0, 1.0, 1e-10);

    std::cout << "∫x² dx from 0 to 1 = " << result.value << std::endl;
    // Output: 0.333333...
}
```

### Trigonometric Function

```cpp
auto f = [](double x) { return std::sin(x); };
auto result = integrate_adaptive(f, 0.0, M_PI, 1e-10);
// result.value ≈ 2.0
```

## Custom Composition

### High-Precision Integration

```cpp
using Acc = klein_accumulator<double>;
gauss_kronrod_integrator<Acc, 15, 31> integrator;

auto f = [](double x) { return std::exp(-x*x); };
auto result = integrator(f, -5.0, 5.0, 1e-12);
// Very accurate result
```

### Fast Integration

```cpp
using Acc = simple_accumulator<float>;
simpson_univariate_integrator<Acc> integrator;

auto f = [](float x) { return x * x + 1.0f; };
auto result = integrator(f, 0.0f, 10.0f, 1000);
```

## ODE Examples

### Population Growth

```cpp
using Acc = kahan_accumulator<double>;
rk4_ode1<Acc> solver;

// Logistic growth: dP/dt = r*P*(1 - P/K)
double r = 0.1;  // Growth rate
double K = 100.0; // Carrying capacity

auto growth = [r, K](double t, double P) {
    return r * P * (1.0 - P/K);
};

// Start with P(0) = 10
auto result = solver(growth, 0.0, 10.0, 50.0, 0.1);
std::cout << "Population at t=50: " << result.value << std::endl;
```

### Pendulum

```cpp
rk4_ode2<simple_accumulator<double>> solver;

double g = 9.81, L = 1.0;
auto pendulum = [g, L](double t, double theta, double omega) {
    return -(g/L) * std::sin(theta);
};

// Start at 30° with zero velocity
auto result = solver(pendulum, 0.0, M_PI/6, 0.0, 10.0, 0.01);
```

## Differentiation Examples

### Derivative of Polynomial

```cpp
auto f = [](double x) { return x*x*x - 2*x*x + x - 1; };
auto f_prime = deriv(f);

double x = 2.0;
std::cout << "f'(2) = " << f_prime(x) << std::endl;
// Output: 7.0 (exact: 3*4 - 4*2 + 1 = 7)
```

### Gradient Descent

```cpp
// Minimize f(x,y) = (x-3)² + (y-4)²
auto f = [](const std::vector<double>& v) {
    return std::pow(v[0]-3, 2) + std::pow(v[1]-4, 2);
};

std::vector<double> x = {0.0, 0.0};
double learning_rate = 0.1;

for (int i = 0; i < 100; i++) {
    auto gradient = grad(f, x);
    x[0] -= learning_rate * gradient[0];
    x[1] -= learning_rate * gradient[1];
}
// x converges to {3.0, 4.0}
```

## Antiderivatives

### Symbolic Antiderivative

```cpp
using Acc = neumaier_accumulator<double>;
simpson_univariate_integrator<Acc> integrator;

auto f = [](double x) { return 2*x; };
auto F = antideriv(f, integrator);

std::cout << "F(5) = " << F(5.0) << std::endl;
// Output: 25.0 (antiderivative is x²)
```

### Initial Value Problem

```cpp
auto f = [](double x) { return std::sin(x); };
auto F = antideriv(f, integrator);

// Solve y' = sin(x), y(0) = 3
auto y = solve(F, 0.0, 3.0);

std::cout << "y(π) = " << y(M_PI) << std::endl;
// Output: 5.0 (exact: 3 + 2 = 5)
```

## Parallel Integration

```cpp
#include "parallel_integrators.hpp"

// Enable OpenMP in CMake with AI_ENABLE_OPENMP=ON
auto f = [](double x) { return std::exp(-x*x) * std::sin(x); };
auto result = integrate_parallel(f, 0.0, 10.0, 1e-8);
// Automatically uses all CPU cores
```

## Transform Examples

### Infinite Interval

```cpp
// Integrate exp(-x) from 0 to ∞
auto f = [](double x) { return std::exp(-x); };
auto result = integrate_adaptive(
    f, 0.0, std::numeric_limits<double>::infinity(), 1e-8
);
// result.value ≈ 1.0
```

### Singularity at Endpoint

```cpp
// Integrate 1/√x from 0 to 1
tanh_sinh_integrator<kahan_accumulator<double>> integrator;
auto f = [](double x) { return 1.0 / std::sqrt(x); };
auto result = integrator(f, 0.0, 1.0, 1e-8);
// result.value ≈ 2.0
```
