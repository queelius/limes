# Algebraic Integrators

Modern C++20 header-only library for numerical integration, differentiation, and ODE solving.

## Key Features

- **Composable Architecture**: Mix and match accumulators, quadrature rules, and integrators
- **Generic Programming**: Templates and concepts for type-safe numerical algorithms
- **High Performance**: SIMD support (AVX2) and OpenMP parallelization
- **Comprehensive**: Integration, differentiation, ODEs, and antiderivatives
- **Well-Tested**: 118+ unit tests with Google Test

## Quick Example

```cpp
#include "algebraic_integrators.hpp"

int main() {
    // Simple integration
    auto f = [](double x) { return x * x; };
    auto result = integrate_adaptive(f, 0.0, 1.0, 1e-8);

    std::cout << "Integral: " << result.value << std::endl;
    // Output: Integral: 0.333333
}
```

## Components

- **Accumulators**: Simple, Kahan, Neumaier, Klein, Pairwise
- **Quadrature Rules**: Gauss-Legendre, Gauss-Kronrod, Tanh-Sinh, Simpson
- **ODE Solvers**: Euler, Runge-Kutta 4
- **Differentiation**: Central finite difference, gradients
- **Transforms**: Coordinate transforms for infinite intervals and singularities

## Installation

Header-only library - just include the headers:

```bash
#include "algebraic_integrators.hpp"
```

Or with CMake:

```cmake
find_package(algebraic_integrators REQUIRED)
target_link_libraries(your_target ai::algebraic_integrators)
```

## Requirements

- C++20 compliant compiler (GCC 10+, Clang 12+, MSVC 2019+)
- CMake 3.20+ (for building tests/examples)
- Google Test (for testing)
