# Getting Started

## Installation

This is a header-only library. Simply include the main header:

```cpp
#include "calckit.hpp"
```

Or with CMake:

```cmake
find_package(calckit REQUIRED)
target_link_libraries(your_target ai::calckit)
```

## First Integration

```cpp
#include "calckit.hpp"

int main() {
    // Define a function to integrate
    auto f = [](double x) { return x * x; };

    // Integrate from 0 to 1 with tolerance 1e-8
    auto result = integrate_adaptive(f, 0.0, 1.0, 1e-8);

    std::cout << "Integral: " << result.value << std::endl;
    // Output: 0.333333
}
```

## Building the Library

```bash
mkdir build && cd build
cmake ..
make
ctest  # Run tests
```

### CMake Options

- `AI_BUILD_TESTS=ON/OFF` - Build test suite (default: ON)
- `AI_BUILD_EXAMPLES=ON/OFF` - Build examples (default: ON)
- `AI_ENABLE_AVX2=ON/OFF` - Enable AVX2 SIMD (default: OFF)
- `AI_ENABLE_OPENMP=ON/OFF` - Enable OpenMP parallelization (default: OFF)

## Requirements

- C++20 compliant compiler (GCC 10+, Clang 12+, MSVC 2019+)
- CMake 3.20+ (for building tests/examples)
- Google Test (for testing)
