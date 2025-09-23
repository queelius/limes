# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

This is a header-only C++ library for numerical integration, differentiation, and solving ordinary differential equations (ODEs). The library uses modern C++ (C++20) and template metaprogramming to provide generic, type-safe numerical algorithms.

## Commands

### Build and Test
```bash
# Build the ODE test executable
g++ -O3 -std=c++2a -Wall -o ode_test include/ode/ode_test.cpp -I.

# Build integration test (if integration_test.cpp exists)
make integration_test
```

### Compilation Settings
- Compiler: g++ with C++20 standard (`-std=c++2a`)
- Optimization: `-O3` for release builds
- Warnings: `-Wall` enabled
- Include path: `-I.` (current directory as include path)

## Architecture

### Namespace Structure
- Primary namespace: `algebraic_integrators` for integration components
- Secondary namespace: `alex::math` for ODE solvers and mathematical utilities

### Core Components

1. **Numerical Integrators** (`include/numerical_integrators/`)
   - `numerical_univariate_integrator.hpp`: Base template for univariate integration with concept modeling
   - `double_exponential_integrator.hpp`: Double exponential quadrature implementation
   - `gauss_quad_integrator.hpp`: Gaussian quadrature methods
   - `simpson_integrator.hpp`: Simpson's rule integration
   - `rectangle_rule_integrator.hpp`: Rectangle rule integration
   - `line_integral.hpp`: Line integral computations

2. **ODE Solvers** (`include/ode/`)
   - `euler_ode1.hpp`: First-order Euler method for ODEs
   - `euler_ode2.hpp`: Second-order Euler method
   - `rk4_ode1.hpp`: Fourth-order Runge-Kutta for first-order ODEs
   - `rk4_ode2.hpp`: Fourth-order Runge-Kutta for second-order ODEs
   - `interval.hpp`: Interval utilities for ODE solving

3. **Differentiation** (`include/`)
   - `numerical_differentiation.hpp`: General numerical differentiation framework
   - `central_finite_difference.hpp`: Central finite difference methods
   - `antiderivative.hpp`: Antiderivative computations

### Design Patterns

- **Type Erasure**: The library provides type erasure wrappers allowing any object modeling integrator concepts to be wrapped
- **Concept Modeling**: Heavy use of template concepts for compile-time interface checking
- **Value/Result Types**: Each integrator defines `value_type` and `result_type` for consistency
- **Generic Programming**: Templates allow working with arbitrary numeric types and function objects