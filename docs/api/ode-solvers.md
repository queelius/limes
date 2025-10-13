# ODE Solvers API

## First-Order ODEs

Solve equations of the form: dy/dt = f(t, y)

### `euler_ode1`

First-order Euler method (simple, less accurate).

```cpp
template <typename Acc>
class euler_ode1;
```

**Example:**
```cpp
using Acc = simple_accumulator<double>;
euler_ode1<Acc> solver;

// Solve dy/dt = -y with y(0) = 1
auto f = [](double t, double y) { return -y; };
auto result = solver(f, 0.0, 1.0, 1.0, 0.01);
// result.value ≈ exp(-1) ≈ 0.368
```

### `rk4_ode1`

Fourth-order Runge-Kutta (more accurate).

```cpp
rk4_ode1<Acc> solver;
auto result = solver(f, t0, y0, t_final, step_size);
```

**Parameters:**
- `f` - Function f(t, y) defining the ODE
- `t0` - Initial time
- `y0` - Initial value y(t0)
- `t_final` - Final time
- `step_size` - Integration step

## Second-Order ODEs

Solve equations of the form: d²y/dt² = f(t, y, dy/dt)

### `euler_ode2`

```cpp
euler_ode2<Acc> solver;

// Solve simple harmonic oscillator: d²y/dt² = -y
auto f = [](double t, double y, double dy) { return -y; };
auto result = solver(f, 0.0, 1.0, 0.0, 1.0, 0.01);
```

### `rk4_ode2`

```cpp
rk4_ode2<Acc> solver;
auto result = solver(f, t0, y0, dy0, t_final, step_size);
```

**Parameters for second-order:**
- `f` - Function f(t, y, dy/dt)
- `t0` - Initial time
- `y0` - Initial position y(t0)
- `dy0` - Initial velocity dy/dt(t0)
- `t_final` - Final time
- `step_size` - Integration step

## Examples

### Exponential Decay

```cpp
// dy/dt = -k*y
auto decay = [k](double t, double y) { return -k * y; };
rk4_ode1<simple_accumulator<double>> solver;
auto result = solver(decay, 0.0, 100.0, 10.0, 0.1);
```

### Simple Harmonic Oscillator

```cpp
// d²y/dt² = -ω²y
auto sho = [omega](double t, double y, double dy) {
    return -omega * omega * y;
};
rk4_ode2<kahan_accumulator<double>> solver;
auto result = solver(sho, 0.0, 1.0, 0.0, 10.0, 0.01);
```

### Logistic Growth

```cpp
// dy/dt = r*y*(1 - y/K)
auto logistic = [r, K](double t, double y) {
    return r * y * (1.0 - y/K);
};
```

## Choosing Step Size

- Smaller step: More accurate, slower
- Larger step: Faster, may be unstable
- Rule of thumb: Start with h = 0.01, adjust based on accuracy needs
- Use RK4 for better accuracy with larger steps

## Result Type

```cpp
template <typename T>
struct ode_result {
    T value;           // Solution at final time
    size_t steps;      // Number of steps taken
};
```
