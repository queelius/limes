# Examples {#examples}

Complete worked examples demonstrating limes features.

## Classic Integrals

### Example 1: The Gaussian Integral

The Gaussian integral is famous: ∫₋∞^∞ e^(-x²) dx = √π

```cpp
#include <limes/limes.hpp>
#include <iostream>
#include <cmath>

int main() {
    using namespace limes::expr;
    using namespace limes::methods;

    auto x = arg<0>;
    auto gaussian = exp(-x * x);

    // Approximate with finite bounds
    auto I = integral(gaussian).over<0>(-10.0, 10.0);

    // Compare methods
    auto r_adaptive = I.eval(adaptive_method(1e-12));
    auto r_gauss = I.eval(gauss<15>());
    auto r_monte = I.eval(monte_carlo_method(1000000).with_seed(42));

    std::cout << "Gaussian integral ∫e^(-x²)dx from -10 to 10:\n";
    std::cout << "  Adaptive:    " << r_adaptive.value() << "\n";
    std::cout << "  Gauss-15:    " << r_gauss.value() << "\n";
    std::cout << "  Monte Carlo: " << r_monte.value() << "\n";
    std::cout << "  Exact (√π):  " << std::sqrt(M_PI) << "\n";
}
```

### Example 2: The Fresnel Integral

The Fresnel sine integral: S(x) = ∫₀ˣ sin(t²) dt

```cpp
#include <limes/limes.hpp>
#include <iostream>

int main() {
    using namespace limes::expr;

    auto t = arg<0>;
    auto fresnel_integrand = sin(t * t);

    // Fresnel S(2) = ∫₀² sin(t²) dt
    auto S = integral(fresnel_integrand).over<0>(0.0, 2.0);
    auto result = S.eval();

    std::cout << "Fresnel S(2) = " << result.value() << "\n";
    std::cout << "Error estimate: " << result.error() << "\n";
    // S(2) ≈ 0.3434
}
```

### Example 3: Elliptic Integral

Complete elliptic integral of the first kind: K(k) = ∫₀^(π/2) 1/√(1-k²sin²θ) dθ

```cpp
#include <limes/limes.hpp>
#include <iostream>
#include <cmath>

int main() {
    using namespace limes::expr;

    auto theta = arg<0>;
    double k = 0.5;  // modulus

    auto integrand = 1.0 / sqrt(1.0 - k*k * sin(theta)*sin(theta));

    auto K = integral(integrand).over<0>(0.0, M_PI / 2);
    auto result = K.eval();

    std::cout << "K(0.5) = " << result.value() << "\n";
    // K(0.5) ≈ 1.6858
}
```

## Multivariate Integration

### Example 4: Double Integral Over a Rectangle

Volume under a paraboloid: ∫₀¹ ∫₀¹ (x² + y²) dx dy = 2/3

```cpp
#include <limes/limes.hpp>
#include <iostream>

int main() {
    using namespace limes::expr;

    auto x = arg<0>;
    auto y = arg<1>;

    auto paraboloid = x*x + y*y;

    // Nested integration
    auto I = integral(paraboloid)
        .over<0>(0.0, 1.0)
        .over<1>(0.0, 1.0);

    auto result = I.eval();
    std::cout << "∫∫(x² + y²) dx dy = " << result.value() << "\n";
    std::cout << "Exact: " << 2.0/3.0 << "\n";
}
```

### Example 5: Triangular Region

Integrate over the triangle 0 ≤ x ≤ 1, 0 ≤ y ≤ x:

```cpp
#include <limes/limes.hpp>
#include <iostream>

int main() {
    using namespace limes::expr;

    auto x = arg<0>;
    auto y = arg<1>;

    // ∫₀¹ ∫₀ˣ xy dy dx
    auto I = integral(x * y)
        .over<1>(0.0, x)     // y from 0 to x (triangular region)
        .over<0>(0.0, 1.0);  // x from 0 to 1

    auto result = I.eval();
    std::cout << "∫∫ xy dy dx over triangle = " << result.value() << "\n";
    std::cout << "Exact: " << 1.0/8.0 << "\n";
    // Result: 1/8 = 0.125
}
```

### Example 6: Monte Carlo Over a Disk

Compute π by integrating over the unit disk:

```cpp
#include <limes/limes.hpp>
#include <iostream>
#include <cmath>

int main() {
    using namespace limes::expr;
    using namespace limes::methods;

    // Area of unit circle = π
    auto I = integral(One<double>{})
        .over_box({{-1.0, 1.0}, {-1.0, 1.0}})
        .where([](double x, double y) {
            return x*x + y*y <= 1.0;
        });

    auto result = I.eval(monte_carlo_method(1000000).with_seed(42));

    std::cout << "Estimated π = " << result.value() << "\n";
    std::cout << "Actual π    = " << M_PI << "\n";
    std::cout << "Error: " << std::abs(result.value() - M_PI) << "\n";
}
```

## Separable Integrals

### Example 7: Product of Independent Integrals

When integrands factor, integration is much faster:

```cpp
#include <limes/limes.hpp>
#include <iostream>
#include <cmath>

int main() {
    using namespace limes::expr;

    auto x = arg<0>;
    auto y = arg<1>;

    // Independent integrals
    auto I_x = integral(exp(-x*x)).over<0>(-5.0, 5.0);  // ≈ √π
    auto I_y = integral(exp(-y*y)).over<1>(-5.0, 5.0);  // ≈ √π

    // Product: evaluates each 1D integral separately
    auto I_product = I_x * I_y;
    auto result = I_product.eval();

    std::cout << "Product integral: " << result.value() << "\n";
    std::cout << "Exact (π):        " << M_PI << "\n";

    // Compare to naive nested integration (much slower for same accuracy)
    auto I_nested = integral(exp(-x*x - y*y))
        .over<0>(-5.0, 5.0)
        .over<1>(-5.0, 5.0);

    auto result_nested = I_nested.eval();
    std::cout << "Nested integral:  " << result_nested.value() << "\n";

    // The product approach: O(n) + O(n) = O(n) evaluations
    // The nested approach:  O(n²) evaluations
}
```

### Example 8: Three-Way Product

Chain multiple independent integrals:

```cpp
#include <limes/limes.hpp>
#include <iostream>

int main() {
    using namespace limes::expr;

    auto x = arg<0>;
    auto y = arg<1>;
    auto z = arg<2>;

    // Three independent 1D integrals
    auto I = integral(x*x).over<0>(0.0, 1.0);     // 1/3
    auto J = integral(y*y).over<1>(0.0, 1.0);     // 1/3
    auto K = integral(z*z).over<2>(0.0, 1.0);     // 1/3

    // Product: I × J × K = (1/3)³ = 1/27
    auto IJK = I * J * K;
    auto result = IJK.eval();

    std::cout << "Product I*J*K = " << result.value() << "\n";
    std::cout << "Exact (1/27): " << 1.0/27.0 << "\n";
}
```

## Differentiation

### Example 9: Gradient of a Scalar Field

Compute the gradient of f(x, y) = x²y + sin(xy):

```cpp
#include <limes/limes.hpp>
#include <iostream>

int main() {
    using namespace limes::expr;

    auto x = var(0, "x");
    auto y = var(1, "y");

    auto f = x*x*y + sin(x*y);

    // Gradient
    auto [df_dx, df_dy] = derivative(f).gradient();

    std::cout << "f = " << f.to_string() << "\n";
    std::cout << "∂f/∂x = " << df_dx.to_string() << "\n";
    std::cout << "∂f/∂y = " << df_dy.to_string() << "\n";

    // Evaluate at (x, y) = (1, 2)
    std::array<double, 2> point{1.0, 2.0};
    std::cout << "\nAt (1, 2):\n";
    std::cout << "f     = " << f.eval(point) << "\n";
    std::cout << "∂f/∂x = " << df_dx.eval(point) << "\n";
    std::cout << "∂f/∂y = " << df_dy.eval(point) << "\n";
}
```

### Example 10: Higher-Order Derivatives

```cpp
#include <limes/limes.hpp>
#include <iostream>

int main() {
    using namespace limes::expr;

    auto x = arg<0>;
    auto f = exp(sin(x));  // e^(sin(x))

    auto df = derivative(f).wrt<0>();        // First derivative
    auto d2f = derivative(f).wrt<0, 0>();    // Second derivative
    auto d3f = derivative(df).wrt<0>().wrt<0>();  // Third derivative (chained)

    std::cout << "f   = " << f.to_string() << "\n";
    std::cout << "f'  = " << df.to_string() << "\n";
    std::cout << "f'' = " << d2f.to_string() << "\n";

    std::array<double, 1> point{0.0};
    std::cout << "\nAt x = 0:\n";
    std::cout << "f(0)   = " << f.eval(point) << "\n";    // e^0 = 1
    std::cout << "f'(0)  = " << df.eval(point) << "\n";   // cos(0)·e^0 = 1
    std::cout << "f''(0) = " << d2f.eval(point) << "\n";  // (cos²(0) - sin(0))·e^0 = 1
}
```

## Fundamental Theorem of Calculus

### Example 11: Verifying FTC

The Fundamental Theorem: ∫ₐˣ f'(t) dt = f(x) - f(a)

```cpp
#include <limes/limes.hpp>
#include <iostream>
#include <cmath>

int main() {
    using namespace limes::expr;

    auto x = arg<0>;
    auto f = sin(x * x);  // f(x) = sin(x²)

    auto df = derivative(f).wrt<0>();  // f'(x) = 2x·cos(x²)

    // Integrate the derivative from 0 to 2
    auto I = integral(df).over<0>(0.0, 2.0);
    auto integral_result = I.eval();

    // Should equal f(2) - f(0)
    std::array<double, 1> at_2{2.0};
    std::array<double, 1> at_0{0.0};
    double ftc_result = f.eval(at_2) - f.eval(at_0);

    std::cout << "∫₀² f'(x) dx = " << integral_result.value() << "\n";
    std::cout << "f(2) - f(0)  = " << ftc_result << "\n";
    std::cout << "Difference:  " << std::abs(integral_result.value() - ftc_result) << "\n";
}
```

## Probability and Statistics

### Example 12: Bayesian Marginal

Computing a marginal distribution from a joint:

```cpp
#include <limes/limes.hpp>
#include <iostream>
#include <cmath>

int main() {
    using namespace limes::expr;

    auto x = arg<0>;
    auto y = arg<1>;

    // Joint PDF: bivariate normal (independent)
    // p(x, y) = (1/2π) exp(-(x² + y²)/2)
    constexpr double two_pi = 2 * M_PI;
    auto joint = exp(-(x*x + y*y) / 2.0) / two_pi;

    // Marginal: p(x) = ∫ p(x, y) dy
    auto marginal = integral(joint).over<1>(-10.0, 10.0);

    // Evaluate marginal at several x values
    std::cout << "Marginal p(x):\n";
    for (double xval : {-2.0, -1.0, 0.0, 1.0, 2.0}) {
        std::array<double, 1> args{xval};
        auto result = marginal.eval(args);
        double exact = std::exp(-xval*xval/2) / std::sqrt(two_pi);
        std::cout << "  p(" << xval << ") = " << result.value()
                  << " (exact: " << exact << ")\n";
    }

    // Verify normalization: ∫ p(x) dx = 1
    auto total = integral(marginal).over<0>(-10.0, 10.0);
    std::cout << "\n∫ p(x) dx = " << total.eval().value() << " (should be 1)\n";
}
```

## Domain Transformations

### Example 13: Removing a Singularity

The integral ∫₀¹ 1/√x dx has a singularity at x = 0.

Using x = t² removes it:

```cpp
#include <limes/limes.hpp>
#include <iostream>
#include <cmath>

int main() {
    using namespace limes::expr;

    auto x = arg<0>;

    // Direct integration (problematic due to singularity)
    auto I_direct = integral(1.0 / sqrt(x)).over<0>(0.001, 1.0);  // Avoid x=0
    auto r_direct = I_direct.eval();
    std::cout << "Direct (avoiding 0): " << r_direct.value() << "\n";

    // Transform: let x = t², then dx = 2t dt
    // ∫ 1/√x dx = ∫ 1/t · 2t dt = ∫ 2 dt
    auto I_transformed = integral(1.0 / sqrt(x)).over<0>(0.0, 1.0)
        .transform(
            [](double t) { return t * t; },  // x = t²
            [](double t) { return 2 * t; },  // |dx/dt| = 2t
            0.0, 1.0
        );
    auto r_transformed = I_transformed.eval();
    std::cout << "Transformed: " << r_transformed.value() << "\n";
    std::cout << "Exact:       " << 2.0 << "\n";
}
```

### Example 14: Domain Splitting

Split at discontinuities for better accuracy:

```cpp
#include <limes/limes.hpp>
#include <iostream>
#include <cmath>

int main() {
    using namespace limes::expr;

    auto x = arg<0>;

    // |x| has a corner at x = 0
    // ∫₋₁¹ |x| dx = 1
    auto I = integral(abs(x)).over<0>(-1.0, 1.0);

    // Direct evaluation
    auto r_direct = I.eval();
    std::cout << "Direct: " << r_direct.value() << "\n";

    // Split at the corner
    auto [left, right] = I.split(0.0);
    auto r_left = left.eval();
    auto r_right = right.eval();
    auto r_split = r_left.value() + r_right.value();
    std::cout << "Split:  " << r_split << "\n";
    std::cout << "Exact:  " << 1.0 << "\n";
}
```

## Method Comparison

### Example 15: Comparing Integration Methods

```cpp
#include <limes/limes.hpp>
#include <iostream>
#include <iomanip>
#include <cmath>

int main() {
    using namespace limes::expr;
    using namespace limes::methods;

    auto x = arg<0>;

    // A smooth oscillatory integrand
    auto f = sin(10 * x) * exp(-x);
    auto I = integral(f).over<0>(0.0, 1.0);

    // Compute with various methods
    struct Result {
        const char* name;
        double value;
        double error;
        std::size_t evals;
    };

    std::vector<Result> results;

    auto r1 = I.eval(gauss<7>());
    results.push_back({"Gauss-7", r1.value(), r1.error(), r1.evaluations()});

    auto r2 = I.eval(gauss<15>());
    results.push_back({"Gauss-15", r2.value(), r2.error(), r2.evaluations()});

    auto r3 = I.eval(simpson_method<100>());
    results.push_back({"Simpson-100", r3.value(), r3.error(), r3.evaluations()});

    auto r4 = I.eval(adaptive_method(1e-10));
    results.push_back({"Adaptive", r4.value(), r4.error(), r4.evaluations()});

    auto r5 = I.eval(monte_carlo_method(100000).with_seed(42));
    results.push_back({"Monte Carlo", r5.value(), r5.error(), r5.evaluations()});

    std::cout << std::setprecision(10);
    std::cout << "∫₀¹ sin(10x)e^(-x) dx:\n\n";
    std::cout << std::setw(15) << "Method" << std::setw(16) << "Value"
              << std::setw(14) << "Error Est" << std::setw(10) << "Evals\n";
    std::cout << std::string(55, '-') << "\n";

    for (const auto& r : results) {
        std::cout << std::setw(15) << r.name
                  << std::setw(16) << r.value
                  << std::setw(14) << r.error
                  << std::setw(10) << r.evals << "\n";
    }
}
```

## Physics Applications

### Example 16: Work Done by a Force

Work = ∫ F·ds along a path:

```cpp
#include <limes/limes.hpp>
#include <iostream>
#include <cmath>

int main() {
    using namespace limes::expr;

    auto t = arg<0>;  // Parameter along path

    // Path: r(t) = (t, t²) for t ∈ [0, 1]
    // dr/dt = (1, 2t)

    // Force field: F(x, y) = (y, -x) = (t², -t)
    // F·dr/dt = t²·1 + (-t)·2t = t² - 2t² = -t²

    auto integrand = -t * t;
    auto work = integral(integrand).over<0>(0.0, 1.0);

    std::cout << "Work done: " << work.eval().value() << "\n";
    std::cout << "Exact:     " << -1.0/3.0 << "\n";
}
```

### Example 17: Center of Mass

For a 2D region with density ρ(x,y) = 1:

```cpp
#include <limes/limes.hpp>
#include <iostream>

int main() {
    using namespace limes::expr;
    using namespace limes::methods;

    auto x = arg<0>;
    auto y = arg<1>;

    // Triangular region: 0 ≤ x ≤ 1, 0 ≤ y ≤ x

    // Total mass (area of triangle = 1/2)
    auto mass = integral(One<double>{})
        .over<1>(0.0, x)
        .over<0>(0.0, 1.0);

    // First moment about y-axis
    auto mx = integral(x)
        .over<1>(0.0, x)
        .over<0>(0.0, 1.0);

    // First moment about x-axis
    auto my = integral(y)
        .over<1>(0.0, x)
        .over<0>(0.0, 1.0);

    double M = mass.eval().value();
    double Mx = mx.eval().value();
    double My = my.eval().value();

    std::cout << "Mass: " << M << " (exact: 0.5)\n";
    std::cout << "Center of mass: (" << Mx/M << ", " << My/M << ")\n";
    std::cout << "Exact: (2/3, 1/3)\n";
}
```

These examples demonstrate the breadth of problems limes can handle. The key insight is that by treating expressions as algebraic objects, we gain composability, inspectability, and the ability to choose the right method for each problem.
