#include <iostream>
#include <iomanip>
#include <cmath>
#include <chrono>
#include <numbers>
#include "../include/calckit.hpp"

using namespace calckit;
using namespace std::chrono;

template<typename T>
void print_result(const char* name, const integration_result<T>& result, T exact = 0) {
    std::cout << std::setw(25) << name << ": "
              << std::setprecision(12) << result.value()
              << " ± " << std::scientific << result.error()
              << " (" << result.evaluations() << " evals)";

    if (exact != 0) {
        T actual_error = std::abs(result.value() - exact);
        std::cout << " | actual error: " << actual_error;
    }
    std::cout << "\n";
}

int main() {
    using T = double;
    constexpr T pi = std::numbers::pi_v<T>;

    std::cout << "=== Algebraic Integrators Demo ===\n\n";

    // 1. Simple integral: ∫₀¹ x² dx = 1/3
    {
        std::cout << "1. Simple polynomial: ∫₀¹ x² dx = 1/3\n";

        auto f = [](T x) { return x * x; };
        T exact = T(1) / T(3);

        // Different methods
        auto r1 = integrate_adaptive(f, T(0), T(1));
        print_result("Adaptive", r1, exact);

        auto r2 = integrate_robust(f, T(0), T(1));
        print_result("Robust", r2, exact);

        auto r3 = integrate_precise(f, T(0), T(1));
        print_result("Precise", r3, exact);

        std::cout << "\n";
    }

    // 2. Trigonometric: ∫₀^π sin(x) dx = 2
    {
        std::cout << "2. Trigonometric: ∫₀^π sin(x) dx = 2\n";

        auto f = [](T x) { return std::sin(x); };
        T exact = T(2);

        auto result = integrate_adaptive(f, T(0), pi);
        print_result("Adaptive", result, exact);

        std::cout << "\n";
    }

    // 3. Gaussian: ∫_{-∞}^{∞} exp(-x²) dx = √π
    {
        std::cout << "3. Gaussian (infinite): ∫_{-∞}^{∞} exp(-x²) dx = √π\n";

        auto f = [](T x) { return std::exp(-x * x); };
        T exact = std::sqrt(pi);

        auto result = integrate<T>::adaptive(
            f, -std::numeric_limits<T>::infinity(),
            std::numeric_limits<T>::infinity(), T(1e-10)
        );
        print_result("Tanh-Sinh", result, exact);

        std::cout << "\n";
    }

    // 4. Singular integrand: ∫₀¹ 1/√x dx = 2
    {
        std::cout << "4. Singular at endpoint: ∫₀¹ 1/√x dx = 2\n";

        auto f = [](T x) { return T(1) / std::sqrt(x + T(1e-15)); };
        T exact = T(2);

        // Using transform to handle singularity
        auto transform = transforms::make_imt<T>(T(0.5));
        auto result = integrate<T>::with_transform(f, T(0), T(1), transform, T(1e-8));
        print_result("With IMT transform", result, exact);

        std::cout << "\n";
    }

    // 5. Oscillatory integrand
    {
        std::cout << "5. Oscillatory: ∫₀^{2π} sin(10x) dx = 0\n";

        auto f = [](T x) { return std::sin(T(10) * x); };
        T exact = T(0);

        auto result = integrate<T>::oscillatory(f, T(0), T(2) * pi, T(10), T(1e-10));
        print_result("Oscillatory", result, exact);

        std::cout << "\n";
    }

    // 6. Monte Carlo for multivariate
    {
        std::cout << "6. Monte Carlo: ∫∫ (x² + y²) over unit square\n";

        auto mc = integrate<T>::monte_carlo(1000000);

        auto f = [](T x, T y) { return x * x + y * y; };
        auto f_multi = [&f](T* begin, T*) { return f(begin[0], begin[1]); };

        std::vector<std::pair<T, T>> bounds = {{0, 1}, {0, 1}};
        auto result = mc.integrate_multivariate(f_multi, bounds);

        T exact = T(2) / T(3);  // ∫₀¹∫₀¹ (x² + y²) dx dy = 2/3
        print_result("Monte Carlo", result, exact);

        std::cout << "\n";
    }

    // 7. Parallel integration benchmark
    {
        std::cout << "7. Parallel integration benchmark\n";

        // Expensive function
        auto f = [](T x) {
            T sum = T(0);
            for (int i = 0; i < 100; ++i) {
                sum += std::sin(x * (i + 1)) / (i + 1);
            }
            return sum;
        };

        // Sequential
        auto start = high_resolution_clock::now();
        auto r1 = integrate_adaptive(f, T(0), T(10), T(1e-6));
        auto seq_time = duration_cast<milliseconds>(high_resolution_clock::now() - start);

        // Parallel
        start = high_resolution_clock::now();
        auto r2 = integrate_parallel(f, T(0), T(10), T(1e-6));
        auto par_time = duration_cast<milliseconds>(high_resolution_clock::now() - start);

        std::cout << "  Sequential: " << seq_time.count() << " ms\n";
        std::cout << "  Parallel:   " << par_time.count() << " ms\n";
        std::cout << "  Speedup:    " << double(seq_time.count()) / par_time.count() << "x\n";
        std::cout << "  Difference: " << std::abs(r1.value() - r2.value()) << "\n";

        std::cout << "\n";
    }

    // 8. Using different accumulators
    {
        std::cout << "8. Accumulator comparison (sum of 1e-16, 1e6 times)\n";

        auto f = [](T x) { return T(1e-16); };
        std::size_t n = 1000000;
        T a = 0, b = n;

        // Simple accumulator
        using simple = accumulators::simple_accumulator<T>;
        quadrature_integrator<T, quadrature::midpoint_rule<T>, simple> i1{
            quadrature::midpoint_rule<T>{}, simple{}
        };
        auto r1 = i1(f, a, b);

        // Kahan accumulator
        using kahan = accumulators::kahan_accumulator<T>;
        quadrature_integrator<T, quadrature::midpoint_rule<T>, kahan> i2{
            quadrature::midpoint_rule<T>{}, kahan{}
        };
        auto r2 = i2(f, a, b);

        // Klein accumulator
        using klein = accumulators::klein_accumulator<T>;
        quadrature_integrator<T, quadrature::midpoint_rule<T>, klein> i3{
            quadrature::midpoint_rule<T>{}, klein{}
        };
        auto r3 = i3(f, a, b);

        T exact = T(1e-16) * n;
        print_result("Simple accumulator", r1, exact);
        print_result("Kahan accumulator", r2, exact);
        print_result("Klein accumulator", r3, exact);

        std::cout << "\n";
    }

    // 9. Custom integrator with builder
    {
        std::cout << "9. Custom integrator with builder pattern\n";

        auto integrator = make_integrator<T>()
            .with_quadrature("gauss-legendre")
            .with_accumulator("neumaier")
            .with_parallel(false)
            .with_tolerance(T(1e-12));

        auto f = [](T x) { return std::exp(-x) * std::cos(x); };
        auto result = integrator.integrate(f, T(0), T(10));

        print_result("Custom integrator", result);

        std::cout << "\n";
    }

    std::cout << "=== Demo Complete ===\n";

    return 0;
}