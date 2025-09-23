#include <iostream>
#include <iomanip>
#include <cmath>
#include <numbers>
#include "../include/algebraic_integrators.hpp"

using namespace algebraic_integrators;

int main() {
    using T = double;
    constexpr T pi = std::numbers::pi_v<T>;

    std::cout << "=== Algebraic Integrators Simple Demo ===\n\n";

    // 1. Simple polynomial: ∫₀¹ x² dx = 1/3
    {
        std::cout << "1. Polynomial: ∫₀¹ x² dx = 1/3\n";

        auto f = [](T x) { return x * x; };
        T exact = T(1) / T(3);

        auto result = integrate_adaptive(f, T(0), T(1), T(1e-10));

        std::cout << "   Result: " << std::setprecision(12) << result.value()
                  << " (error: " << std::scientific << result.error() << ")\n";
        std::cout << "   Exact:  " << std::fixed << std::setprecision(12) << exact << "\n";
        std::cout << "   Diff:   " << std::scientific
                  << std::abs(result.value() - exact) << "\n\n";
    }

    // 2. Trigonometric: ∫₀^π sin(x) dx = 2
    {
        std::cout << "2. Trigonometric: ∫₀^π sin(x) dx = 2\n";

        auto f = [](T x) { return std::sin(x); };
        T exact = T(2);

        auto result = integrate_adaptive(f, T(0), pi, T(1e-10));

        std::cout << "   Result: " << std::setprecision(12) << result.value()
                  << " (error: " << std::scientific << result.error() << ")\n";
        std::cout << "   Exact:  " << std::fixed << std::setprecision(12) << exact << "\n";
        std::cout << "   Diff:   " << std::scientific
                  << std::abs(result.value() - exact) << "\n\n";
    }

    // 3. Exponential decay: ∫₀^∞ e^(-x) dx = 1
    {
        std::cout << "3. Semi-infinite: ∫₀^∞ e^(-x) dx = 1\n";

        auto f = [](T x) { return std::exp(-x); };
        T exact = T(1);

        // Use large upper bound as approximation
        auto result = integrate_adaptive(f, T(0), T(20), T(1e-10));

        std::cout << "   Result: " << std::setprecision(12) << result.value()
                  << " (error: " << std::scientific << result.error() << ")\n";
        std::cout << "   Exact:  " << std::fixed << std::setprecision(12) << exact << "\n";
        std::cout << "   Diff:   " << std::scientific
                  << std::abs(result.value() - exact) << "\n\n";
    }

    // 4. Test different accumulators
    {
        std::cout << "4. Accumulator comparison\n";

        // Sum many small numbers
        constexpr T small = T(1e-15);
        constexpr std::size_t n = 1000000;

        accumulators::simple_accumulator<T> simple;
        accumulators::kahan_accumulator<T> kahan;
        accumulators::neumaier_accumulator<T> neumaier;

        for (std::size_t i = 0; i < n; ++i) {
            simple += small;
            kahan += small;
            neumaier += small;
        }

        T exact = small * n;

        std::cout << "   Summing " << n << " values of " << small << "\n";
        std::cout << "   Simple:   " << std::scientific << simple()
                  << " (error: " << std::abs(simple() - exact) << ")\n";
        std::cout << "   Kahan:    " << kahan()
                  << " (error: " << std::abs(kahan() - exact) << ")\n";
        std::cout << "   Neumaier: " << neumaier()
                  << " (error: " << std::abs(neumaier() - exact) << ")\n\n";
    }

    // 5. Test quadrature rules directly
    {
        std::cout << "5. Direct quadrature rule test\n";

        auto f = [](T x) { return x * x; };

        // Create a 3-point Gauss-Legendre rule
        quadrature::gauss_legendre<T, 3> gl3;
        quadrature_integrator<T, decltype(gl3)> integrator{gl3};

        auto result = integrator(f, T(-1), T(1));

        T exact = T(2) / T(3);  // ∫_{-1}^{1} x² dx = 2/3

        std::cout << "   ∫_{-1}^{1} x² dx using GL3\n";
        std::cout << "   Result: " << std::setprecision(12) << result.value() << "\n";
        std::cout << "   Exact:  " << exact << "\n";
        std::cout << "   Error:  " << std::scientific
                  << std::abs(result.value() - exact) << "\n\n";
    }

    // 6. Romberg integration
    {
        std::cout << "6. Romberg integration\n";

        auto f = [](T x) { return std::exp(-x * x); };

        auto romberg = make_romberg<T>();
        auto result = romberg(f, T(0), T(2), T(1e-10));

        T exact = std::sqrt(pi) * std::erf(T(2)) / T(2);

        std::cout << "   ∫₀² e^(-x²) dx\n";
        std::cout << "   Result: " << std::setprecision(12) << result.value() << "\n";
        std::cout << "   Exact:  " << exact << "\n";
        std::cout << "   Error:  " << std::scientific
                  << std::abs(result.value() - exact) << "\n\n";
    }

    std::cout << "=== Demo Complete ===\n";

    return 0;
}