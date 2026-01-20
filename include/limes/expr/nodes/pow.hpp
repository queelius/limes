#pragma once

#include <span>
#include <string>
#include <cstddef>
#include <cmath>
#include "binary.hpp"

namespace limes::expr {

// Forward declaration for is_pow_v
template<typename E, int N> struct Pow;

// Type trait for power detection
template<typename E> inline constexpr bool is_pow_v = false;
template<typename E, int N> inline constexpr bool is_pow_v<Pow<E, N>> = true;

// Helper to extract base type from Pow
template<typename E> struct pow_base { using type = void; };
template<typename E, int N> struct pow_base<Pow<E, N>> { using type = E; };
template<typename E> using pow_base_t = typename pow_base<E>::type;

// Helper to extract exponent from Pow
template<typename E> inline constexpr int pow_exponent_v = 0;
template<typename E, int N> inline constexpr int pow_exponent_v<Pow<E, N>> = N;

// Pow<E, N>: Expression raised to compile-time integer power
// Represents E^N where N is a compile-time constant
template<typename E, int N>
struct Pow {
    using value_type = typename E::value_type;
    using base_type = E;
    static constexpr std::size_t arity_v = E::arity_v;
    static constexpr int exponent = N;

    E base;

    constexpr explicit Pow(E b) noexcept : base{b} {}

    [[nodiscard]] constexpr value_type eval(std::span<value_type const> args) const {
        if constexpr (N == 0) {
            return value_type(1);
        } else if constexpr (N == 1) {
            return base.eval(args);
        } else if constexpr (N == 2) {
            auto v = base.eval(args);
            return v * v;
        } else if constexpr (N == 3) {
            auto v = base.eval(args);
            return v * v * v;
        } else if constexpr (N > 0) {
            return std::pow(base.eval(args), N);
        } else {
            // Negative exponent
            return std::pow(base.eval(args), N);
        }
    }

    // Deprecated: use eval() instead
    [[nodiscard]] [[deprecated("use eval() instead")]]
    constexpr value_type evaluate(std::span<value_type const> args) const {
        return eval(args);
    }

    // Power rule: d/dx[u^n] = n * u^(n-1) * du
    template<std::size_t Dim>
    [[nodiscard]] constexpr auto derivative() const {
        auto du = base.template derivative<Dim>();

        if constexpr (N == 0) {
            return Zero<value_type>{};           // d/dx[x^0] = 0
        } else if constexpr (N == 1) {
            return du;                           // d/dx[x^1] = dx
        } else if constexpr (N == 2) {
            // d/dx[x^2] = 2x * dx
            return Const<value_type>{value_type(2)} * base * du;
        } else {
            // d/dx[x^n] = n * x^(n-1) * dx
            auto coeff = Const<value_type>{value_type(N)};
            auto reduced = Pow<E, N-1>{base};
            return coeff * reduced * du;
        }
    }

    [[nodiscard]] std::string to_string() const {
        return "(^ " + base.to_string() + " " + std::to_string(N) + ")";
    }
};

// Factory function: pow<N>(expr)
template<int N, typename E>
    requires is_expr_node_v<E>
[[nodiscard]] constexpr auto pow(E e) {
    if constexpr (N == 0) {
        return One<typename E::value_type>{};
    } else if constexpr (N == 1) {
        return e;
    } else {
        return Pow<E, N>{e};
    }
}

// Convenience: square(expr) = pow<2>(expr)
template<typename E>
    requires is_expr_node_v<E>
[[nodiscard]] constexpr auto square(E e) {
    return Pow<E, 2>{e};
}

// Convenience: cube(expr) = pow<3>(expr)
template<typename E>
    requires is_expr_node_v<E>
[[nodiscard]] constexpr auto cube(E e) {
    return Pow<E, 3>{e};
}

} // namespace limes::expr
