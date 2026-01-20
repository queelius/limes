#pragma once

#include <span>
#include <string>
#include <cstddef>
#include <sstream>

namespace limes::expr {

// Forward declarations for marker types
template<typename T> struct Zero;
template<typename T> struct One;

// =============================================================================
// Type traits for marker constants (defined early for use by all types)
// =============================================================================

template<typename E> inline constexpr bool is_zero_v = false;
template<typename T> inline constexpr bool is_zero_v<Zero<T>> = true;

template<typename E> inline constexpr bool is_one_v = false;
template<typename T> inline constexpr bool is_one_v<One<T>> = true;

// =============================================================================
// Zero<T>: Compile-time zero constant for algebraic simplification
// =============================================================================

template<typename T>
struct Zero {
    using value_type = T;
    static constexpr std::size_t arity_v = 0;

    constexpr Zero() noexcept = default;

    [[nodiscard]] constexpr T eval(std::span<T const> /*args*/) const noexcept {
        return T(0);
    }

    // Deprecated: use eval() instead
    [[nodiscard]] [[deprecated("use eval() instead")]]
    constexpr T evaluate(std::span<T const> args) const noexcept {
        return eval(args);
    }

    template<std::size_t Dim>
    [[nodiscard]] constexpr auto derivative() const noexcept {
        return Zero<T>{};
    }

    [[nodiscard]] std::string to_string() const {
        return "0";
    }
};

// =============================================================================
// One<T>: Compile-time one constant for algebraic simplification
// =============================================================================

template<typename T>
struct One {
    using value_type = T;
    static constexpr std::size_t arity_v = 0;

    constexpr One() noexcept = default;

    [[nodiscard]] constexpr T eval(std::span<T const> /*args*/) const noexcept {
        return T(1);
    }

    // Deprecated: use eval() instead
    [[nodiscard]] [[deprecated("use eval() instead")]]
    constexpr T evaluate(std::span<T const> args) const noexcept {
        return eval(args);
    }

    template<std::size_t Dim>
    [[nodiscard]] constexpr auto derivative() const noexcept {
        return Zero<T>{};
    }

    [[nodiscard]] std::string to_string() const {
        return "1";
    }
};

// =============================================================================
// Const<T>: A constant value expression node
// =============================================================================

// Const<T>: A constant value expression node
// Arity 0 - requires no arguments to evaluate
template<typename T>
struct Const {
    using value_type = T;
    static constexpr std::size_t arity_v = 0;

    T value;

    constexpr Const() noexcept : value{} {}
    constexpr explicit Const(T v) noexcept : value{v} {}

    // Evaluate: just return the constant value
    [[nodiscard]] constexpr T eval(std::span<T const> /*args*/) const noexcept {
        return value;
    }

    // Deprecated: use eval() instead
    [[nodiscard]] [[deprecated("use eval() instead")]]
    constexpr T evaluate(std::span<T const> args) const noexcept {
        return eval(args);
    }

    // Derivative: derivative of a constant is zero
    template<std::size_t Dim>
    [[nodiscard]] constexpr auto derivative() const noexcept {
        return Zero<T>{};
    }

    // String representation
    [[nodiscard]] std::string to_string() const {
        std::ostringstream oss;
        oss << value;
        return oss.str();
    }
};

// Deduction guide
template<typename T>
Const(T) -> Const<T>;

// Factory function for cleaner syntax
template<typename T>
[[nodiscard]] constexpr auto constant(T value) noexcept {
    return Const<T>{value};
}

} // namespace limes::expr
