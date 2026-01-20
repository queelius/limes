#pragma once

#include <span>
#include <string>
#include <string_view>
#include <cstddef>
#include <cassert>
#include "const.hpp"

namespace limes::expr {

// =============================================================================
// NamedVar<T>: A named variable expression for better debugging and inspection
// =============================================================================

/// NamedVar is like Var<N, T> but carries a name for debugging and to_string()
/// The name is a runtime string_view, allowing for descriptive variable names.
///
/// Example:
///   auto x = var(0, "x");
///   auto y = var(1, "y");
///   auto f = sin(x) * cos(y);
///   f.to_string();  // "(* (sin x) (cos y))" instead of "(* (sin x0) (cos x1))"
///
template<typename T = double>
struct NamedVar {
    using value_type = T;
    static constexpr std::size_t arity_v = 1;  // Set per instance

    std::size_t dim;
    std::string_view name;

    constexpr NamedVar(std::size_t dimension, std::string_view n) noexcept
        : dim{dimension}, name{n} {}

    /// Get the effective arity (dim + 1)
    [[nodiscard]] constexpr std::size_t arity() const noexcept {
        return dim + 1;
    }

    // Evaluate: return the dim-th argument
    [[nodiscard]] constexpr T eval(std::span<T const> args) const noexcept {
        assert(args.size() > dim && "Not enough arguments for named variable");
        return args[dim];
    }

    // Deprecated: use eval() instead
    [[nodiscard]] [[deprecated("use eval() instead")]]
    constexpr T evaluate(std::span<T const> args) const noexcept {
        return eval(args);
    }

    // Derivative: 1 if differentiating with respect to this variable, 0 otherwise
    template<std::size_t Dim>
    [[nodiscard]] constexpr auto derivative() const noexcept {
        if constexpr (Dim == 0) {
            // Assuming this is effectively a Var<0> if called this way
            // We need runtime check
            return (dim == Dim) ? One<T>{} : Zero<T>{};
        } else {
            // Static dimension check not directly possible with runtime dim
            // Return marker based on compile-time dimension
            // Note: This is a limitation - runtime dims can't be checked at compile time
            return Zero<T>{};
        }
    }

    // Runtime derivative with respect to dimension
    [[nodiscard]] constexpr auto derivative_rt(std::size_t d) const noexcept {
        return (d == dim) ? One<T>{} : Zero<T>{};
    }

    // String representation using the name
    [[nodiscard]] std::string to_string() const {
        return std::string(name);
    }

    // Access the dimension index
    [[nodiscard]] constexpr std::size_t dimension() const noexcept {
        return dim;
    }
};

// Type trait for NamedVar detection
template<typename T>
struct is_named_var : std::false_type {};

template<typename T>
struct is_named_var<NamedVar<T>> : std::true_type {};

template<typename T>
inline constexpr bool is_named_var_v = is_named_var<T>::value;

/// Factory function for creating named variables
/// Usage: auto x = var(0, "x");
template<typename T = double>
[[nodiscard]] constexpr NamedVar<T> var(std::size_t dim, std::string_view name) noexcept {
    return NamedVar<T>{dim, name};
}

// =============================================================================
// StaticNamedVar<N, T>: Named variable with compile-time dimension
// =============================================================================

/// StaticNamedVar combines compile-time dimension with runtime name
/// This provides both the type-safety of Var<N, T> and the readability of named vars.
///
/// Example:
///   auto x = named<0>("x");
///   auto y = named<1>("y");
///
template<std::size_t N, typename T = double>
struct StaticNamedVar {
    using value_type = T;
    static constexpr std::size_t arity_v = N + 1;
    static constexpr std::size_t index_v = N;

    std::string_view name;

    constexpr explicit StaticNamedVar(std::string_view n) noexcept : name{n} {}

    // Evaluate: return the N-th argument
    [[nodiscard]] constexpr T eval(std::span<T const> args) const noexcept {
        assert(args.size() > N && "Not enough arguments for named variable");
        return args[N];
    }

    // Deprecated: use eval() instead
    [[nodiscard]] [[deprecated("use eval() instead")]]
    constexpr T evaluate(std::span<T const> args) const noexcept {
        return eval(args);
    }

    // Derivative: 1 if differentiating with respect to this variable, 0 otherwise
    template<std::size_t Dim>
    [[nodiscard]] constexpr auto derivative() const noexcept {
        if constexpr (Dim == N) {
            return One<T>{};
        } else {
            return Zero<T>{};
        }
    }

    // String representation using the name
    [[nodiscard]] std::string to_string() const {
        return std::string(name);
    }

    // Access the dimension index
    [[nodiscard]] static constexpr std::size_t dimension() noexcept {
        return N;
    }
};

// Type trait for StaticNamedVar detection
template<typename T>
struct is_static_named_var : std::false_type {};

template<std::size_t N, typename T>
struct is_static_named_var<StaticNamedVar<N, T>> : std::true_type {};

template<typename T>
inline constexpr bool is_static_named_var_v = is_static_named_var<T>::value;

/// Factory function for creating static named variables
/// Usage: auto x = named<0>("x");
template<std::size_t N, typename T = double>
[[nodiscard]] constexpr StaticNamedVar<N, T> named(std::string_view name) noexcept {
    return StaticNamedVar<N, T>{name};
}

// =============================================================================
// Variable declarations helper
// =============================================================================

/// Helper struct to declare multiple named variables at once
/// Usage:
///   auto [x, y, z] = declare_vars<3>("x", "y", "z");
///
template<std::size_t N, typename T, std::size_t... Is>
constexpr auto declare_vars_impl(std::array<std::string_view, N> const& names,
                                   std::index_sequence<Is...>) {
    return std::make_tuple(StaticNamedVar<Is, T>{names[Is]}...);
}

template<std::size_t N, typename T = double, typename... Names>
    requires (sizeof...(Names) == N)
[[nodiscard]] constexpr auto declare_vars(Names... names) {
    std::array<std::string_view, N> name_array{names...};
    return declare_vars_impl<N, T>(name_array, std::make_index_sequence<N>{});
}

/// Convenience for 2 variables
template<typename T = double>
[[nodiscard]] constexpr auto vars_xy(std::string_view x_name = "x",
                                      std::string_view y_name = "y") {
    return std::make_pair(StaticNamedVar<0, T>{x_name},
                          StaticNamedVar<1, T>{y_name});
}

/// Convenience for 3 variables
template<typename T = double>
[[nodiscard]] constexpr auto vars_xyz(std::string_view x_name = "x",
                                       std::string_view y_name = "y",
                                       std::string_view z_name = "z") {
    return std::make_tuple(StaticNamedVar<0, T>{x_name},
                           StaticNamedVar<1, T>{y_name},
                           StaticNamedVar<2, T>{z_name});
}

} // namespace limes::expr
