#pragma once

#include <span>
#include <string>
#include <cstddef>
#include <cassert>
#include "const.hpp"

namespace limes::expr {

// Var<N, T>: A variable reference at position N
// Arity is N+1 (requires at least N+1 arguments to evaluate)
template<std::size_t N, typename T = double>
struct Var {
    using value_type = T;
    static constexpr std::size_t arity_v = N + 1;
    static constexpr std::size_t index_v = N;

    constexpr Var() noexcept = default;

    // Evaluate: return the N-th argument
    [[nodiscard]] constexpr T eval(std::span<T const> args) const noexcept {
        assert(args.size() > N && "Not enough arguments for variable");
        return args[N];
    }

    // Deprecated: use eval() instead
    [[nodiscard]] [[deprecated("use eval() instead")]]
    constexpr T evaluate(std::span<T const> args) const noexcept {
        return eval(args);
    }

    // Derivative: 1 if differentiating with respect to this variable, 0 otherwise
    // Returns One<T>/Zero<T> marker types for compile-time simplification
    template<std::size_t Dim>
    [[nodiscard]] constexpr auto derivative() const noexcept {
        if constexpr (Dim == N) {
            return One<T>{};
        } else {
            return Zero<T>{};
        }
    }

    // String representation
    [[nodiscard]] std::string to_string() const {
        return "x" + std::to_string(N);
    }
};

// Convenient variable aliases for common use cases
template<typename T = double>
inline constexpr Var<0, T> x{};

template<typename T = double>
inline constexpr Var<1, T> y{};

template<typename T = double>
inline constexpr Var<2, T> z{};

// arg<N> for positional variable access
template<std::size_t N, typename T = double>
inline constexpr Var<N, T> arg{};

} // namespace limes::expr
