#pragma once

#include <span>
#include <string>
#include <cstddef>
#include "binary.hpp"

namespace limes::expr {

// Operation tags for unary operations
struct Neg {};

// Forward declaration of Unary for is_negation_v
template<typename Op, typename E> struct Unary;

// Type trait to detect Unary<Neg, E> (negation expressions)
template<typename E> inline constexpr bool is_negation_v = false;
template<typename E> inline constexpr bool is_negation_v<Unary<Neg, E>> = true;

// Helper to extract the inner type from a negation
template<typename E> struct negation_inner { using type = void; };
template<typename E> struct negation_inner<Unary<Neg, E>> { using type = E; };
template<typename E> using negation_inner_t = typename negation_inner<E>::type;

// Unary<Op, E>: A unary operation node wrapping a child expression
// Arity is same as child's arity
template<typename Op, typename E>
struct Unary {
    using value_type = typename E::value_type;
    using op_type = Op;
    using child_type = E;

    static constexpr std::size_t arity_v = E::arity_v;

    E child;

    constexpr explicit Unary(E c) noexcept : child{c} {}

    // Evaluate: compute child and apply operation
    [[nodiscard]] constexpr value_type eval(std::span<value_type const> args) const {
        value_type c_val = child.eval(args);

        if constexpr (std::is_same_v<Op, Neg>) {
            return -c_val;
        }
    }

    // Deprecated: use eval() instead
    [[nodiscard]] [[deprecated("use eval() instead")]]
    constexpr value_type evaluate(std::span<value_type const> args) const {
        return eval(args);
    }

    // Derivative
    // d/dx[-e] = -de
    template<std::size_t Dim>
    [[nodiscard]] constexpr auto derivative() const {
        auto dc = child.template derivative<Dim>();

        if constexpr (std::is_same_v<Op, Neg>) {
            // Use operator- for simplification (-Zero = Zero)
            return -dc;
        }
    }

    // String representation
    [[nodiscard]] std::string to_string() const {
        if constexpr (std::is_same_v<Op, Neg>) {
            return "(- " + child.to_string() + ")";
        }
    }
};

// Unary negation operator with simplification
template<typename E>
    requires is_expr_node_v<E>
[[nodiscard]] constexpr auto operator-(E e) {
    if constexpr (is_zero_v<E>) {
        return e;           // -0 = 0
    } else if constexpr (is_negation_v<E>) {
        return e.child;     // --x = x
    } else {
        return Unary<Neg, E>{e};
    }
}

// Unary plus (identity) operator
template<typename E>
    requires is_expr_node_v<E>
[[nodiscard]] constexpr auto operator+(E e) {
    return e;
}

} // namespace limes::expr
