#pragma once

#include <span>
#include <string>
#include <cstddef>
#include <algorithm>
#include "../concepts.hpp"

namespace limes::expr {

// Forward declarations
template<typename T> struct Const;
template<typename T> struct Zero;
template<typename T> struct One;

// =============================================================================
// Conditional<Cond, Then, Else>: Piecewise/conditional expressions
// =============================================================================

// Conditional expression: if (condition > 0) then then_branch else else_branch
// This enables piecewise functions like Heaviside, ReLU, etc.
template<typename Cond, typename Then, typename Else>
struct Conditional {
    using value_type = typename Then::value_type;
    using condition_type = Cond;
    using then_type = Then;
    using else_type = Else;

    static constexpr std::size_t arity_v =
        std::max({Cond::arity_v, Then::arity_v, Else::arity_v});

    Cond condition;
    Then then_branch;
    Else else_branch;

    constexpr Conditional(Cond c, Then t, Else e) noexcept
        : condition{c}, then_branch{t}, else_branch{e} {}

    // Evaluate: check condition and return appropriate branch
    [[nodiscard]] constexpr value_type eval(std::span<value_type const> args) const {
        value_type cond_val = condition.eval(args);
        if (cond_val > value_type(0)) {
            return then_branch.eval(args);
        } else {
            return else_branch.eval(args);
        }
    }

    // Deprecated: use eval() instead
    [[nodiscard]] [[deprecated("use eval() instead")]]
    constexpr value_type evaluate(std::span<value_type const> args) const {
        return eval(args);
    }

    // Derivative using indicator function approach
    // d/dx[if c>0 then f else g] = indicator(c>0) * df + indicator(c<=0) * dg
    // At the boundary (c=0), we use the subgradient convention:
    // We treat it as: (df if c>0, dg otherwise)
    template<std::size_t Dim>
    [[nodiscard]] constexpr auto derivative() const {
        auto df = then_branch.template derivative<Dim>();
        auto dg = else_branch.template derivative<Dim>();

        // Return a conditional derivative
        return Conditional<Cond, decltype(df), decltype(dg)>{condition, df, dg};
    }

    // String representation (Lisp-style conditional)
    [[nodiscard]] std::string to_string() const {
        return "(if " + condition.to_string() + " "
             + then_branch.to_string() + " "
             + else_branch.to_string() + ")";
    }
};

// =============================================================================
// Type traits for Conditional
// =============================================================================

template<typename T>
struct is_conditional : std::false_type {};

template<typename C, typename T, typename E>
struct is_conditional<Conditional<C, T, E>> : std::true_type {};

template<typename T>
inline constexpr bool is_conditional_v = is_conditional<T>::value;

// =============================================================================
// Factory functions
// =============================================================================

// if_then_else(condition, then_expr, else_expr)
// Returns then_expr if condition > 0, else_expr otherwise
template<typename C, typename T, typename E>
    requires (is_expr_node_v<C> && is_expr_node_v<T> && is_expr_node_v<E>)
[[nodiscard]] constexpr auto if_then_else(C cond, T then_expr, E else_expr) {
    return Conditional<C, T, E>{cond, then_expr, else_expr};
}

// Overload: condition is expression, then/else are scalars
template<typename C, typename T>
    requires (is_expr_node_v<C> && std::is_arithmetic_v<T>)
[[nodiscard]] constexpr auto if_then_else(C cond, T then_val, T else_val) {
    using VT = typename C::value_type;
    return Conditional<C, Const<VT>, Const<VT>>{
        cond,
        Const<VT>{static_cast<VT>(then_val)},
        Const<VT>{static_cast<VT>(else_val)}
    };
}

// =============================================================================
// Convenience functions for common piecewise functions
// =============================================================================

// heaviside(e): H(x) = 1 if x > 0, 0 otherwise
// Note: At x=0, returns 0 (convention choice)
template<typename E>
    requires is_expr_node_v<E>
[[nodiscard]] constexpr auto heaviside(E e) {
    using T = typename E::value_type;
    return if_then_else(e, One<T>{}, Zero<T>{});
}

// ramp(e): max(e, 0) - the positive part function
// This is equivalent to ReLU in neural networks
template<typename E>
    requires is_expr_node_v<E>
[[nodiscard]] constexpr auto ramp(E e) {
    using T = typename E::value_type;
    return if_then_else(e, e, Zero<T>{});
}

// sign(e): sign(x) = 1 if x > 0, -1 if x < 0, 0 if x = 0
template<typename E>
    requires is_expr_node_v<E>
[[nodiscard]] constexpr auto sign(E e) {
    using T = typename E::value_type;
    // sign(x) = if x > 0 then 1 else (if -x > 0 then -1 else 0)
    auto neg_e = -e;
    auto neg_one = Const<T>{T(-1)};
    auto inner = if_then_else(neg_e, neg_one, Zero<T>{});
    return if_then_else(e, One<T>{}, inner);
}

// clamp(e, lo, hi): clamps e to [lo, hi]
template<typename E, typename T>
    requires (is_expr_node_v<E> && std::is_arithmetic_v<T>)
[[nodiscard]] constexpr auto clamp(E e, T lo, T hi) {
    using VT = typename E::value_type;
    auto lo_const = Const<VT>{static_cast<VT>(lo)};
    auto hi_const = Const<VT>{static_cast<VT>(hi)};

    // if (e - lo) > 0 then (if (hi - e) > 0 then e else hi) else lo
    auto above_lo = e - lo_const;
    auto below_hi = hi_const - e;
    auto upper_clamped = if_then_else(below_hi, e, hi_const);
    return if_then_else(above_lo, upper_clamped, lo_const);
}

// indicator(e): 1 if e > 0, 0 otherwise (same as heaviside but emphasizes indicator function semantics)
template<typename E>
    requires is_expr_node_v<E>
[[nodiscard]] constexpr auto indicator(E e) {
    return heaviside(e);
}

} // namespace limes::expr
