#pragma once

#include <span>
#include <string>
#include <cstddef>
#include <cmath>
#include "binary.hpp"
#include "unary.hpp"

namespace limes::expr {

// Primitive function tags
struct ExpTag {};
struct LogTag {};
struct SinTag {};
struct CosTag {};
struct SqrtTag {};
struct AbsTag {};
struct TanTag {};
struct SinhTag {};
struct CoshTag {};
struct TanhTag {};
struct AsinTag {};
struct AcosTag {};
struct AtanTag {};
struct AsinhTag {};
struct AcoshTag {};
struct AtanhTag {};

// UnaryFunc<Tag, E>: A unary function applied to a child expression
// This is the base template for primitive functions like exp, sin, cos, etc.
template<typename Tag, typename E>
struct UnaryFunc {
    using value_type = typename E::value_type;
    using tag_type = Tag;
    using child_type = E;

    static constexpr std::size_t arity_v = E::arity_v;

    E child;

    constexpr explicit UnaryFunc(E c) noexcept : child{c} {}

    // Evaluate: compute child and apply the function
    [[nodiscard]] constexpr value_type eval(std::span<value_type const> args) const {
        value_type c_val = child.eval(args);

        if constexpr (std::is_same_v<Tag, ExpTag>) {
            return std::exp(c_val);
        } else if constexpr (std::is_same_v<Tag, LogTag>) {
            return std::log(c_val);
        } else if constexpr (std::is_same_v<Tag, SinTag>) {
            return std::sin(c_val);
        } else if constexpr (std::is_same_v<Tag, CosTag>) {
            return std::cos(c_val);
        } else if constexpr (std::is_same_v<Tag, SqrtTag>) {
            return std::sqrt(c_val);
        } else if constexpr (std::is_same_v<Tag, AbsTag>) {
            return std::abs(c_val);
        } else if constexpr (std::is_same_v<Tag, TanTag>) {
            return std::tan(c_val);
        } else if constexpr (std::is_same_v<Tag, SinhTag>) {
            return std::sinh(c_val);
        } else if constexpr (std::is_same_v<Tag, CoshTag>) {
            return std::cosh(c_val);
        } else if constexpr (std::is_same_v<Tag, TanhTag>) {
            return std::tanh(c_val);
        } else if constexpr (std::is_same_v<Tag, AsinTag>) {
            return std::asin(c_val);
        } else if constexpr (std::is_same_v<Tag, AcosTag>) {
            return std::acos(c_val);
        } else if constexpr (std::is_same_v<Tag, AtanTag>) {
            return std::atan(c_val);
        } else if constexpr (std::is_same_v<Tag, AsinhTag>) {
            return std::asinh(c_val);
        } else if constexpr (std::is_same_v<Tag, AcoshTag>) {
            return std::acosh(c_val);
        } else if constexpr (std::is_same_v<Tag, AtanhTag>) {
            return std::atanh(c_val);
        }
    }

    // Deprecated: use eval() instead
    [[nodiscard]] [[deprecated("use eval() instead")]]
    constexpr value_type evaluate(std::span<value_type const> args) const {
        return eval(args);
    }

    // Derivative using chain rule: d/dx[f(g(x))] = f'(g(x)) * g'(x)
    // Uses operator overloads for automatic simplification with Zero/One types
    template<std::size_t Dim>
    [[nodiscard]] constexpr auto derivative() const {
        auto dc = child.template derivative<Dim>();

        if constexpr (std::is_same_v<Tag, ExpTag>) {
            // d/dx[exp(u)] = exp(u) * du
            auto this_copy = *this;
            return this_copy * dc;  // Uses operator* for simplification
        } else if constexpr (std::is_same_v<Tag, LogTag>) {
            // d/dx[log(u)] = (1/u) * du
            auto one = One<value_type>{};
            auto recip = one / child;
            return recip * dc;
        } else if constexpr (std::is_same_v<Tag, SinTag>) {
            // d/dx[sin(u)] = cos(u) * du
            auto cos_child = UnaryFunc<CosTag, E>{child};
            return cos_child * dc;
        } else if constexpr (std::is_same_v<Tag, CosTag>) {
            // d/dx[cos(u)] = -sin(u) * du
            auto sin_child = UnaryFunc<SinTag, E>{child};
            auto neg_sin = -sin_child;
            return neg_sin * dc;
        } else if constexpr (std::is_same_v<Tag, SqrtTag>) {
            // d/dx[sqrt(u)] = (1 / (2*sqrt(u))) * du
            auto two = Const<value_type>{value_type(2)};
            auto sqrt_child = UnaryFunc<SqrtTag, E>{child};
            auto two_sqrt = two * sqrt_child;
            auto one = One<value_type>{};
            auto recip = one / two_sqrt;
            return recip * dc;
        } else if constexpr (std::is_same_v<Tag, AbsTag>) {
            // d/dx[|u|] = sign(u) * du
            // Note: sign is implemented as u / |u|
            auto abs_child = UnaryFunc<AbsTag, E>{child};
            auto sign_child = child / abs_child;
            return sign_child * dc;
        } else if constexpr (std::is_same_v<Tag, TanTag>) {
            // d/dx[tan(u)] = sec²(u) * du = (1/cos²(u)) * du
            auto cos_child = UnaryFunc<CosTag, E>{child};
            auto cos_sq = cos_child * cos_child;
            auto one = One<value_type>{};
            return (one / cos_sq) * dc;
        } else if constexpr (std::is_same_v<Tag, SinhTag>) {
            // d/dx[sinh(u)] = cosh(u) * du
            auto cosh_child = UnaryFunc<CoshTag, E>{child};
            return cosh_child * dc;
        } else if constexpr (std::is_same_v<Tag, CoshTag>) {
            // d/dx[cosh(u)] = sinh(u) * du
            auto sinh_child = UnaryFunc<SinhTag, E>{child};
            return sinh_child * dc;
        } else if constexpr (std::is_same_v<Tag, TanhTag>) {
            // d/dx[tanh(u)] = sech²(u) * du = (1 - tanh²(u)) * du
            auto tanh_child = *this;
            auto tanh_sq = tanh_child * tanh_child;
            auto one = One<value_type>{};
            return (one - tanh_sq) * dc;
        } else if constexpr (std::is_same_v<Tag, AsinTag>) {
            // d/dx[asin(u)] = 1/√(1-u²) * du
            auto one = One<value_type>{};
            auto u_sq = child * child;
            auto denom = UnaryFunc<SqrtTag, decltype(one - u_sq)>{one - u_sq};
            return (one / denom) * dc;
        } else if constexpr (std::is_same_v<Tag, AcosTag>) {
            // d/dx[acos(u)] = -1/√(1-u²) * du
            auto one = One<value_type>{};
            auto u_sq = child * child;
            auto denom = UnaryFunc<SqrtTag, decltype(one - u_sq)>{one - u_sq};
            return -(one / denom) * dc;
        } else if constexpr (std::is_same_v<Tag, AtanTag>) {
            // d/dx[atan(u)] = 1/(1+u²) * du
            auto one = One<value_type>{};
            auto u_sq = child * child;
            return (one / (one + u_sq)) * dc;
        } else if constexpr (std::is_same_v<Tag, AsinhTag>) {
            // d/dx[asinh(u)] = 1/√(1+u²) * du
            auto one = One<value_type>{};
            auto u_sq = child * child;
            auto denom = UnaryFunc<SqrtTag, decltype(one + u_sq)>{one + u_sq};
            return (one / denom) * dc;
        } else if constexpr (std::is_same_v<Tag, AcoshTag>) {
            // d/dx[acosh(u)] = 1/√(u²-1) * du
            auto one = One<value_type>{};
            auto u_sq = child * child;
            auto denom = UnaryFunc<SqrtTag, decltype(u_sq - one)>{u_sq - one};
            return (one / denom) * dc;
        } else if constexpr (std::is_same_v<Tag, AtanhTag>) {
            // d/dx[atanh(u)] = 1/(1-u²) * du
            auto one = One<value_type>{};
            auto u_sq = child * child;
            return (one / (one - u_sq)) * dc;
        }
    }

    // String representation
    [[nodiscard]] std::string to_string() const {
        std::string func_name;
        if constexpr (std::is_same_v<Tag, ExpTag>) {
            func_name = "exp";
        } else if constexpr (std::is_same_v<Tag, LogTag>) {
            func_name = "log";
        } else if constexpr (std::is_same_v<Tag, SinTag>) {
            func_name = "sin";
        } else if constexpr (std::is_same_v<Tag, CosTag>) {
            func_name = "cos";
        } else if constexpr (std::is_same_v<Tag, SqrtTag>) {
            func_name = "sqrt";
        } else if constexpr (std::is_same_v<Tag, AbsTag>) {
            func_name = "abs";
        } else if constexpr (std::is_same_v<Tag, TanTag>) {
            func_name = "tan";
        } else if constexpr (std::is_same_v<Tag, SinhTag>) {
            func_name = "sinh";
        } else if constexpr (std::is_same_v<Tag, CoshTag>) {
            func_name = "cosh";
        } else if constexpr (std::is_same_v<Tag, TanhTag>) {
            func_name = "tanh";
        } else if constexpr (std::is_same_v<Tag, AsinTag>) {
            func_name = "asin";
        } else if constexpr (std::is_same_v<Tag, AcosTag>) {
            func_name = "acos";
        } else if constexpr (std::is_same_v<Tag, AtanTag>) {
            func_name = "atan";
        } else if constexpr (std::is_same_v<Tag, AsinhTag>) {
            func_name = "asinh";
        } else if constexpr (std::is_same_v<Tag, AcoshTag>) {
            func_name = "acosh";
        } else if constexpr (std::is_same_v<Tag, AtanhTag>) {
            func_name = "atanh";
        }
        return "(" + func_name + " " + child.to_string() + ")";
    }
};

// Convenience functions for creating primitive function expressions

template<typename E>
    requires is_expr_node_v<E>
[[nodiscard]] constexpr auto exp(E e) {
    return UnaryFunc<ExpTag, E>{e};
}

template<typename E>
    requires is_expr_node_v<E>
[[nodiscard]] constexpr auto log(E e) {
    return UnaryFunc<LogTag, E>{e};
}

template<typename E>
    requires is_expr_node_v<E>
[[nodiscard]] constexpr auto sin(E e) {
    return UnaryFunc<SinTag, E>{e};
}

template<typename E>
    requires is_expr_node_v<E>
[[nodiscard]] constexpr auto cos(E e) {
    return UnaryFunc<CosTag, E>{e};
}

template<typename E>
    requires is_expr_node_v<E>
[[nodiscard]] constexpr auto sqrt(E e) {
    return UnaryFunc<SqrtTag, E>{e};
}

template<typename E>
    requires is_expr_node_v<E>
[[nodiscard]] constexpr auto abs(E e) {
    return UnaryFunc<AbsTag, E>{e};
}

template<typename E>
    requires is_expr_node_v<E>
[[nodiscard]] constexpr auto tan(E e) {
    return UnaryFunc<TanTag, E>{e};
}

template<typename E>
    requires is_expr_node_v<E>
[[nodiscard]] constexpr auto sinh(E e) {
    return UnaryFunc<SinhTag, E>{e};
}

template<typename E>
    requires is_expr_node_v<E>
[[nodiscard]] constexpr auto cosh(E e) {
    return UnaryFunc<CoshTag, E>{e};
}

template<typename E>
    requires is_expr_node_v<E>
[[nodiscard]] constexpr auto tanh(E e) {
    return UnaryFunc<TanhTag, E>{e};
}

template<typename E>
    requires is_expr_node_v<E>
[[nodiscard]] constexpr auto asin(E e) {
    return UnaryFunc<AsinTag, E>{e};
}

template<typename E>
    requires is_expr_node_v<E>
[[nodiscard]] constexpr auto acos(E e) {
    return UnaryFunc<AcosTag, E>{e};
}

template<typename E>
    requires is_expr_node_v<E>
[[nodiscard]] constexpr auto atan(E e) {
    return UnaryFunc<AtanTag, E>{e};
}

template<typename E>
    requires is_expr_node_v<E>
[[nodiscard]] constexpr auto asinh(E e) {
    return UnaryFunc<AsinhTag, E>{e};
}

template<typename E>
    requires is_expr_node_v<E>
[[nodiscard]] constexpr auto acosh(E e) {
    return UnaryFunc<AcoshTag, E>{e};
}

template<typename E>
    requires is_expr_node_v<E>
[[nodiscard]] constexpr auto atanh(E e) {
    return UnaryFunc<AtanhTag, E>{e};
}

// Type aliases for common use
template<typename E>
using Exp = UnaryFunc<ExpTag, E>;

template<typename E>
using Log = UnaryFunc<LogTag, E>;

template<typename E>
using Sin = UnaryFunc<SinTag, E>;

template<typename E>
using Cos = UnaryFunc<CosTag, E>;

template<typename E>
using Sqrt = UnaryFunc<SqrtTag, E>;

template<typename E>
using Abs = UnaryFunc<AbsTag, E>;

template<typename E>
using Tan = UnaryFunc<TanTag, E>;

template<typename E>
using Sinh = UnaryFunc<SinhTag, E>;

template<typename E>
using Cosh = UnaryFunc<CoshTag, E>;

template<typename E>
using Tanh = UnaryFunc<TanhTag, E>;

template<typename E>
using Asin = UnaryFunc<AsinTag, E>;

template<typename E>
using Acos = UnaryFunc<AcosTag, E>;

template<typename E>
using Atan = UnaryFunc<AtanTag, E>;

template<typename E>
using Asinh = UnaryFunc<AsinhTag, E>;

template<typename E>
using Acosh = UnaryFunc<AcoshTag, E>;

template<typename E>
using Atanh = UnaryFunc<AtanhTag, E>;

} // namespace limes::expr
