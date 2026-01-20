#pragma once

#include <span>
#include <string>
#include <vector>
#include <cstddef>
#include <sstream>
#include "const.hpp"
#include "binary.hpp"

namespace limes::expr {

// =============================================================================
// BoundExpr: Expression with one dimension fixed to a constant value
// =============================================================================

/// BoundExpr<E, Dim, BoundValue>: An expression E with dimension Dim bound to a value
///
/// When evaluated, this expression injects the bound value at position Dim
/// and uses the remaining arguments for other dimensions.
///
/// Example:
///   auto f = x * y;                    // f(x, y) = x * y, arity = 2
///   auto g = bind<0>(f, 3.0);          // g(y) = 3 * y, arity = 1
///   g.evaluate({4.0}) == 12.0          // 3 * 4 = 12
///
template<typename E, std::size_t Dim, typename BoundValue>
struct BoundExpr {
    using value_type = typename E::value_type;
    using expr_type = E;
    using bound_type = BoundValue;

    static constexpr std::size_t dim_v = Dim;
    static constexpr std::size_t original_arity_v = E::arity_v;

    // Arity decreases by 1 (bound dimension is consumed)
    // But never goes below 0
    static constexpr std::size_t arity_v = (E::arity_v > 0) ? (E::arity_v - 1) : 0;

    E expr;
    BoundValue bound_value;

    constexpr BoundExpr(E e, BoundValue v) noexcept
        : expr{e}, bound_value{v} {}

    /// Evaluate by injecting bound_value at position Dim
    [[nodiscard]] constexpr value_type eval(std::span<value_type const> args) const {
        // Build full argument vector with bound value at Dim
        std::vector<value_type> full_args;
        full_args.reserve(original_arity_v);

        std::size_t args_idx = 0;
        for (std::size_t i = 0; i < original_arity_v; ++i) {
            if (i == Dim) {
                // Insert the bound value
                full_args.push_back(get_bound_value());
            } else {
                // Take from input args (shifted)
                if (args_idx < args.size()) {
                    full_args.push_back(args[args_idx++]);
                } else {
                    full_args.push_back(value_type{0});
                }
            }
        }

        return expr.eval(std::span<value_type const>{full_args});
    }

    /// Convenience: eval with no arguments (for fully bound expressions)
    [[nodiscard]] constexpr value_type eval() const {
        return eval(std::span<value_type const>{});
    }

    /// Deprecated: use eval() instead
    [[nodiscard]] [[deprecated("use eval() instead")]]
    constexpr value_type evaluate(std::span<value_type const> args) const {
        return eval(args);
    }

    /// Deprecated: use eval() instead
    [[nodiscard]] [[deprecated("use eval() instead")]]
    constexpr value_type evaluate() const {
        return eval();
    }

    /// Derivative with respect to a dimension
    /// Note: The derivative dimension must account for the shift caused by binding
    template<std::size_t DerivDim>
    [[nodiscard]] constexpr auto derivative() const {
        // Map DerivDim to the original dimension in the underlying expression
        // If DerivDim >= Dim, we need to look at DerivDim + 1 in the original
        constexpr std::size_t original_dim = (DerivDim >= Dim) ? (DerivDim + 1) : DerivDim;

        // Get derivative of underlying expression
        auto d_expr = expr.template derivative<original_dim>();

        // Return bound version of the derivative
        return BoundExpr<decltype(d_expr), Dim, BoundValue>{d_expr, bound_value};
    }

    /// String representation
    [[nodiscard]] std::string to_string() const {
        std::ostringstream oss;
        oss << "(bind " << Dim << " " << get_bound_value() << " " << expr.to_string() << ")";
        return oss.str();
    }

private:
    /// Get the bound value (handles both constant and expression bounds)
    [[nodiscard]] constexpr value_type get_bound_value() const {
        if constexpr (std::is_arithmetic_v<BoundValue>) {
            return static_cast<value_type>(bound_value);
        } else {
            // Assume it's a Const<T> or similar
            return bound_value.value;
        }
    }
};

// =============================================================================
// bind() factory function
// =============================================================================

/// Bind dimension Dim of expression E to a constant value
///
/// Example:
///   auto f = x * y;            // f(x, y)
///   auto g = bind<0>(f, 3.0);  // g(y) = 3 * y
///   auto h = bind<1>(f, 4.0);  // h(x) = x * 4
///
template<std::size_t Dim, typename E, typename T>
    requires is_expr_node_v<E>
[[nodiscard]] constexpr auto bind(E expr, T value) {
    static_assert(Dim < E::arity_v || E::arity_v == 0,
        "Cannot bind a dimension that doesn't exist in the expression");

    using ValueType = typename E::value_type;

    if constexpr (std::is_arithmetic_v<T>) {
        // Store as the expression's value type
        return BoundExpr<E, Dim, Const<ValueType>>{expr, Const<ValueType>{static_cast<ValueType>(value)}};
    } else {
        // Assume it's already a Const or expression type
        return BoundExpr<E, Dim, T>{expr, value};
    }
}

/// Bind multiple dimensions at once
/// bind_all<D1, D2, ...>(expr, v1, v2, ...) binds D1 to v1, D2 to v2, etc.
/// Dimensions are bound in order, so later bindings see the reduced arity

/// Helper for binding multiple dimensions
template<std::size_t Dim, std::size_t... Dims, typename E, typename T, typename... Ts>
    requires is_expr_node_v<E>
[[nodiscard]] constexpr auto bind_all(E expr, T value, Ts... values) {
    auto bound_once = bind<Dim>(expr, value);
    if constexpr (sizeof...(Dims) == 0) {
        return bound_once;
    } else {
        return bind_all<Dims...>(bound_once, values...);
    }
}

// =============================================================================
// Partial application: curry-style binding
// =============================================================================

/// Partially apply an expression from the left (bind dimension 0)
template<typename E, typename T>
    requires is_expr_node_v<E>
[[nodiscard]] constexpr auto partial(E expr, T value) {
    return bind<0>(expr, value);
}

/// Partially apply an expression from the right (bind last dimension)
template<typename E, typename T>
    requires is_expr_node_v<E>
[[nodiscard]] constexpr auto partial_right(E expr, T value) {
    constexpr std::size_t last_dim = (E::arity_v > 0) ? (E::arity_v - 1) : 0;
    return bind<last_dim>(expr, value);
}

} // namespace limes::expr
