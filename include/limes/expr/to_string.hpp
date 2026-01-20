#pragma once

#include <string>
#include <sstream>
#include "nodes/const.hpp"
#include "nodes/var.hpp"
#include "nodes/binary.hpp"
#include "nodes/unary.hpp"
#include "nodes/primitives.hpp"
#include "integral.hpp"

namespace limes::expr {

// Free function to_string for any expression node
// Uses ADL to call the node's to_string() method
template<typename E>
    requires is_expr_node_v<E>
[[nodiscard]] std::string to_string(E const& expr) {
    return expr.to_string();
}

// Pretty print an expression (indented format)
namespace detail {

template<typename E>
void pretty_print_impl(std::ostream& os, E const& expr, int indent) {
    std::string padding(indent * 2, ' ');
    os << padding << expr.to_string();
}

} // namespace detail

template<typename E>
    requires is_expr_node_v<E>
[[nodiscard]] std::string pretty_print(E const& expr) {
    std::ostringstream oss;
    detail::pretty_print_impl(oss, expr, 0);
    return oss.str();
}

// Expression info: arity, type name, etc.
template<typename E>
    requires is_expr_node_v<E>
struct expr_info {
    static constexpr std::size_t arity = E::arity_v;

    static std::string type_name() {
        // Basic type name extraction (simplified)
        return typeid(E).name();
    }
};

// Stream output operator
template<typename E>
    requires is_expr_node_v<E>
std::ostream& operator<<(std::ostream& os, E const& expr) {
    return os << expr.to_string();
}

} // namespace limes::expr
