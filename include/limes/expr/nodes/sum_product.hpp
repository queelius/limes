#pragma once

#include <span>
#include <string>
#include <cstddef>
#include <algorithm>
#include <vector>
#include <sstream>
#include "../concepts.hpp"

namespace limes::expr {

// Forward declarations
template<typename T> struct Const;
template<typename T> struct Zero;
template<typename T> struct One;

// =============================================================================
// FiniteSum<Expr, IndexDim>: Summation over an integer index
// =============================================================================

// FiniteSum represents Σᵢ₌ₗₒʰⁱ body(x, i)
// The IndexDim specifies which variable dimension is used for the summation index.
// The bounds lo and hi are inclusive: sum over i = lo, lo+1, ..., hi
template<typename Expr, std::size_t IndexDim>
struct FiniteSum {
    using value_type = typename Expr::value_type;
    using body_type = Expr;

    // Arity is one less than body's arity (the index variable is consumed)
    static constexpr std::size_t arity_v =
        (Expr::arity_v > IndexDim) ? (Expr::arity_v - 1) : Expr::arity_v;

    Expr body;
    int lo;
    int hi;

    constexpr FiniteSum(Expr b, int l, int h) noexcept : body{b}, lo{l}, hi{h} {}

    // Evaluate: sum body over index from lo to hi
    [[nodiscard]] constexpr value_type eval(std::span<value_type const> args) const {
        value_type sum = value_type(0);

        // Create extended args array with space for the index
        std::vector<value_type> ext_args(args.begin(), args.end());

        // Insert placeholder for index at IndexDim
        if (ext_args.size() <= IndexDim) {
            ext_args.resize(IndexDim + 1);
        } else {
            ext_args.insert(ext_args.begin() + IndexDim, value_type(0));
        }

        for (int i = lo; i <= hi; ++i) {
            ext_args[IndexDim] = static_cast<value_type>(i);
            sum += body.eval(std::span<value_type const>(ext_args));
        }

        return sum;
    }

    // Deprecated: use eval() instead
    [[nodiscard]] [[deprecated("use eval() instead")]]
    constexpr value_type evaluate(std::span<value_type const> args) const {
        return eval(args);
    }

    // Derivative: d/dx[Σf(x,i)] = Σ(df/dx)
    // The derivative passes through the sum
    template<std::size_t Dim>
    [[nodiscard]] constexpr auto derivative() const {
        // Adjust dimension if it's >= IndexDim (shift up to account for index)
        constexpr std::size_t AdjustedDim = (Dim >= IndexDim) ? (Dim + 1) : Dim;
        auto dbody = body.template derivative<AdjustedDim>();
        return FiniteSum<decltype(dbody), IndexDim>{dbody, lo, hi};
    }

    // String representation
    [[nodiscard]] std::string to_string() const {
        std::ostringstream oss;
        oss << "(sum[i" << IndexDim << "=" << lo << ".." << hi << "] " << body.to_string() << ")";
        return oss.str();
    }
};

// =============================================================================
// FiniteProduct<Expr, IndexDim>: Product over an integer index
// =============================================================================

// FiniteProduct represents Πᵢ₌ₗₒʰⁱ body(x, i)
template<typename Expr, std::size_t IndexDim>
struct FiniteProduct {
    using value_type = typename Expr::value_type;
    using body_type = Expr;

    // Arity is one less than body's arity (the index variable is consumed)
    static constexpr std::size_t arity_v =
        (Expr::arity_v > IndexDim) ? (Expr::arity_v - 1) : Expr::arity_v;

    Expr body;
    int lo;
    int hi;

    constexpr FiniteProduct(Expr b, int l, int h) noexcept : body{b}, lo{l}, hi{h} {}

    // Evaluate: multiply body over index from lo to hi
    [[nodiscard]] constexpr value_type eval(std::span<value_type const> args) const {
        value_type prod = value_type(1);

        // Create extended args array with space for the index
        std::vector<value_type> ext_args(args.begin(), args.end());

        // Insert placeholder for index at IndexDim
        if (ext_args.size() <= IndexDim) {
            ext_args.resize(IndexDim + 1);
        } else {
            ext_args.insert(ext_args.begin() + IndexDim, value_type(0));
        }

        for (int i = lo; i <= hi; ++i) {
            ext_args[IndexDim] = static_cast<value_type>(i);
            prod *= body.eval(std::span<value_type const>(ext_args));
        }

        return prod;
    }

    // Deprecated: use eval() instead
    [[nodiscard]] [[deprecated("use eval() instead")]]
    constexpr value_type evaluate(std::span<value_type const> args) const {
        return eval(args);
    }

    // Derivative: d/dx[Πf(x,i)] is more complex
    // Using logarithmic derivative: d/dx[Π f_i] = (Π f_i) * Σ (f_i'/f_i)
    // For simplicity, we return a product-rule style expansion
    template<std::size_t Dim>
    [[nodiscard]] constexpr auto derivative() const {
        // This is the general derivative for products: Σⱼ (∏ᵢ≠ⱼ fᵢ) * fⱼ'
        // Implemented as: (Π fᵢ) * Σⱼ (fⱼ'/fⱼ)
        constexpr std::size_t AdjustedDim = (Dim >= IndexDim) ? (Dim + 1) : Dim;
        auto dbody = body.template derivative<AdjustedDim>();

        // d(log(f))/dx = f'/f, so d(Πf)/dx = Πf * Σ(f'/f)
        // Return a structure that computes this
        auto term = dbody / body;
        auto sum_of_terms = FiniteSum<decltype(term), IndexDim>{term, lo, hi};
        auto self_copy = *this;
        return self_copy * sum_of_terms;
    }

    // String representation
    [[nodiscard]] std::string to_string() const {
        std::ostringstream oss;
        oss << "(prod[i" << IndexDim << "=" << lo << ".." << hi << "] " << body.to_string() << ")";
        return oss.str();
    }
};

// =============================================================================
// Type traits
// =============================================================================

template<typename T>
struct is_finite_sum : std::false_type {};

template<typename E, std::size_t I>
struct is_finite_sum<FiniteSum<E, I>> : std::true_type {};

template<typename T>
inline constexpr bool is_finite_sum_v = is_finite_sum<T>::value;

template<typename T>
struct is_finite_product : std::false_type {};

template<typename E, std::size_t I>
struct is_finite_product<FiniteProduct<E, I>> : std::true_type {};

template<typename T>
inline constexpr bool is_finite_product_v = is_finite_product<T>::value;

// =============================================================================
// Factory functions
// =============================================================================

// sum<IndexDim>(body, lo, hi): Σᵢ₌ₗₒʰⁱ body
// IndexDim specifies which variable slot is the summation index
template<std::size_t IndexDim, typename E>
    requires is_expr_node_v<E>
[[nodiscard]] constexpr auto sum(E body, int lo, int hi) {
    return FiniteSum<E, IndexDim>{body, lo, hi};
}

// product<IndexDim>(body, lo, hi): Πᵢ₌ₗₒʰⁱ body
template<std::size_t IndexDim, typename E>
    requires is_expr_node_v<E>
[[nodiscard]] constexpr auto product(E body, int lo, int hi) {
    return FiniteProduct<E, IndexDim>{body, lo, hi};
}

// =============================================================================
// Index variable helper
// =============================================================================

// Index variable: use arg<N> or Var<N, T> as the index in sum/product expressions
// Example: sum<1>(arg<1>, 1, 10) computes Σᵢ₌₁¹⁰ i = 55

} // namespace limes::expr
