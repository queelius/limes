#pragma once

#include <array>
#include <cmath>
#include <span>
#include <numbers>
#include <ranges>
#include "../concepts/integrator_concepts.hpp"

namespace calckit::quadrature {

// Base quadrature rule interface
template<concepts::Field T, std::size_t N>
struct quadrature_rule {
    using value_type = T;
    using size_type = std::size_t;

    static constexpr size_type size() noexcept { return N; }

    std::array<T, N> weights;
    std::array<T, N> abscissas;

    constexpr T weight(size_type i) const noexcept { return weights[i]; }
    constexpr T abscissa(size_type i) const noexcept { return abscissas[i]; }

    // Note: std::views::zip requires C++23, so we provide a simpler interface
    struct node_view {
        const T& x;
        const T& w;
    };

    constexpr auto node(size_type i) const noexcept {
        return node_view{abscissas[i], weights[i]};
    }
};

// Gauss-Legendre quadrature
template<concepts::Field T, std::size_t N>
struct gauss_legendre : quadrature_rule<T, N> {
    constexpr gauss_legendre() noexcept {
        compute_nodes();
    }

private:
    constexpr void compute_nodes() noexcept;
};

// Specialized implementations for common orders
template<concepts::Field T>
struct gauss_legendre<T, 2> : quadrature_rule<T, 2> {
    constexpr gauss_legendre() noexcept {
        constexpr T sqrt3 = T(0.5773502691896257645091487805019574556);
        this->abscissas = {-sqrt3, sqrt3};
        this->weights = {T(1), T(1)};
    }
};

template<concepts::Field T>
struct gauss_legendre<T, 3> : quadrature_rule<T, 3> {
    constexpr gauss_legendre() noexcept {
        constexpr T sqrt35 = T(0.7745966692414833770358530799564799221);
        this->abscissas = {-sqrt35, T(0), sqrt35};
        this->weights = {T(5)/T(9), T(8)/T(9), T(5)/T(9)};
    }
};

template<concepts::Field T>
struct gauss_legendre<T, 5> : quadrature_rule<T, 5> {
    constexpr gauss_legendre() noexcept {
        constexpr T x1 = T(0.9061798459386639927976268782993929651);
        constexpr T x2 = T(0.5384693101056830910363144207002088049);
        constexpr T w1 = T(0.2369268850561890875142640407199173626);
        constexpr T w2 = T(0.4786286704993664680412915148356381929);
        constexpr T w3 = T(0.5688888888888888888888888888888888889);

        this->abscissas = {-x1, -x2, T(0), x2, x1};
        this->weights = {w1, w2, w3, w2, w1};
    }
};

// Gauss-Kronrod quadrature (extends Gauss-Legendre)
template<concepts::Field T>
struct gauss_kronrod_15 : quadrature_rule<T, 15> {
    constexpr gauss_kronrod_15() noexcept {
        // Gauss-Kronrod 7-15 rule
        this->abscissas = {
            T(-0.9914553711208126), T(-0.9491079123427585), T(-0.8648644233597691),
            T(-0.7415311855993944), T(-0.5860872354676911), T(-0.4058451513773972),
            T(-0.2077849550078985), T(0),
            T(0.2077849550078985), T(0.4058451513773972), T(0.5860872354676911),
            T(0.7415311855993944), T(0.8648644233597691), T(0.9491079123427585),
            T(0.9914553711208126)
        };

        this->weights = {
            T(0.0229353220105292), T(0.0630920926299785), T(0.1047900103222502),
            T(0.1406532597155259), T(0.1690047266392679), T(0.1903505780647854),
            T(0.2044329400752989), T(0.2094821410847278),
            T(0.2044329400752989), T(0.1903505780647854), T(0.1690047266392679),
            T(0.1406532597155259), T(0.1047900103222502), T(0.0630920926299785),
            T(0.0229353220105292)
        };
    }

    // Embedded Gauss rule for error estimation
    static constexpr std::size_t gauss_size = 7;
    std::array<T, gauss_size> gauss_weights = {
        T(0.1294849661688697), T(0.2797053914892767), T(0.3818300505051189),
        T(0.4179591836734694),
        T(0.3818300505051189), T(0.2797053914892767), T(0.1294849661688697)
    };
    std::array<std::size_t, gauss_size> gauss_indices = {1, 3, 5, 7, 9, 11, 13};
};

// Clenshaw-Curtis quadrature
template<concepts::Field T, std::size_t N>
struct clenshaw_curtis : quadrature_rule<T, N> {
    constexpr clenshaw_curtis() noexcept {
        compute_nodes();
    }

private:
    constexpr void compute_nodes() noexcept {
        constexpr T pi = std::numbers::pi_v<T>;

        for (std::size_t i = 0; i < N; ++i) {
            T theta = pi * T(i) / T(N - 1);
            this->abscissas[i] = -std::cos(theta);

            // Compute weights using FFT-based approach for efficiency
            T w = T(2) / T(N - 1);
            if (i == 0 || i == N - 1) {
                w /= T(2);
            }

            for (std::size_t j = 1; j < (N - 1) / 2; ++j) {
                w -= T(2) * std::cos(T(2) * j * theta) / (T(4) * j * j - T(1));
            }

            this->weights[i] = w;
        }
    }
};

// Simpson's rule nodes
template<concepts::Field T>
struct simpson_rule : quadrature_rule<T, 3> {
    constexpr simpson_rule() noexcept {
        this->abscissas = {T(-1), T(0), T(1)};
        this->weights = {T(1)/T(3), T(4)/T(3), T(1)/T(3)};
    }
};

// Trapezoidal rule nodes
template<concepts::Field T>
struct trapezoidal_rule : quadrature_rule<T, 2> {
    constexpr trapezoidal_rule() noexcept {
        this->abscissas = {T(-1), T(1)};
        this->weights = {T(1), T(1)};
    }
};

// Midpoint rule
template<concepts::Field T>
struct midpoint_rule : quadrature_rule<T, 1> {
    constexpr midpoint_rule() noexcept {
        this->abscissas = {T(0)};
        this->weights = {T(2)};
    }
};

// Tanh-sinh (double exponential) quadrature nodes
template<concepts::Field T>
class tanh_sinh_nodes {
public:
    using value_type = T;

    struct node {
        T abscissa;
        T weight;
    };

    static constexpr std::size_t max_level = 10;

    constexpr T abscissa(std::size_t level, std::size_t index) const noexcept {
        T h = T(1) / std::pow(T(2), level);
        T t = (index == 0) ? T(0) : h * T(index);
        return transform_abscissa(t);
    }

    constexpr T weight(std::size_t level, std::size_t index) const noexcept {
        T h = T(1) / std::pow(T(2), level);
        T t = (index == 0) ? T(0) : h * T(index);
        return h * transform_weight(t);
    }

private:
    static constexpr T transform_abscissa(T t) noexcept {
        T expmt = std::exp(-t);
        T u = std::exp(t - expmt);
        return (u - T(1)/u) / (u + T(1)/u);
    }

    static constexpr T transform_weight(T t) noexcept {
        T expmt = std::exp(-t);
        T u = std::exp(t - expmt);
        T cosh_val = (u + T(1)/u) / T(2);
        return (T(1) + expmt) / (cosh_val * cosh_val);
    }
};

// Factory functions
template<concepts::Field T, std::size_t N>
constexpr auto make_gauss_legendre() { return gauss_legendre<T, N>{}; }

template<concepts::Field T>
constexpr auto make_gauss_kronrod() { return gauss_kronrod_15<T>{}; }

template<concepts::Field T, std::size_t N>
constexpr auto make_clenshaw_curtis() { return clenshaw_curtis<T, N>{}; }

template<concepts::Field T>
constexpr auto make_simpson() { return simpson_rule<T>{}; }

template<concepts::Field T>
constexpr auto make_trapezoidal() { return trapezoidal_rule<T>{}; }

template<concepts::Field T>
constexpr auto make_midpoint() { return midpoint_rule<T>{}; }

} // namespace calckit::quadrature