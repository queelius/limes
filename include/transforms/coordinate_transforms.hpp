#pragma once

#include <cmath>
#include <numbers>
#include <limits>
#include <memory>
#include <stdexcept>
#include "../concepts/integrator_concepts.hpp"

namespace algebraic_integrators::transforms {

// Base transform interface
template<concepts::Field T>
struct coordinate_transform {
    using value_type = T;

    virtual ~coordinate_transform() = default;
    virtual T forward(T x) const = 0;
    virtual T jacobian(T x) const = 0;
    virtual T inverse(T y) const = 0;
};

// Identity transform
template<concepts::Field T>
struct identity_transform : coordinate_transform<T> {
    constexpr T forward(T x) const override { return x; }
    constexpr T jacobian(T) const override { return T(1); }
    constexpr T inverse(T y) const override { return y; }
};

// Linear transform: y = ax + b
template<concepts::Field T>
class linear_transform : coordinate_transform<T> {
public:
    constexpr linear_transform(T scale, T shift = T(0))
        : scale_{scale}, shift_{shift} {}

    constexpr T forward(T x) const override {
        return scale_ * x + shift_;
    }

    constexpr T jacobian(T) const override {
        return scale_;
    }

    constexpr T inverse(T y) const override {
        return (y - shift_) / scale_;
    }

    constexpr T operator()(T x) const { return forward(x); }

private:
    T scale_;
    T shift_;
};

// Maps [a, b] to [-1, 1]
template<concepts::Field T>
class interval_map : coordinate_transform<T> {
public:
    constexpr interval_map(T a, T b)
        : a_{a}, b_{b}, center_{(a + b) / T(2)}, half_width_{(b - a) / T(2)} {}

    constexpr T forward(T t) const override {
        return center_ + half_width_ * t;
    }

    constexpr T jacobian(T) const override {
        return half_width_;
    }

    constexpr T inverse(T x) const override {
        return (x - center_) / half_width_;
    }

    constexpr T operator()(T t) const { return forward(t); }

private:
    T a_, b_;
    T center_;
    T half_width_;
};

// Sigmoid transform for infinite intervals
template<concepts::Field T>
class sigmoid_transform : coordinate_transform<T> {
public:
    constexpr T forward(T t) const override {
        return t / (T(1) + std::abs(t));
    }

    constexpr T jacobian(T t) const override {
        T denom = T(1) + std::abs(t);
        return T(1) / (denom * denom);
    }

    constexpr T inverse(T y) const override {
        return y / (T(1) - std::abs(y));
    }

    constexpr T operator()(T t) const { return forward(t); }
};

// Tanh transform: maps (-∞, ∞) to (-1, 1)
template<concepts::Field T>
class tanh_transform : coordinate_transform<T> {
public:
    constexpr tanh_transform(T scale = T(1)) : scale_{scale} {}

    constexpr T forward(T t) const override {
        return std::tanh(scale_ * t);
    }

    constexpr T jacobian(T t) const override {
        T sech = T(1) / std::cosh(scale_ * t);
        return scale_ * sech * sech;
    }

    constexpr T inverse(T y) const override {
        return std::atanh(y) / scale_;
    }

    constexpr T operator()(T t) const { return forward(t); }

private:
    T scale_;
};

// Exponential transform for (0, ∞) to (-∞, ∞)
template<concepts::Field T>
class log_transform : coordinate_transform<T> {
public:
    constexpr T forward(T t) const override {
        return std::log(t);
    }

    constexpr T jacobian(T t) const override {
        return T(1) / t;
    }

    constexpr T inverse(T y) const override {
        return std::exp(y);
    }

    constexpr T operator()(T t) const { return forward(t); }
};

// Sinh transform for heavy-tailed integrands
template<concepts::Field T>
class sinh_transform : coordinate_transform<T> {
public:
    constexpr sinh_transform(T scale = T(1)) : scale_{scale} {}

    constexpr T forward(T t) const override {
        return std::sinh(scale_ * t) / scale_;
    }

    constexpr T jacobian(T t) const override {
        return std::cosh(scale_ * t);
    }

    constexpr T inverse(T y) const override {
        return std::asinh(scale_ * y) / scale_;
    }

    constexpr T operator()(T t) const { return forward(t); }

private:
    T scale_;
};

// IMT (Iri-Mori-Takahasi) transform for endpoint singularities
template<concepts::Field T>
class imt_transform : coordinate_transform<T> {
public:
    constexpr imt_transform(T alpha = T(1)) : alpha_{alpha} {}

    constexpr T forward(T t) const override {
        if (std::abs(t) >= T(1)) return t;

        T t2 = t * t;
        T factor = std::pow(T(1) - t2, alpha_);
        return t * factor;
    }

    constexpr T jacobian(T t) const override {
        if (std::abs(t) >= T(1)) return T(1);

        T t2 = t * t;
        T one_minus_t2 = T(1) - t2;
        T factor = std::pow(one_minus_t2, alpha_ - T(1));

        return factor * (one_minus_t2 - T(2) * alpha_ * t2);
    }

    constexpr T inverse(T) const override {
        // IMT inverse is complex to compute analytically
        // In practice, we'd use numerical root finding
        throw std::runtime_error("IMT inverse not implemented");
    }

    constexpr T operator()(T t) const { return forward(t); }

private:
    T alpha_;
};

// Double exponential (tanh-sinh) transform
template<concepts::Field T>
class double_exponential_transform : coordinate_transform<T> {
public:
    constexpr T forward(T t) const override {
        T exp_t = std::exp(t);
        T exp_neg_exp_t = std::exp(-exp_t);
        T u = std::exp(t - exp_t);
        return (u - T(1)/u) / (u + T(1)/u);
    }

    constexpr T jacobian(T t) const override {
        T exp_t = std::exp(t);
        T exp_neg_exp_t = std::exp(-exp_t);
        T u = std::exp(t - exp_t);
        T cosh_val = (u + T(1)/u) / T(2);

        return (T(1) - exp_neg_exp_t) / (cosh_val * cosh_val);
    }

    constexpr T inverse(T y) const override {
        // Complex to compute; typically not needed
        throw std::runtime_error("DE inverse not implemented");
    }

    constexpr T operator()(T t) const { return forward(t); }
};

// Composite transform
template<concepts::Field T>
class composite_transform : coordinate_transform<T> {
public:
    template<typename Transform1, typename Transform2>
    composite_transform(Transform1 t1, Transform2 t2)
        : first_{std::make_unique<transform_impl<Transform1>>(std::move(t1))},
          second_{std::make_unique<transform_impl<Transform2>>(std::move(t2))} {}

    T forward(T x) const override {
        return second_->forward(first_->forward(x));
    }

    T jacobian(T x) const override {
        T y = first_->forward(x);
        return first_->jacobian(x) * second_->jacobian(y);
    }

    T inverse(T z) const override {
        return first_->inverse(second_->inverse(z));
    }

    T operator()(T x) const { return forward(x); }

private:
    struct transform_base {
        virtual ~transform_base() = default;
        virtual T forward(T x) const = 0;
        virtual T jacobian(T x) const = 0;
        virtual T inverse(T y) const = 0;
    };

    template<typename Transform>
    struct transform_impl : transform_base {
        Transform t;

        explicit transform_impl(Transform transform) : t{std::move(transform)} {}

        T forward(T x) const override { return t.forward(x); }
        T jacobian(T x) const override { return t.jacobian(x); }
        T inverse(T y) const override { return t.inverse(y); }
    };

    std::unique_ptr<transform_base> first_;
    std::unique_ptr<transform_base> second_;
};

// Transform factories
template<concepts::Field T>
constexpr auto make_identity() { return identity_transform<T>{}; }

template<concepts::Field T>
constexpr auto make_linear(T scale, T shift = T(0)) {
    return linear_transform<T>{scale, shift};
}

template<concepts::Field T>
constexpr auto make_interval_map(T a, T b) {
    return interval_map<T>{a, b};
}

template<concepts::Field T>
constexpr auto make_sigmoid() { return sigmoid_transform<T>{}; }

template<concepts::Field T>
constexpr auto make_tanh(T scale = T(1)) {
    return tanh_transform<T>{scale};
}

template<concepts::Field T>
constexpr auto make_log() { return log_transform<T>{}; }

template<concepts::Field T>
constexpr auto make_sinh(T scale = T(1)) {
    return sinh_transform<T>{scale};
}

template<concepts::Field T>
constexpr auto make_imt(T alpha = T(1)) {
    return imt_transform<T>{alpha};
}

template<concepts::Field T>
constexpr auto make_double_exponential() {
    return double_exponential_transform<T>{};
}

// Compose transforms
template<typename T1, typename T2>
auto compose(T1&& t1, T2&& t2) {
    using T = typename std::decay_t<T1>::value_type;
    return composite_transform<T>{std::forward<T1>(t1), std::forward<T2>(t2)};
}

} // namespace algebraic_integrators::transforms