#pragma once

#include <array>
#include <vector>
#include <random>
#include <cmath>
#include <algorithm>
#include <numeric>
#include "../concepts/multivariate_concepts.hpp"

namespace calckit::multivariate {

// Hyperrectangle: [a1, b1] × [a2, b2] × ... × [aN, bN]
template <typename T, int Dim>
class hyperrectangle {
public:
    using value_type = T;
    static constexpr int dimension_value = Dim;

    std::array<T, Dim> lower;
    std::array<T, Dim> upper;

    // Constructors
    constexpr hyperrectangle() = default;

    constexpr hyperrectangle(std::array<T, Dim> lower, std::array<T, Dim> upper)
        : lower(lower), upper(upper) {
        for (int i = 0; i < Dim; ++i) {
            if (lower[i] > upper[i]) {
                std::swap(this->lower[i], this->upper[i]);
            }
        }
    }

    // Concept requirements
    constexpr int dimension() const { return Dim; }

    constexpr T volume() const {
        T vol = T(1);
        for (int i = 0; i < Dim; ++i) {
            vol *= (upper[i] - lower[i]);
        }
        return vol;
    }

    constexpr std::array<T, Dim> lower_bounds() const { return lower; }
    constexpr std::array<T, Dim> upper_bounds() const { return upper; }

    constexpr bool contains(const std::array<T, Dim>& point) const {
        for (int i = 0; i < Dim; ++i) {
            if (point[i] < lower[i] || point[i] > upper[i]) {
                return false;
            }
        }
        return true;
    }

    // Sampling support
    template <typename RNG>
    std::array<T, Dim> sample_uniform(RNG& rng) const {
        std::uniform_real_distribution<T> dist(0, 1);
        std::array<T, Dim> point;
        for (int i = 0; i < Dim; ++i) {
            T u = dist(rng);
            point[i] = lower[i] + u * (upper[i] - lower[i]);
        }
        return point;
    }

    // Convenience: unit hypercube [0,1]^Dim
    static constexpr hyperrectangle unit() {
        std::array<T, Dim> zeros{}, ones{};
        for (int i = 0; i < Dim; ++i) {
            zeros[i] = T(0);
            ones[i] = T(1);
        }
        return hyperrectangle(zeros, ones);
    }

    // Center point
    constexpr std::array<T, Dim> center() const {
        std::array<T, Dim> c;
        for (int i = 0; i < Dim; ++i) {
            c[i] = (lower[i] + upper[i]) / T(2);
        }
        return c;
    }

    // Diameter (longest diagonal)
    constexpr T diameter() const {
        T diam_sq = T(0);
        for (int i = 0; i < Dim; ++i) {
            T diff = upper[i] - lower[i];
            diam_sq += diff * diff;
        }
        return std::sqrt(diam_sq);
    }

    // Subdivide along longest dimension
    std::pair<hyperrectangle, hyperrectangle> subdivide() const {
        // Find longest dimension
        int max_dim = 0;
        T max_length = upper[0] - lower[0];
        for (int i = 1; i < Dim; ++i) {
            T length = upper[i] - lower[i];
            if (length > max_length) {
                max_length = length;
                max_dim = i;
            }
        }

        // Split at midpoint of longest dimension
        T midpoint = (lower[max_dim] + upper[max_dim]) / T(2);

        hyperrectangle left = *this;
        hyperrectangle right = *this;
        left.upper[max_dim] = midpoint;
        right.lower[max_dim] = midpoint;

        return {left, right};
    }

    // Subdivide along specific dimension
    std::pair<hyperrectangle, hyperrectangle> subdivide_dim(int dim) const {
        T midpoint = (lower[dim] + upper[dim]) / T(2);

        hyperrectangle left = *this;
        hyperrectangle right = *this;
        left.upper[dim] = midpoint;
        right.lower[dim] = midpoint;

        return {left, right};
    }
};

// Simplex: convex hull of Dim+1 vertices
template <typename T, int Dim>
class simplex {
public:
    using value_type = T;
    static constexpr int dimension_value = Dim;

    std::array<std::array<T, Dim>, Dim + 1> vertices;

    // Constructors
    constexpr simplex() = default;

    constexpr simplex(std::array<std::array<T, Dim>, Dim + 1> vertices)
        : vertices(vertices) {}

    // Concept requirements
    constexpr int dimension() const { return Dim; }

    T volume() const {
        // Use determinant formula for simplex volume
        // Volume = |det(v1-v0, v2-v0, ..., vN-v0)| / N!
        if constexpr (Dim == 1) {
            return std::abs(vertices[1][0] - vertices[0][0]);
        } else if constexpr (Dim == 2) {
            // Triangle area
            T x0 = vertices[0][0], y0 = vertices[0][1];
            T x1 = vertices[1][0], y1 = vertices[1][1];
            T x2 = vertices[2][0], y2 = vertices[2][1];
            return std::abs((x1 - x0) * (y2 - y0) - (x2 - x0) * (y1 - y0)) / T(2);
        } else if constexpr (Dim == 3) {
            // Tetrahedron volume
            auto v1 = subtract(vertices[1], vertices[0]);
            auto v2 = subtract(vertices[2], vertices[0]);
            auto v3 = subtract(vertices[3], vertices[0]);
            return std::abs(dot(v1, cross(v2, v3))) / T(6);
        } else {
            // General case: use recursive determinant (not implemented)
            // For now, return bounding box volume as approximation
            return bounding_box().volume();
        }
    }

    hyperrectangle<T, Dim> bounding_box() const {
        std::array<T, Dim> lower, upper;
        for (int d = 0; d < Dim; ++d) {
            lower[d] = vertices[0][d];
            upper[d] = vertices[0][d];
            for (int v = 1; v <= Dim; ++v) {
                lower[d] = std::min(lower[d], vertices[v][d]);
                upper[d] = std::max(upper[d], vertices[v][d]);
            }
        }
        return hyperrectangle<T, Dim>(lower, upper);
    }

    std::array<T, Dim> lower_bounds() const { return bounding_box().lower; }
    std::array<T, Dim> upper_bounds() const { return bounding_box().upper; }

    bool contains(const std::array<T, Dim>& point) const {
        // Point in simplex using barycentric coordinates
        // All barycentric coords >= 0 and sum to 1
        auto bary = barycentric_coordinates(point);
        for (auto lambda : bary) {
            if (lambda < T(0) || lambda > T(1)) return false;
        }
        return true;
    }

    // Sampling support (uniform in simplex)
    template <typename RNG>
    std::array<T, Dim> sample_uniform(RNG& rng) const {
        // Use Dirichlet distribution sampling
        std::gamma_distribution<T> gamma(1.0, 1.0);
        std::array<T, Dim + 1> weights;
        T sum = T(0);
        for (int i = 0; i <= Dim; ++i) {
            weights[i] = gamma(rng);
            sum += weights[i];
        }

        // Normalize and combine vertices
        std::array<T, Dim> point{};
        for (int i = 0; i <= Dim; ++i) {
            T w = weights[i] / sum;
            for (int d = 0; d < Dim; ++d) {
                point[d] += w * vertices[i][d];
            }
        }
        return point;
    }

    // Unit simplex: vertices at origin and unit axes
    static constexpr simplex unit() {
        std::array<std::array<T, Dim>, Dim + 1> verts{};
        // First vertex at origin
        for (int d = 0; d < Dim; ++d) {
            verts[0][d] = T(0);
        }
        // Remaining vertices at unit axes
        for (int v = 1; v <= Dim; ++v) {
            for (int d = 0; d < Dim; ++d) {
                verts[v][d] = (d == v - 1) ? T(1) : T(0);
            }
        }
        return simplex(verts);
    }

private:
    std::array<T, Dim> subtract(const std::array<T, Dim>& a, const std::array<T, Dim>& b) const {
        std::array<T, Dim> result;
        for (int i = 0; i < Dim; ++i) {
            result[i] = a[i] - b[i];
        }
        return result;
    }

    T dot(const std::array<T, Dim>& a, const std::array<T, Dim>& b) const {
        T sum = T(0);
        for (int i = 0; i < Dim; ++i) {
            sum += a[i] * b[i];
        }
        return sum;
    }

    std::array<T, 3> cross(const std::array<T, 3>& a, const std::array<T, 3>& b) const
        requires (Dim == 3) {
        return {
            a[1] * b[2] - a[2] * b[1],
            a[2] * b[0] - a[0] * b[2],
            a[0] * b[1] - a[1] * b[0]
        };
    }

    std::array<T, Dim + 1> barycentric_coordinates(const std::array<T, Dim>& point) const {
        // Solve linear system for barycentric coordinates
        // Simplified for now - return approximate
        std::array<T, Dim + 1> bary;
        bary.fill(T(1) / T(Dim + 1));  // Placeholder
        return bary;
    }
};

// Ball: {x : ||x - center|| <= radius}
template <typename T, int Dim>
class ball {
public:
    using value_type = T;
    static constexpr int dimension_value = Dim;

    std::array<T, Dim> center_point;
    T radius_value;

    // Constructors
    constexpr ball() : radius_value(T(1)) {
        center_point.fill(T(0));
    }

    constexpr ball(std::array<T, Dim> center, T radius)
        : center_point(center), radius_value(radius) {}

    // Concept requirements
    constexpr int dimension() const { return Dim; }

    T volume() const {
        // Volume of unit ball times radius^Dim
        T r_pow_d = std::pow(radius_value, Dim);

        if constexpr (Dim == 1) {
            return T(2) * radius_value;  // Line segment
        } else if constexpr (Dim == 2) {
            return std::numbers::pi_v<T> * radius_value * radius_value;
        } else if constexpr (Dim == 3) {
            return T(4.0 / 3.0) * std::numbers::pi_v<T> * r_pow_d;
        } else {
            // General formula: V_n(r) = π^(n/2) / Γ(n/2 + 1) * r^n
            // Approximate for now
            return std::pow(std::numbers::pi_v<T>, Dim / 2.0) * r_pow_d /
                   std::tgamma(Dim / 2.0 + 1);
        }
    }

    hyperrectangle<T, Dim> bounding_box() const {
        std::array<T, Dim> lower, upper;
        for (int i = 0; i < Dim; ++i) {
            lower[i] = center_point[i] - radius_value;
            upper[i] = center_point[i] + radius_value;
        }
        return hyperrectangle<T, Dim>(lower, upper);
    }

    std::array<T, Dim> lower_bounds() const { return bounding_box().lower; }
    std::array<T, Dim> upper_bounds() const { return bounding_box().upper; }

    bool contains(const std::array<T, Dim>& point) const {
        T dist_sq = T(0);
        for (int i = 0; i < Dim; ++i) {
            T diff = point[i] - center_point[i];
            dist_sq += diff * diff;
        }
        return dist_sq <= radius_value * radius_value;
    }

    // Sampling support (rejection sampling for now)
    template <typename RNG>
    std::array<T, Dim> sample_uniform(RNG& rng) const {
        // Use rejection sampling from bounding box
        auto bbox = bounding_box();
        std::array<T, Dim> point;
        do {
            point = bbox.sample_uniform(rng);
        } while (!contains(point));
        return point;
    }

    // Unit ball at origin
    static constexpr ball unit() {
        std::array<T, Dim> origin{};
        origin.fill(T(0));
        return ball(origin, T(1));
    }
};

// Arbitrary region (using characteristic function)
template <typename T, int Dim, typename CharFunc>
class arbitrary_region {
public:
    using value_type = T;
    static constexpr int dimension_value = Dim;

    hyperrectangle<T, Dim> bbox;
    CharFunc characteristic;  // Returns true if point in region

    // Constructor
    arbitrary_region(hyperrectangle<T, Dim> bounding_box, CharFunc char_func)
        : bbox(bounding_box), characteristic(char_func) {}

    // Concept requirements
    constexpr int dimension() const { return Dim; }

    T volume() const {
        // Estimate via Monte Carlo (could cache this)
        // For now, return bounding box volume (overestimate)
        return bbox.volume();
    }

    std::array<T, Dim> lower_bounds() const { return bbox.lower; }
    std::array<T, Dim> upper_bounds() const { return bbox.upper; }

    bool contains(const std::array<T, Dim>& point) const {
        return bbox.contains(point) && characteristic(point);
    }

    // Sampling support (rejection sampling)
    template <typename RNG>
    std::array<T, Dim> sample_uniform(RNG& rng) const {
        std::array<T, Dim> point;
        do {
            point = bbox.sample_uniform(rng);
        } while (!characteristic(point));
        return point;
    }
};

// Type aliases for common dimensions
template <typename T>
using rectangle = hyperrectangle<T, 2>;

template <typename T>
using box = hyperrectangle<T, 3>;

template <typename T>
using triangle = simplex<T, 2>;

template <typename T>
using tetrahedron = simplex<T, 3>;

template <typename T>
using circle = ball<T, 2>;

template <typename T>
using sphere = ball<T, 3>;

} // namespace calckit::multivariate
