#pragma once

#include <limits>
#include <cmath>

namespace math::integration
{
    // N models a numerical integration method
	template <typename N>
	struct numerical_univariate_integrator
	{
        using value_type = typename N::value_type;
        using result_type = typename N::result_type;
        using integrator_type = N;

		explicit numerical_univariate_integrator(
            N nint,
            value_type h) :
                nint(std::move(nint)),
				h(h) {}

		numerical_univariate_integrator() :
			nint(),
			h(nint.default_step_size()) {}

		numerical_univariate_integrator(
			numerical_univariate_integrator const & ) = default;

		auto step_size() const { return h; }
		auto default_step_size() const { return nint.default_step_size(); }
		auto integrator() const { return nint; }

		template <typename F>
		auto left_bounded(F&& f, value_type t0) const
		{
			return left_bounded(std::forward<F>(f), t0, _h);
		}

		template <typename F>
		auto left_bounded(F&& f, value_type t0, value_type h) const
		{
			return nint([f = std::forward<F>(f), t0](value_type t) -> value_type
                {
				    const value_type c = value_type(1) / (value_type(1) - t);
				    return f(t0 + c * t) * c * c;
			    }, value_type(0),
                value_type(1) - std::numeric_limits<value_type>::epsilon(),
			    std::move(h));
		}

		template <typename F>
		auto right_bounded(F&& f, value_type t1) const
		{
			return right_bounded(std::forward<F>(f), t1, h);
		}

		template <typename F>
		auto right_bounded(F&& f, value_type t1, value_type h) const
		{
			return nint([f = std::forward<F>(f), t1](value_type t) -> value_type
                {
                    const value_type c = value_type(1) / t;
                    return f(t1 + value_type(1) - c) * c * c;
                },
                value_type(0) + std::numeric_limits<value_type>::epsilon(),
                value_type(1),
                std::move(h));
		}

		template <typename F>
		auto unbounded(F&& f) const
		{
			return unbounded(std::forward<F>(f), _h);
		}

		template <typename F>
		auto unbounded(F&& f, value_type h) const
		{
			const static auto HALF_PI = value_type(1.5707963267949);
			return nint([f = std::forward<F>(f)](value_type t) -> value_type
                {
                    const value_type c = value_type(1) / std::cos(t);
                    return f(std::sin(t) * c) * c * c;
                }, -HALF_PI, HALF_PI, std::move(h));
		}

		template <typename F>
		auto bounded(F&& f, value_type t0, value_type t1) const
		{
			return bounded(std::forward<F>(f), t0, t1, h);
		}

		template <typename F>
		auto bounded(F&& f, value_type t0, value_type t1, value_type h) const
		{
			return operator()(std::forward<F>(f), t0, t1, h);
		}

		template <typename F>
		auto operator()(F&& f, value_type t0, value_type t1, value_type h) const
		{
			return nint(std::forward<F>(f), t0, t1, h);			
		}

		template <typename F>
		auto operator()(F&& f, value_type t0, value_type t1) const
		{
			return nint(std::forward<F>(f), t0, t1, _h);			
		}
		
		N nint;
        value_type h;
	};
}
