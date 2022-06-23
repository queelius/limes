#pragma once

namespace alex::math
{
	template <typename Accumulator>
	struct rk4_ode1
	{
		using value_type = typename Accumulator::value_type;

		struct result_type
		{
			value_type operator()() const { return value; };
			value_type value;
			value_type error_estimate;
			int iterations;
		};

		value_type default_step_size() const { return 1e-3; };

		template <typename F>
        result_type operator()(F&& f, value_type t0, value_type y0, value_type t) const
        {
            return operator()(std::forward<F>(f), t0, y0, t, default_step_size());
       	};

		template <typename F>
		result_type operator()(F&& f, value_type t0, value_type y0, value_type t, value_type h)
		{
			int i = 0;
			Accumulator yn(y0), tn(t0);
			value_type k1, k2, k3, k4;
			
			while (value_type(tn) < t)
			{
				++i;
				if (tn + h > t)
					h = t - tn;

				k1 = f(tn, yn);
				k2 = f(tn + value_type(.5) * h, yn + value_type(.5) * h * k1);
				k3 = f(tn + value_type(.5) * h, yn + value_type(.5) * h * k2);
				k4 = f(tn + h, yn + h * k3);
				yn += (h * (std::move(k1) + value_type(2) * std::move(k2) +
					value_type(2) * std::move(k3) + std::move(k4)) / value_type(6));
				tn += h;
			}

			return result_type{value_type(yn), h * h * h * h, i};
		};
	};
}