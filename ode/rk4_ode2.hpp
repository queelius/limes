#pragma once

#include <fstream>

namespace alex::math
{
	template <typename Accumulator>
	struct rk4_ode2
	{
		using value_type = typename Accumulator::value_type;

		struct result_type
		{
			value_type operator()() const { return y; };
			value_type y;
            value_type m;
			value_type error_estimate;
			int iterations;
		};

		value_type default_step_size() const { return 1e-3; };

        template <typename F>
        result_type linear_bvp(F&& f, value_type t0, value_type y0, value_type t1, value_type y1, value_type lb, value_type ub, value_type h) const
        {
            // implemented via the interpolation method
            auto s = ivp(f, t0, y0, lb, t1, h);
            auto q0 = s.y;
            int i = s.iterations;

            s = ivp(f, t0, y0, ub, t1, h);
            auto q1 = s.y;
            i += s.iterations;

            auto m = lb + (ub - lb) / (q1 - q0) * (y1 - q0);

            s = ivp(f, std::move(t0), std::move(y0), m, std::move(t1), std::move(h));
            i += s.iterations;

            return result_type{std::move(s.y), std::move(m), std::move(s.error_estimate), i};
        };

		template <typename F>
        result_type ivp(F&& f, value_type t0, value_type y0, value_type m0, value_type t, value_type h) const
        {
            Accumulator yn(y0), tn(t0), mn(m0);
			value_type k1, k2, k3, k4;

            std::ofstream out("graph.dat");

            int i = 0;
			while (value_type(tn) < t)
			{
				++i;
				if (tn + h > t)
					h = t - tn;
                auto y_old = value_type(yn);
                yn += mn * h;

				k1 = f(tn, y_old, mn);
				k2 = f(tn + value_type(.5) * h, y_old + value_type(.5) * h * k1, mn);
				k3 = f(tn + value_type(.5) * h, y_old + value_type(.5) * h * k2, mn);
				k4 = f(tn + h, std::move(y_old) + h * k3, mn);
				mn += (h * (std::move(k1) + value_type(2) * std::move(k2) +
					value_type(2) * std::move(k3) + std::move(k4)) / value_type(6));

                out << tn << "\t" << yn << "\n";

				tn += h;
			}
			return result_type{value_type(yn), value_type(mn), h * h * h * h, i};
       	};
	};
}