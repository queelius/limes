#pragma once

namespace alex::math
{
	template <typename Accumulator>
	struct gauss_quad
	{
		using value_type = typename Accumulator::value_type;

		struct result_type
		{
			value_type operator()() { return sum; };

			value_type sum;
			value_type error_estimate;
			int iterations;
		};

		value_type default_step_size() const { return 1e-3; };

		template <typename F>
        result_type operator()(F&& f, value_type t0, value_type t1) const
		{
			return operator()(std::forward<F>(f), t0, t1, default_step_size());
		};

		template <typename F>
		result_type operator()(F&& f, value_type t0, value_type t1, value_type h) const
		{
			constexpr value_type w[] =
			{
				0.2025782419255613, 0.1984314853271116,
				0.1984314853271116, 0.1861610000155622,
				0.1861610000155622,	0.1662692058169939,
				0.1662692058169939,	0.1395706779261543,
				0.1395706779261543,	0.1071592204671719,
				0.1071592204671719,	0.0703660474881081,
				0.0703660474881081,	0.0307532419961173,
				0.0307532419961173
			};

			constexpr value_type x[] =
			{
				0.0000000000000000, -0.2011940939974345,
				0.2011940939974345, -0.3941513470775634,
				0.3941513470775634, -0.5709721726085388,
				0.5709721726085388, -0.7244177313601701,
				0.7244177313601701, -0.8482065834104272,
				0.8482065834104272, -0.9372733924007060,
				0.9372733924007060, -0.9879925180204854,
				0.9879925180204854
			};

			auto i = 0;
			auto tn = std::move(t0);
			value_type yn = 0;
			bool not_done = true;
			while (not_done)
			{
				++i;
				if (tn + h >= t1)
				{
					h = tn - t1;
					not_done = false;
				}

				value_type
					c = h / value_type(2),
					d = (value_type(2) * tn + h) / value_type(2),
					s = 0;

				for (auto i = 0; i < 15; ++i)
					s += w[i] * f(c * x[i] + d);

				yn += c * s;
				tn += h;
			}

			return result_type{yn, h * h * h * h, i};
		};
	};
}