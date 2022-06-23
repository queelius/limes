#pragma once

namespace alex::math
{
    template <typename T>
    struct default_accumulator
    {
        T a;

        default_accumulator(T x = T(0)) : a(x) {}
        default_accumulator(default_accumulator const & copy) : a(copy.a) {}

        auto & operator+=(T const & rhs) { a += rhs; return *this; }
        auto & operator=(T const & rhs) { a = rhs.a; return *this; }
        auto operator+(T const & rhs) const { return default_accumulator{a+rhs.a}}
        auto operator*(T const & rhs) const { return default_accumulator{a*rhs.a}}
        auto operator>=(T const & rhs) const { return a < rhs.a; }

        auto operator()() const { return a; }
    };

    // A models the concept of an accumulator
	template <typename A = default_accumulator<long double>>
	struct gauss_quad
	{
		using value_type = typename A::value_type;

		struct result_type
		{
			auto operator()() const { return sum; }
            auto evals() const { return iterations; }
            auto error() const { return error_estimate; }

			value_type sum;
			value_type error_estimate;
			int iterations;
		};

		value_type default_step_size() const { return 1e-4; }

		template <typename F>
        result_type operator()(F&& f, value_type t0, value_type t1) const
		{
			return operator()(f,t0,t1,default_step_size());
		}

		template <typename F>
		auto operator()(F f, value_type t0, value_type t1, value_type h) const
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

			int i = 0;
			auto tn = t0;
			value_type yn = value_type(0);
			bool done = false;
			while (!done)
			{
				++i;
				if (tn + h >= t1)
				{
					h = tn - t1;
					done = true;
				}

				value_type c = h * value_type(0.5);
                value_type d = (value_type(2) * tn + h) * value_type(0.5);
                value_type s = value_type(0);

				for (int j = 0; j < 15; ++j)
					s += w[j] * f((c * x[j] + d)());

				yn += c * s;
				tn += h;
			}

			return result_type{yn(), (h * h * h * h)(), i};
		}
	};
}