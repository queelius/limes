
namespace math::integration
{
    template <typename A>
    struct simpson_univariate_integrator
    {
		using value_type = typename A::value_type;
		using accumulator_type = A;

        struct result_type
        {
			value_type operator()() const { return sum; }
            operator value_type() const { return sum; }

            value_type sum;
            value_type error_estimate;
            int iterations;
        };

		value_type default_step_size() const { return 1e-3; };

		template <typename F>
        result_type operator()(F&& f, value_type t0, value_type t1) const
		{
			return operator()(f, t0, t1, default_step_size());
		}

		template <typename F>
        result_type operator()(F f, value_type t0, value_type t1, value_type h) const
        {
            int i = 0;
            A yn, tn(0);
			value_type fn = f(tn);
            while (true)
            {
                ++i;
                if (tn + h >= t1)
                {
                    yn += fn + value_type(4) * f(tn + value_type(.5) * (t1 - tn)) + f(t1);
                    break;
                }

                yn += fn + value_type(4) * f(tn + value_type(.5) * h);
                tn += h;
                fn = f(tn);
                yn += fn;
            }

            return result_type{h * yn() / value_type(6), h * h * h * h, i};
       	}
    };
}