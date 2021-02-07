
namespace math::integration
{
    template <typename Accumulator>
    struct simpson_univariate_integrator
    {
		using value_type = typename Accumulator::value_type;
		using accumulator_type = Accumulator;

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
			return operator()(std::forward<F>(f), t0, t1, default_step_size());
		}

		template <typename F>
        result_type operator()(F&& f, value_type t0, value_type t1, value_type h) const
        {
            int i = 0;
            Accumulator yn, tn(std::move(t0));
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