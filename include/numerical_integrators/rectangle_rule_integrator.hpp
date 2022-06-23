namespace math::integration
{
	template <typename A>
    struct rectangle_rule_integrator
    {
		using value_type = typename A::value_type;

        struct result_type
        {
			auto operator()() const { return sum; }
            operator value_type() const { return sum; }

            value_type sum;
            int iterations;
            value_type error_estimate;
        };

		auto default_step_size() const { return value_type(1e-3); }

		template <typename F>
        auto operator()(F f, value_type t0,
            value_type t1, value_type h = default_step_size()) const
        {
			A y = T(0);
            int j = 0;
            while (t0 < t1 - h)
            {
                ++j;
                y += f(t0) * h;
                t0 += h;
            }
			y += f(t1 - t0) * (t1 - t0);
            return result_type{y(),j,h};
       	}
    };
}